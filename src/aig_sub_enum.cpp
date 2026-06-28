#include "gene_multi_aig/aag.hpp"
#include "gene_multi_aig/sub_aig.hpp"

#include <chrono>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace gene_multi_aig;

namespace {

struct Options {
    fs::path input;
    fs::path vars;
    fs::path out_dir;
    fs::path tmp_dir;
    std::string aigtoaig = "aigtoaig";
    bool keep_temp = false;
    bool verbose = false;
};

void print_usage(std::ostream& os) {
    os << "Usage: aig_sub_enum --in input.aig --vars vars.txt --out-dir DIR [options]\n"
       << "\n"
       << "Options:\n"
       << "  --aigtoaig PATH       Path to external aigtoaig command [default: aigtoaig]\n"
       << "  --tmp-dir DIR         Directory for temporary AAG files\n"
       << "  --keep-temp           Keep temporary AAG files\n"
       << "  --verbose             Print progress details\n"
       << "  --help                Show this help\n";
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw AigError("Missing value for " + name);
            }
            ++i;
            return argv[i];
        };
        auto option_value = [&](const std::string& name) -> std::optional<std::string> {
            if (arg == name) {
                return require_value(name);
            }
            const std::string prefix = name + "=";
            if (arg.rfind(prefix, 0) == 0U) {
                return arg.substr(prefix.size());
            }
            return std::nullopt;
        };

        if (arg == "--help" || arg == "-h") {
            print_usage(std::cout);
            std::exit(0);
        } else if (const auto value = option_value("--in")) {
            options.input = *value;
        } else if (const auto value = option_value("--vars")) {
            options.vars = *value;
        } else if (const auto value = option_value("--out-dir")) {
            options.out_dir = *value;
        } else if (const auto value = option_value("--aigtoaig")) {
            options.aigtoaig = *value;
        } else if (const auto value = option_value("--tmp-dir")) {
            options.tmp_dir = *value;
        } else if (arg == "--keep-temp") {
            options.keep_temp = true;
        } else if (arg == "--verbose") {
            options.verbose = true;
        } else {
            throw AigError("Unknown option: " + arg);
        }
    }

    if (options.input.empty()) {
        throw AigError("Missing required --in");
    }
    if (options.vars.empty()) {
        throw AigError("Missing required --vars");
    }
    if (options.out_dir.empty()) {
        throw AigError("Missing required --out-dir");
    }
    if (options.tmp_dir.empty()) {
        options.tmp_dir = fs::temp_directory_path();
    }
    return options;
}

std::string shell_quote(const fs::path& path) {
#ifdef _WIN32
    std::string value = path.string();
    std::string out = "\"";
    for (const char ch : value) {
        if (ch == '"') {
            out += "\\\"";
        } else {
            out += ch;
        }
    }
    out += "\"";
    return out;
#else
    std::string value = path.string();
    std::string out = "'";
    for (const char ch : value) {
        if (ch == '\'') {
            out += "'\\''";
        } else {
            out += ch;
        }
    }
    out += "'";
    return out;
#endif
}

std::string shell_quote_command(const std::string& command) {
#ifdef _WIN32
    if (command.find_first_of(" \t\"") == std::string::npos) {
        return command;
    }
    std::string out = "\"";
    for (const char ch : command) {
        if (ch == '"') {
            out += "\\\"";
        } else {
            out += ch;
        }
    }
    out += "\"";
    return out;
#else
    if (command.find('/') != std::string::npos || command.find(' ') != std::string::npos) {
        std::string out = "'";
        for (const char ch : command) {
            if (ch == '\'') {
                out += "'\\''";
            } else {
                out += ch;
            }
        }
        out += "'";
        return out;
    }
    return command;
#endif
}

void run_aigtoaig(const std::string& aigtoaig,
                  const fs::path& input,
                  const fs::path& output,
                  bool verbose) {
    const std::string command = shell_quote_command(aigtoaig) + " " +
                                shell_quote(input) + " " + shell_quote(output);
    if (verbose) {
        std::cerr << "[aig_sub_enum] " << command << '\n';
    }
    const int status = std::system(command.c_str());
    if (status != 0) {
        throw AigError("aigtoaig failed with status " + std::to_string(status) +
                       ": " + command);
    }
}

fs::path make_temp_path(const fs::path& dir, const std::string& suffix) {
    const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::random_device rd;
    std::mt19937_64 rng(static_cast<std::uint64_t>(now) ^ rd());
    for (int i = 0; i < 100; ++i) {
        std::ostringstream name;
        name << "aig_sub_enum_" << std::hex << rng() << suffix;
        fs::path candidate = dir / name.str();
        if (!fs::exists(candidate)) {
            return candidate;
        }
    }
    throw AigError("Failed to allocate a temporary path in " + dir.string());
}

bool has_extension(const fs::path& path, const std::string& extension) {
    std::string actual = path.extension().string();
    for (char& ch : actual) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return actual == extension;
}

void remove_if_exists(const fs::path& path) {
    std::error_code ec;
    fs::remove(path, ec);
}

void write_sub_aigs(const Aag& base,
                    const std::vector<Lit>& variables,
                    const Options& options) {
    if (variables.size() > 20U) {
        throw AigError("Refusing to generate more than 2^20 sub-AIG files");
    }

    fs::create_directories(options.out_dir);

    const std::uint64_t total = 1ULL << variables.size();
    for (std::uint64_t mask = 0; mask < total; ++mask) {
        std::vector<Lit> units;
        units.reserve(variables.size());
        for (std::size_t i = 0; i < variables.size(); ++i) {
            const bool value = (mask & (1ULL << i)) != 0U;
            units.push_back(value ? variables[i] : invert_lit(variables[i]));
        }

        const std::string bits = assignment_bits(mask, variables.size());
        const fs::path sub_aig = sub_aig_path(options.input, options.out_dir, bits);
        const fs::path sub_aag = make_temp_path(options.tmp_dir, "_sub.aag");
        write_aag(constrain_output_with_units(base, units), sub_aag);
        run_aigtoaig(options.aigtoaig, sub_aag, sub_aig, options.verbose);
        if (!options.keep_temp) {
            remove_if_exists(sub_aag);
        } else if (options.verbose) {
            std::cerr << "[aig_sub_enum] kept temp sub AAG: " << sub_aag << '\n';
        }
    }

    std::cerr << "[aig_sub_enum] wrote " << total
              << " sub-AIG files to " << options.out_dir << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        fs::create_directories(options.tmp_dir);

        const fs::path input_aag = has_extension(options.input, ".aag")
                                       ? options.input
                                       : make_temp_path(options.tmp_dir, "_in.aag");

        if (input_aag != options.input) {
            run_aigtoaig(options.aigtoaig, options.input, input_aag, options.verbose);
        }

        Aag base = read_aag(input_aag);
        const std::vector<Lit> variables = read_sub_variables(base, options.vars);
        write_sub_aigs(base, variables, options);

        if (!options.keep_temp && input_aag != options.input) {
            remove_if_exists(input_aag);
        } else if (options.keep_temp && input_aag != options.input) {
            std::cerr << "[aig_sub_enum] kept temp input AAG: " << input_aag << '\n';
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "aig_sub_enum: " << ex.what() << '\n';
        return 1;
    }
}
