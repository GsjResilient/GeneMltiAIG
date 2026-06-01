#include "gene_multi_aig/aag.hpp"
#include "gene_multi_aig/filters.hpp"
#include "gene_multi_aig/validation.hpp"

#include <chrono>
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
    fs::path output;
    fs::path tmp_dir;
    std::string aigtoaig = "aigtoaig";
    fs::path fastlec;
    unsigned fastlec_cores = 8;
    unsigned fastlec_timeout = 100;
    std::optional<Lit> lhs_lit;
    std::optional<Lit> rhs_lit;
    bool keep_temp = false;
    bool verbose = false;
};

struct ExpandedCircuit {
    Aag aag;
    Lit left = 0;
    Lit right = 0;
};

struct MaterializedRootCircuit {
    Lit group1_root = 0;
    Lit group2_root = 0;
};

void print_usage(std::ostream& os) {
    os << "Usage: ec_expand --in input.aig --out expanded.aig [options]\n"
       << "\n"
       << "Options:\n"
       << "  --aigtoaig PATH       Path to external aigtoaig command\n"
       << "  --fastlec PATH        Path to external fastLEC command for required validation\n"
       << "  --fastlec-cores N     CPU threads for fastLEC validation [default: 8]\n"
       << "  --fastlec-timeout S   Timeout seconds for fastLEC validation [default: 100]\n"
       << "  --tmp-dir DIR         Directory for temporary AAG files\n"
       << "  --keep-temp           Keep temporary AAG files and validation logs\n"
       << "  --lhs-lit L           Explicit left miter root literal\n"
       << "  --rhs-lit R           Explicit right miter root literal\n"
       << "  --verbose             Print progress details\n"
       << "  --help                Show this help\n";
}

Lit parse_lit_arg(const std::string& value, const std::string& option) {
    std::size_t consumed = 0;
    unsigned long parsed = 0;
    try {
        parsed = std::stoul(value, &consumed, 10);
    } catch (const std::exception&) {
        throw AigError("Invalid literal for " + option + ": " + value);
    }
    if (consumed != value.size() || parsed > UINT32_MAX) {
        throw AigError("Invalid literal for " + option + ": " + value);
    }
    return static_cast<Lit>(parsed);
}

unsigned parse_unsigned_arg(const std::string& value, const std::string& option) {
    std::size_t consumed = 0;
    unsigned long parsed = 0;
    try {
        parsed = std::stoul(value, &consumed, 10);
    } catch (const std::exception&) {
        throw AigError("Invalid unsigned value for " + option + ": " + value);
    }
    if (consumed != value.size() || parsed > UINT32_MAX) {
        throw AigError("Invalid unsigned value for " + option + ": " + value);
    }
    return static_cast<unsigned>(parsed);
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

        if (arg == "--help" || arg == "-h") {
            print_usage(std::cout);
            std::exit(0);
        } else if (arg == "--in") {
            options.input = require_value(arg);
        } else if (arg == "--out") {
            options.output = require_value(arg);
        } else if (arg == "--aigtoaig") {
            options.aigtoaig = require_value(arg);
        } else if (arg == "--fastlec") {
            options.fastlec = require_value(arg);
        } else if (arg == "--fastlec-cores") {
            options.fastlec_cores = parse_unsigned_arg(require_value(arg), arg);
        } else if (arg == "--fastlec-timeout") {
            options.fastlec_timeout = parse_unsigned_arg(require_value(arg), arg);
        } else if (arg == "--tmp-dir") {
            options.tmp_dir = require_value(arg);
        } else if (arg == "--keep-temp") {
            options.keep_temp = true;
        } else if (arg == "--lhs-lit") {
            options.lhs_lit = parse_lit_arg(require_value(arg), arg);
        } else if (arg == "--rhs-lit") {
            options.rhs_lit = parse_lit_arg(require_value(arg), arg);
        } else if (arg == "--verbose") {
            options.verbose = true;
        } else {
            throw AigError("Unknown option: " + arg);
        }
    }

    if (options.input.empty()) {
        throw AigError("Missing required --in");
    }
    if (options.output.empty()) {
        throw AigError("Missing required --out");
    }
    if (options.fastlec.empty()) {
        throw AigError("Missing required --fastlec for unique equivalence-pair validation");
    }
    if (static_cast<bool>(options.lhs_lit) != static_cast<bool>(options.rhs_lit)) {
        throw AigError("--lhs-lit and --rhs-lit must be provided together");
    }
    if (options.tmp_dir.empty()) {
        options.tmp_dir = fs::temp_directory_path();
    }
    if (options.fastlec_cores == 0U) {
        throw AigError("--fastlec-cores must be greater than zero");
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
        std::cerr << "[ec_expand] " << command << '\n';
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
        name << "ec_expand_" << std::hex << rng() << suffix;
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

void ensure_combinational_single_output(const Aag& graph) {
    if (graph.latch_count != 0U) {
        throw AigError("Sequential AIGs are not supported: latch count is " +
                       std::to_string(graph.latch_count));
    }
    if (graph.outputs.size() != 1U) {
        throw AigError("Expected exactly one miter output, found " +
                       std::to_string(graph.outputs.size()));
    }
}

void ensure_literal_exists(const Aag& graph, Lit lit, const std::string& name) {
    if (lit == 0U || lit == 1U) {
        return;
    }
    if (!graph.has_var(lit_var(lit))) {
        throw AigError(name + " references unknown variable: " + std::to_string(lit));
    }
}

XorMatch resolve_roots(const Aag& graph, const Options& options) {
    if (options.lhs_lit && options.rhs_lit) {
        ensure_literal_exists(graph, *options.lhs_lit, "--lhs-lit");
        ensure_literal_exists(graph, *options.rhs_lit, "--rhs-lit");
        return XorMatch{*options.lhs_lit, *options.rhs_lit};
    }

    const auto match = match_xor(graph, graph.outputs.front());
    if (!match) {
        throw AigError("Could not recognize the single PO as an XOR miter. "
                       "Use --lhs-lit and --rhs-lit to specify the compared roots.");
    }
    return *match;
}

MaterializedRootCircuit materialize_root_circuit(NetworkCopier& copier, Lit root) {
    return MaterializedRootCircuit{copier.copy_lit(root, 1), copier.copy_lit(root, 2)};
}

Lit build_difference_filter(AigBuilder& builder, const MaterializedRootCircuit& circuit) {
    return builder.add_and(circuit.group1_root, invert_lit(circuit.group2_root));
}

ExpandedCircuit expand_circuit(const Aag& source, Lit lhs_root, Lit rhs_root) {
    AigBuilder builder; //创建一个新电路构造器。后面所有新 input、新 gate 都会通过它生成，而不是直接修改原始 source
    std::vector<Lit> group1_inputs;
    std::vector<Lit> group2_inputs;
    group1_inputs.reserve(source.inputs.size());
    group2_inputs.reserve(source.inputs.size());

    for (std::size_t i = 0; i < source.inputs.size(); ++i) {
        group1_inputs.push_back(builder.add_input());
    }
    for (std::size_t i = 0; i < source.inputs.size(); ++i) {
        group2_inputs.push_back(builder.add_input());
    }

    NetworkCopier copier(source, builder, group1_inputs, group2_inputs);
    const MaterializedRootCircuit lhs_circuit = materialize_root_circuit(copier, lhs_root);
    const MaterializedRootCircuit rhs_circuit = materialize_root_circuit(copier, rhs_root);

    const Lit left = build_difference_filter(builder, lhs_circuit);
    const Lit right = build_difference_filter(builder, rhs_circuit);
    const Lit output = builder.add_xor(left, right);

    return ExpandedCircuit{builder.to_aag({output}), left, right};
}

void remove_if_exists(const fs::path& path) {
    std::error_code ec;
    fs::remove(path, ec);
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        fs::create_directories(options.tmp_dir);

        const fs::path input_aag = has_extension(options.input, ".aag")
                                       ? options.input
                                       : make_temp_path(options.tmp_dir, "_in.aag");
        const fs::path output_aag = has_extension(options.output, ".aag")
                                        ? options.output
                                        : make_temp_path(options.tmp_dir, "_out.aag");
        const fs::path fastlec_input_log = options.fastlec.empty()
                                               ? fs::path{}
                                               : make_temp_path(options.tmp_dir, "_fastlec_input.log");
        const fs::path fastlec_output_log = make_temp_path(options.tmp_dir, "_fastlec_output.log");
        const FastLecValidationConfig fastlec_input_config{options.fastlec,
                                                          fastlec_input_log,
                                                          options.fastlec_cores,
                                                          options.fastlec_timeout,
                                                          options.verbose};
        const FastLecValidationConfig fastlec_output_config{options.fastlec,
                                                           fastlec_output_log,
                                                           options.fastlec_cores,
                                                           options.fastlec_timeout,
                                                           options.verbose};

        if (input_aag != options.input) {
            run_aigtoaig(options.aigtoaig, options.input, input_aag, options.verbose);
        }

        Aag source = read_aag(input_aag);
        ensure_combinational_single_output(source);

        const XorMatch roots = resolve_roots(source, options);//找到a和b电路的最终输出。
        if (options.verbose) {
            std::cerr << "[ec_expand] miter roots: " << roots.lhs
                      << " XOR " << roots.rhs << '\n';
        }

        const ValidationResult input_validation =
            validate_input_single_eq_pair(
                options.input,
                roots.lhs,
                roots.rhs,
                &fastlec_input_config);
        if (!input_validation.ok) {
            throw AigError("Input validation failed: " + input_validation.message);
        }
        std::cerr << "[ec_expand] input validation passed: "
                  << input_validation.message << '\n';
        //开始扩展
        const ExpandedCircuit expanded = expand_circuit(source, roots.lhs, roots.rhs);
        write_aag(expanded.aag, output_aag);

        if (output_aag != options.output) {
            run_aigtoaig(options.aigtoaig, output_aag, options.output, options.verbose);
        }

        const ValidationResult output_validation =
            validate_output_single_eq_pair(options.output,
                                           expanded.left,
                                           expanded.right,
                                           &fastlec_output_config);
        if (!output_validation.ok) {
            std::cerr << "[ec_expand] output validation warning: "
                      << output_validation.message << '\n';
        } else {
            std::cerr << "[ec_expand] output validation passed: "
                      << output_validation.message << '\n';
        }

        if (options.verbose) {
            std::cerr << "[ec_expand] wrote " << options.output << '\n';
        }

        if (!options.keep_temp) {
            if (input_aag != options.input) {
                remove_if_exists(input_aag);
            }
            if (output_aag != options.output) {
                remove_if_exists(output_aag);
            }
            remove_if_exists(fastlec_input_log);
            if (output_validation.ok) {
                remove_if_exists(fastlec_output_log);
            }
        } else {
            std::cerr << "[ec_expand] kept temp input AAG: " << input_aag << '\n'
                      << "[ec_expand] kept temp output AAG: " << output_aag << '\n';
            if (!fastlec_input_log.empty()) {
                std::cerr << "[ec_expand] kept fastLEC input log: "
                          << fastlec_input_log << '\n';
            }
            std::cerr << "[ec_expand] kept fastLEC output log: "
                      << fastlec_output_log << '\n';
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ec_expand: " << ex.what() << '\n';
        return 1;
    }
}
