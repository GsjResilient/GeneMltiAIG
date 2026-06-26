#include "gene_multi_aig/aag.hpp"
#include "gene_multi_aig/validation.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using namespace gene_multi_aig;

namespace {

struct Options {
    fs::path input;
    fs::path output;
    fs::path tmp_dir;
    fs::path sub_out_dir;
    std::string aigtoaig = "aigtoaig";
    fs::path fastlec;
    unsigned fastlec_cores = 8;
    unsigned fastlec_timeout = 100;
    unsigned copies = 2;
    std::optional<Lit> lhs_lit;
    std::optional<Lit> rhs_lit;
    bool keep_temp = false;
    bool verbose = false;
    bool direct_copy_or = false;
};

struct ExpandedCircuit {
    Aag aag;
    Lit left = 0;
    Lit right = 0;
    std::vector<Lit> final_xor_internal_vars;
};

void print_usage(std::ostream& os) {
    os << "Usage: ec_expand --in input.aig --out expanded.aig [options]\n"
       << "\n"
       << "Options:\n"
       << "  --aigtoaig PATH       Path to external aigtoaig command\n"
       << "  --fastlec PATH        Path to external fastLEC command for required validation\n"
       << "  --fastlec-cores N     CPU threads for fastLEC validation [default: 8]\n"
       << "  --fastlec-timeout S   Timeout seconds for fastLEC validation [default: 100]\n"
       << "  --copies N            Number of independent copies [default: 2]\n"
       << "  --tmp-dir DIR         Directory for temporary AAG files\n"
       << "  --sub-out-dir DIR     Generate sub-AIGs by enumerating final XOR internals\n"
       << "  --keep-temp           Keep temporary AAG files and validation logs\n"
       << "  --direct-copy-or      Copy the input N times and OR the copied outputs\n"
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
        } else if (const auto value = option_value("--out")) {
            options.output = *value;
        } else if (const auto value = option_value("--aigtoaig")) {
            options.aigtoaig = *value;
        } else if (const auto value = option_value("--fastlec")) {
            options.fastlec = *value;
        } else if (const auto value = option_value("--fastlec-cores")) {
            options.fastlec_cores = parse_unsigned_arg(*value, "--fastlec-cores");
        } else if (const auto value = option_value("--fastlec-timeout")) {
            options.fastlec_timeout = parse_unsigned_arg(*value, "--fastlec-timeout");
        } else if (const auto value = option_value("--copies")) {
            options.copies = parse_unsigned_arg(*value, "--copies");
        } else if (const auto value = option_value("--tmp-dir")) {
            options.tmp_dir = *value;
        } else if (const auto value = option_value("--sub-out-dir")) {
            options.sub_out_dir = *value;
        } else if (arg == "--keep-temp") {
            options.keep_temp = true;
        } else if (arg == "--direct-copy-or") {
            options.direct_copy_or = true;
        } else if (const auto value = option_value("--lhs-lit")) {
            options.lhs_lit = parse_lit_arg(*value, "--lhs-lit");
        } else if (const auto value = option_value("--rhs-lit")) {
            options.rhs_lit = parse_lit_arg(*value, "--rhs-lit");
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
        throw AigError("Missing required --fastlec for validation/reporting");
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
    if (options.copies < 2U) {
        throw AigError("--copies must be at least 2");
    }
    if (!options.sub_out_dir.empty() && options.direct_copy_or) {
        throw AigError("--sub-out-dir is only supported in normal XOR expansion mode");
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

std::vector<std::vector<Lit>> create_input_copies(AigBuilder& builder,
                                                  std::size_t input_count,
                                                  unsigned copy_count) {
    std::vector<std::vector<Lit>> copies;
    copies.reserve(copy_count);

    for (unsigned copy = 0; copy < copy_count; ++copy) {
        std::vector<Lit> group;
        group.reserve(input_count);
        for (std::size_t i = 0; i < input_count; ++i) {
            group.push_back(builder.add_input());
        }
        copies.push_back(std::move(group));
    }

    return copies;
}

Lit build_balanced_or(AigBuilder& builder, std::vector<Lit> leaves) {
    if (leaves.empty()) {
        return 0;
    }

    while (leaves.size() > 1U) {
        std::vector<Lit> next;
        next.reserve((leaves.size() + 1U) / 2U);
        for (std::size_t i = 0; i < leaves.size(); i += 2U) {
            if (i + 1U < leaves.size()) {
                next.push_back(builder.add_or(leaves[i], leaves[i + 1U]));
            } else {
                next.push_back(leaves[i]);
            }
        }
        leaves = std::move(next);
    }

    return leaves.front();
}

Lit build_left_replacement(AigBuilder& builder, Lit a11, Lit a12, Lit b21, Lit b22) {
    const Lit l1 = builder.add_and(a12, invert_lit(b21));
    const Lit l2 = builder.add_and(a11, l1);
    const Lit l3 = builder.add_and(a12, invert_lit(b22));
    const Lit l4 = builder.add_and(a11, l3);
    return builder.add_or(l2, l4);
}

Lit build_right_replacement(AigBuilder& builder, Lit a21, Lit a22, Lit b11, Lit b12) {
    const Lit r1 = builder.add_and(a22, invert_lit(b11));
    const Lit r2 = builder.add_and(a21, r1);
    const Lit r3 = builder.add_and(a22, invert_lit(b12));
    const Lit r4 = builder.add_and(a21, r3);
    return builder.add_or(r2, r4);
}

std::pair<Lit, Lit> build_generic_replacements(NetworkCopier& copier,
                                                AigBuilder& builder,
                                                Lit lhs_root,
                                                Lit rhs_root) {
    const Lit lhs_group1 = copier.copy_lit(lhs_root, 1);
    const Lit rhs_group1 = copier.copy_lit(rhs_root, 1);
    const Lit lhs_group2 = copier.copy_lit(lhs_root, 2);
    const Lit rhs_group2 = copier.copy_lit(rhs_root, 2);
    const Lit left = builder.add_and(lhs_group1, invert_lit(rhs_group2));
    const Lit right = builder.add_and(rhs_group1, invert_lit(lhs_group2));
    return {left, right};
}

Lit append_and_gate(Aag& graph, Lit lhs, Lit rhs) {
    if (graph.max_var == UINT32_MAX) {
        throw AigError("AAG variable limit exceeded while adding sub-AIG constraints");
    }
    ++graph.max_var;
    const Lit out = graph.max_var * 2U;
    graph.ands.push_back(AndGate{out, lhs, rhs});
    return out;
}

Aag constrain_output_with_units(const Aag& base, const std::vector<Lit>& unit_literals) {
    if (base.outputs.size() != 1U) {
        throw AigError("Sub-AIG generation expects a single-output base AIG");
    }

    Aag constrained = base;
    Lit output = constrained.outputs.front();
    for (const Lit unit : unit_literals) {
        output = append_and_gate(constrained, output, unit);
    }
    constrained.outputs = {output};
    constrained.build_index();
    return constrained;
}

std::string assignment_bits(std::uint64_t mask, std::size_t width) {
    std::string bits(width, '0');
    for (std::size_t i = 0; i < width; ++i) {
        if ((mask & (1ULL << i)) != 0U) {
            bits[i] = '1';
        }
    }
    return bits;
}

fs::path sub_aig_path(const fs::path& output,
                      const fs::path& sub_out_dir,
                      const std::string& bits) {
    std::string stem = output.stem().string();
    if (stem.empty()) {
        stem = "expanded";
    }
    return sub_out_dir / (stem + "_sub_" + bits + ".aig");
}

void remove_if_exists(const fs::path& path);

void write_sub_aigs(const Aag& base,
                    const std::vector<Lit>& internal_vars,
                    const Options& options) {
    if (options.sub_out_dir.empty()) {
        return;
    }
    if (internal_vars.size() > 20U) {
        throw AigError("Refusing to generate more than 2^20 sub-AIG files");
    }

    fs::create_directories(options.sub_out_dir);

    const std::uint64_t total = 1ULL << internal_vars.size();
    for (std::uint64_t mask = 0; mask < total; ++mask) {
        std::vector<Lit> units;
        units.reserve(internal_vars.size());
        for (std::size_t i = 0; i < internal_vars.size(); ++i) {
            const bool value = (mask & (1ULL << i)) != 0U;
            units.push_back(value ? internal_vars[i] : invert_lit(internal_vars[i]));
        }

        const std::string bits = assignment_bits(mask, internal_vars.size());
        const fs::path sub_aig = sub_aig_path(options.output, options.sub_out_dir, bits);
        const fs::path sub_aag = make_temp_path(options.tmp_dir, "_sub.aag");
        write_aag(constrain_output_with_units(base, units), sub_aag);
        run_aigtoaig(options.aigtoaig, sub_aag, sub_aig, options.verbose);
        if (!options.keep_temp) {
            remove_if_exists(sub_aag);
        } else if (options.verbose) {
            std::cerr << "[ec_expand] kept temp sub AAG: " << sub_aag << '\n';
        }
    }

    std::cerr << "[ec_expand] wrote " << total
              << " sub-AIG files to " << options.sub_out_dir << '\n';
}

ExpandedCircuit expand_circuit(const Aag& source,
                               Lit lhs_root,
                               Lit rhs_root,
                               unsigned copy_count) {
    if (copy_count != 2U) {
        throw AigError("Normal replacement mode currently supports exactly 2 copies");
    }

    AigBuilder builder; //创建一个新电路构造器。后面所有新 input、新 gate 都会通过它生成，而不是直接修改原始 source
    std::vector<std::vector<Lit>> inputs =
        create_input_copies(builder, source.inputs.size(), copy_count);
    NetworkCopier copier(source, builder, std::move(inputs));

    Lit left = 0;
    Lit right = 0;
    const AndGate* lhs_gate = source.positive_and(lhs_root);
    const AndGate* rhs_gate = source.positive_and(rhs_root);
    if (lhs_gate != nullptr && rhs_gate != nullptr) {
        const Lit a11 = copier.copy_lit(lhs_gate->rhs0, 1);
        const Lit a12 = copier.copy_lit(lhs_gate->rhs1, 1);
        const Lit a21 = copier.copy_lit(rhs_gate->rhs0, 1);
        const Lit a22 = copier.copy_lit(rhs_gate->rhs1, 1);

        const Lit b11 = copier.copy_lit(lhs_gate->rhs0, 2);
        const Lit b12 = copier.copy_lit(lhs_gate->rhs1, 2);
        const Lit b21 = copier.copy_lit(rhs_gate->rhs0, 2);
        const Lit b22 = copier.copy_lit(rhs_gate->rhs1, 2);

        left = build_left_replacement(builder, a11, a12, b21, b22);
        right = build_right_replacement(builder, a21, a22, b11, b12);
    } else {
        const auto replacements =
            build_generic_replacements(copier, builder, lhs_root, rhs_root);
        left = replacements.first;
        right = replacements.second;
    }
    const XorBuildResult output = builder.add_xor_with_trace(left, right);

    return ExpandedCircuit{
        builder.to_aag({output.output}),
        left,
        right,
        output.internal_vars,
    };
}

ExpandedCircuit expand_direct_copy_or(const Aag& source, unsigned copy_count) {
    AigBuilder builder;
    std::vector<std::vector<Lit>> inputs =
        create_input_copies(builder, source.inputs.size(), copy_count);
    NetworkCopier copier(source, builder, std::move(inputs));

    std::vector<Lit> outputs;
    outputs.reserve(copy_count);
    for (unsigned copy = 1; copy <= copy_count; ++copy) {
        outputs.push_back(copier.copy_lit(source.outputs.front(), static_cast<int>(copy)));
    }

    const Lit output = build_balanced_or(builder, outputs);
    const Lit left = outputs.size() > 0U ? outputs.front() : 0;
    const Lit right = outputs.size() > 1U ? outputs[1] : left;
    return ExpandedCircuit{builder.to_aag({output}), left, right, {}};
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
        const fs::path fastlec_input_log = options.direct_copy_or
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

        ExpandedCircuit expanded;
        if (options.direct_copy_or) {
            if (options.verbose) {
                std::cerr << "[ec_expand] using direct-copy-or mode\n";
            }
            expanded = expand_direct_copy_or(source, options.copies);
        } else {
            if (options.copies != 2U) {
                throw AigError("Normal replacement mode currently supports exactly 2 copies");
            }

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
            expanded = expand_circuit(source, roots.lhs, roots.rhs, options.copies);
        }
        write_aag(expanded.aag, output_aag);

        if (output_aag != options.output) {
            run_aigtoaig(options.aigtoaig, output_aag, options.output, options.verbose);
        }

        bool keep_output_log = false;
        if (options.direct_copy_or) {
            const PotentialEquivClassCount output_count =
                count_output_potential_equiv_classes(options.output,
                                                     &fastlec_output_config);
            if (!output_count.ok) {
                keep_output_log = true;
                std::cerr << "[ec_expand] direct-copy-or fastLEC warning: "
                          << output_count.message << '\n';
            } else {
                std::cerr << "[ec_expand] direct-copy-or output potential-equivalence classes: "
                          << output_count.count << '\n';
            }
        } else {
            const PotentialEquivClassCount output_count =
                count_output_potential_equiv_classes(options.output,
                                                     &fastlec_output_config,
                                                     false);
            if (output_count.saw_bug || output_count.not_equivalent) {
                keep_output_log = true;
                throw AigError("Output validation failed: " + output_count.message);
            }
            if (!output_count.ok) {
                keep_output_log = true;
                std::cerr << "[ec_expand] output fastLEC warning: "
                          << output_count.message << '\n';
            } else {
                std::cerr << "[ec_expand] output potential-equivalence classes: "
                          << output_count.count << '\n';
            }
        }

        write_sub_aigs(expanded.aag, expanded.final_xor_internal_vars, options);

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
            if (!fastlec_input_log.empty()) {
                remove_if_exists(fastlec_input_log);
            }
            if (!keep_output_log) {
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
