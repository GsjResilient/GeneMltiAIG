#include "gene_multi_aig/sub_aig.hpp"

#include <fstream>
#include <unordered_set>
#include <utility>

namespace gene_multi_aig {
namespace {

std::uint32_t parse_var_token(const std::string& value,
                              const std::filesystem::path& path,
                              std::size_t line_number) {
    std::size_t consumed = 0;
    unsigned long parsed = 0;
    try {
        parsed = std::stoul(value, &consumed, 10);
    } catch (const std::exception&) {
        throw AigError("Invalid variable in " + path.string() + ":" +
                       std::to_string(line_number) + ": " + value);
    }
    if (consumed != value.size() || parsed > UINT32_MAX) {
        throw AigError("Invalid variable in " + path.string() + ":" +
                       std::to_string(line_number) + ": " + value);
    }
    return static_cast<std::uint32_t>(parsed);
}

std::string trim_ws(std::string value) {
    const char* whitespace = " \t\r\n";
    const std::size_t first = value.find_first_not_of(whitespace);
    if (first == std::string::npos) {
        return {};
    }
    const std::size_t last = value.find_last_not_of(whitespace);
    return value.substr(first, last - first + 1U);
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

} // namespace

std::vector<Lit> read_sub_variables(const Aag& base, const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw AigError("Failed to open sub-variable file: " + path.string());
    }

    std::vector<Lit> variables;
    std::unordered_set<std::uint32_t> seen;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }
        line = trim_ws(std::move(line));
        if (line.empty()) {
            continue;
        }

        const std::uint32_t var = parse_var_token(line, path, line_number);
        if (var == 0U) {
            throw AigError("Sub-variable file contains constant variable 0 at " +
                           path.string() + ":" + std::to_string(line_number));
        }
        if (var > UINT32_MAX / 2U) {
            throw AigError("Sub-variable is too large to encode as an AIGER literal at " +
                           path.string() + ":" + std::to_string(line_number));
        }
        if (!base.has_var(var)) {
            throw AigError("Sub-variable does not exist in base AIG at " +
                           path.string() + ":" + std::to_string(line_number) +
                           ": " + std::to_string(var));
        }
        if (!seen.insert(var).second) {
            throw AigError("Duplicate sub-variable in " + path.string() + ":" +
                           std::to_string(line_number) + ": " + std::to_string(var));
        }
        variables.push_back(var * 2U);
    }

    if (variables.empty()) {
        throw AigError("Sub-variable file did not contain any variables: " + path.string());
    }

    return variables;
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

std::filesystem::path sub_aig_path(const std::filesystem::path& input,
                                   const std::filesystem::path& out_dir,
                                   const std::string& bits) {
    std::string stem = input.stem().string();
    if (stem.empty()) {
        stem = "input";
    }
    return out_dir / (stem + "_sub_" + bits + ".aig");
}

} // namespace gene_multi_aig
