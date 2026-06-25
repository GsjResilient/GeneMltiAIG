#include "gene_multi_aig/aag.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <utility>

namespace gene_multi_aig {
namespace {

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

Lit parse_lit_token(const std::string& token, const std::string& context) {
    std::size_t consumed = 0;
    unsigned long value = 0;
    try {
        value = std::stoul(token, &consumed, 10);
    } catch (const std::exception&) {
        throw AigError("Invalid literal in " + context + ": " + token);
    }
    if (consumed != token.size()) {
        throw AigError("Invalid literal in " + context + ": " + token);
    }
    if (value > UINT32_MAX) {
        throw AigError("Literal out of range in " + context + ": " + token);
    }
    return static_cast<Lit>(value);
}

bool same_lit_set(Lit a0, Lit a1, Lit b0, Lit b1) {
    return (a0 == b0 && a1 == b1) || (a0 == b1 && a1 == b0);
}

std::optional<std::pair<Lit, Lit>> match_and2(const Aag& graph, Lit lit) {
    const AndGate* gate = graph.positive_and(lit);
    if (gate == nullptr) {
        return std::nullopt;
    }
    return std::make_pair(gate->rhs0, gate->rhs1);
}

std::optional<std::pair<Lit, Lit>> match_or2(const Aag& graph, Lit lit) {
    if (!lit_is_inverted(lit)) {
        return std::nullopt;
    }
    const AndGate* gate = graph.and_by_var(lit_var(lit));
    if (gate == nullptr) {
        return std::nullopt;
    }
    return std::make_pair(invert_lit(gate->rhs0), invert_lit(gate->rhs1));
}

std::optional<XorMatch> match_xor_from_sop(const Aag& graph, Lit lit) {
    const auto outer_or = match_or2(graph, lit);
    if (!outer_or) {
        return std::nullopt;
    }

    const auto t0 = match_and2(graph, outer_or->first);
    const auto t1 = match_and2(graph, outer_or->second);
    if (!t0 || !t1) {
        return std::nullopt;
    }

    const Lit a[2] = {t0->first, t0->second};
    const Lit b[2] = {t1->first, t1->second};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            const Lit u = a[i];
            const Lit v = a[1 - i];
            if (invert_lit(u) == b[j] && invert_lit(v) == b[1 - j]) {
                return XorMatch{u, invert_lit(v)};
            }
        }
    }

    return std::nullopt;
}

std::optional<XorMatch> match_xor_from_and_form(const Aag& graph, Lit lit) {
    const auto outer_and = match_and2(graph, lit);
    if (!outer_and) {
        return std::nullopt;
    }

    const Lit candidates[2] = {outer_and->first, outer_and->second};
    for (int or_index = 0; or_index < 2; ++or_index) {
        const Lit or_lit = candidates[or_index];
        const Lit not_and_lit = candidates[1 - or_index];
        const auto or_terms = match_or2(graph, or_lit);
        if (!or_terms || !lit_is_inverted(not_and_lit)) {
            continue;
        }

        const auto and_terms = match_and2(graph, invert_lit(not_and_lit));
        if (!and_terms) {
            continue;
        }

        if (same_lit_set(or_terms->first, or_terms->second,
                         and_terms->first, and_terms->second)) {
            return XorMatch{or_terms->first, or_terms->second};
        }
    }

    return std::nullopt;
}

} // namespace

void Aag::build_index() {
    and_index_.clear();
    for (std::size_t i = 0; i < ands.size(); ++i) {
        and_index_[lit_var(ands[i].lhs)] = i;
    }
}

const AndGate* Aag::and_by_var(std::uint32_t var) const {
    const auto it = and_index_.find(var);
    if (it == and_index_.end()) {
        return nullptr;
    }
    return &ands[it->second];
}

const AndGate* Aag::positive_and(Lit lit) const {
    if (lit_is_inverted(lit) || lit == 0U || lit == 1U) {
        return nullptr;
    }
    return and_by_var(lit_var(lit));
}

bool Aag::has_var(std::uint32_t var) const {
    if (var == 0U || var > max_var) {
        return false;
    }
    for (const Lit input : inputs) {
        if (lit_var(input) == var) {
            return true;
        }
    }
    return and_by_var(var) != nullptr;
}

Aag read_aag(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw AigError("Failed to open AAG file: " + path.string());
    }

    std::string line;
    if (!std::getline(input, line)) {
        throw AigError("AAG file is empty: " + path.string());
    }

    const auto header = split_ws(line);
    if (header.size() < 6 || header[0] != "aag") {
        throw AigError("Unsupported or invalid AAG header in: " + path.string());
    }

    Aag result;
    result.max_var = parse_lit_token(header[1], "AAG header");
    const auto input_count = parse_lit_token(header[2], "AAG header");
    result.latch_count = parse_lit_token(header[3], "AAG header");
    const auto output_count = parse_lit_token(header[4], "AAG header");
    const auto and_count = parse_lit_token(header[5], "AAG header");

    result.inputs.reserve(input_count);
    result.outputs.reserve(output_count);
    result.ands.reserve(and_count);

    for (std::uint32_t i = 0; i < input_count; ++i) {
        if (!std::getline(input, line)) {
            throw AigError("Unexpected EOF while reading inputs");
        }
        const auto tokens = split_ws(line);
        if (tokens.size() != 1U) {
            throw AigError("Invalid input line: " + line);
        }
        result.inputs.push_back(parse_lit_token(tokens[0], "input line"));
    }

    for (std::uint32_t i = 0; i < result.latch_count; ++i) {
        if (!std::getline(input, line)) {
            throw AigError("Unexpected EOF while reading latches");
        }
    }

    for (std::uint32_t i = 0; i < output_count; ++i) {
        if (!std::getline(input, line)) {
            throw AigError("Unexpected EOF while reading outputs");
        }
        const auto tokens = split_ws(line);
        if (tokens.size() != 1U) {
            throw AigError("Invalid output line: " + line);
        }
        result.outputs.push_back(parse_lit_token(tokens[0], "output line"));
    }

    for (std::uint32_t i = 0; i < and_count; ++i) {
        if (!std::getline(input, line)) {
            throw AigError("Unexpected EOF while reading AND gates");
        }
        const auto tokens = split_ws(line);
        if (tokens.size() != 3U) {
            throw AigError("Invalid AND line: " + line);
        }
        AndGate gate;
        gate.lhs = parse_lit_token(tokens[0], "AND line");
        gate.rhs0 = parse_lit_token(tokens[1], "AND line");
        gate.rhs1 = parse_lit_token(tokens[2], "AND line");
        result.ands.push_back(gate);
    }

    result.build_index();
    return result;
}

void write_aag(const Aag& aag, const std::filesystem::path& path) {
    std::ofstream output(path);
    if (!output) {
        throw AigError("Failed to create AAG file: " + path.string());
    }

    output << "aag " << aag.max_var << ' ' << aag.inputs.size() << ' '
           << aag.latch_count << ' ' << aag.outputs.size() << ' '
           << aag.ands.size() << '\n';

    for (const Lit input : aag.inputs) {
        output << input << '\n';
    }
    for (const Lit out : aag.outputs) {
        output << out << '\n';
    }
    for (const AndGate& gate : aag.ands) {
        output << gate.lhs << ' ' << gate.rhs0 << ' ' << gate.rhs1 << '\n';
    }
    output << "c\n";
    output << "generated by GeneMltiAIG ec_expand\n";
}

std::optional<XorMatch> match_xor(const Aag& graph, Lit lit) {
    if (auto match = match_xor_from_and_form(graph, lit)) {
        return match;
    }
    if (auto match = match_xor_from_sop(graph, lit)) {
        return match;
    }
    return std::nullopt;
}

Lit AigBuilder::add_input() {
    const Lit lit = next_var_ * 2U;
    ++next_var_;
    inputs_.push_back(lit);
    return lit;
}

Lit AigBuilder::add_and(Lit lhs, Lit rhs) {
    const Lit out = next_var_ * 2U;
    ++next_var_;
    ands_.push_back(AndGate{out, lhs, rhs});
    return out;
}

Lit AigBuilder::add_or(Lit lhs, Lit rhs) {
    return invert_lit(add_and(invert_lit(lhs), invert_lit(rhs)));
}

Lit AigBuilder::add_xor(Lit lhs, Lit rhs) {
    return add_xor_with_trace(lhs, rhs).output;
}

XorBuildResult AigBuilder::add_xor_with_trace(Lit lhs, Lit rhs) {
    const Lit lhs_and_not_rhs = add_and(lhs, invert_lit(rhs));
    const Lit not_lhs_and_rhs = add_and(invert_lit(lhs), rhs);
    const Lit or_and = add_and(invert_lit(lhs_and_not_rhs), invert_lit(not_lhs_and_rhs));
    return XorBuildResult{
        invert_lit(or_and),
        {lhs_and_not_rhs, not_lhs_and_rhs, or_and},
    };
}

Aag AigBuilder::to_aag(const std::vector<Lit>& outputs) const {
    Aag aag;
    aag.max_var = max_var();
    aag.latch_count = 0;
    aag.inputs = inputs_;
    aag.outputs = outputs;
    aag.ands = ands_;
    aag.build_index();
    return aag;
}

NetworkCopier::NetworkCopier(const Aag& source,
                             AigBuilder& builder,
                             std::vector<std::vector<Lit>> group_inputs)
    : source_(source), builder_(builder), group_maps_(group_inputs.size()) {
    if (group_inputs.empty()) {
        throw AigError("At least one copy group is required");
    }

    for (std::size_t group = 0; group < group_inputs.size(); ++group) {
        if (group_inputs[group].size() != source.inputs.size()) {
            throw AigError("Input copy maps do not match source input count");
        }

        for (std::size_t i = 0; i < source.inputs.size(); ++i) {
            group_maps_[group][lit_var(source.inputs[i])] = group_inputs[group][i];
        }
    }
}

NetworkCopier::NetworkCopier(const Aag& source,
                             AigBuilder& builder,
                             std::vector<Lit> group1_inputs,
                             std::vector<Lit> group2_inputs)
    : NetworkCopier(source,
                    builder,
                    std::vector<std::vector<Lit>>{std::move(group1_inputs),
                                                  std::move(group2_inputs)}) {}

Lit NetworkCopier::copy_lit(Lit source_lit, int group) {
    if (source_lit == 0U || source_lit == 1U) {
        return source_lit;
    }
    if (group < 1 || static_cast<std::size_t>(group) > group_maps_.size()) {
        throw AigError("Copy group is out of range: " + std::to_string(group));
    }

    const Lit positive = copy_positive_var(lit_var(source_lit), group);
    return lit_is_inverted(source_lit) ? invert_lit(positive) : positive;
}

Lit NetworkCopier::copy_positive_var(std::uint32_t source_var, int group) {
    auto& map = group_maps_[static_cast<std::size_t>(group - 1)];
    const auto existing = map.find(source_var);
    if (existing != map.end()) {
        return existing->second;
    }

    const AndGate* source_gate = source_.and_by_var(source_var);
    if (source_gate == nullptr) {
        throw AigError("Literal references an unknown source variable: " +
                       std::to_string(source_var));
    }

    const Lit rhs0 = copy_lit(source_gate->rhs0, group);
    const Lit rhs1 = copy_lit(source_gate->rhs1, group);
    const Lit copied = builder_.add_and(rhs0, rhs1);
    map[source_var] = copied;
    return copied;
}

} // namespace gene_multi_aig
