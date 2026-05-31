#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace gene_multi_aig {

using Lit = std::uint32_t;

inline Lit invert_lit(Lit lit) {
    return lit ^ 1U;
}

inline std::uint32_t lit_var(Lit lit) {
    return lit >> 1U;
}

inline bool lit_is_inverted(Lit lit) {
    return (lit & 1U) != 0U;
}

struct AndGate {
    Lit lhs = 0;
    Lit rhs0 = 0;
    Lit rhs1 = 0;
};

struct Aag {
    std::uint32_t max_var = 0;
    std::uint32_t latch_count = 0;
    std::vector<Lit> inputs;
    std::vector<Lit> outputs;
    std::vector<AndGate> ands;

    void build_index();
    const AndGate* and_by_var(std::uint32_t var) const;
    const AndGate* positive_and(Lit lit) const;
    bool has_var(std::uint32_t var) const;

private:
    std::unordered_map<std::uint32_t, std::size_t> and_index_;
};

class AigError : public std::runtime_error {
public:
    explicit AigError(const std::string& message) : std::runtime_error(message) {}
};

Aag read_aag(const std::filesystem::path& path);
void write_aag(const Aag& aag, const std::filesystem::path& path);

struct XorMatch {
    Lit lhs = 0;
    Lit rhs = 0;
};

std::optional<XorMatch> match_xor(const Aag& graph, Lit lit);

class AigBuilder {
public:
    Lit add_input();
    Lit add_and(Lit lhs, Lit rhs);
    Lit add_or(Lit lhs, Lit rhs);
    Lit add_xor(Lit lhs, Lit rhs);

    Aag to_aag(const std::vector<Lit>& outputs) const;
    std::uint32_t max_var() const { return next_var_ - 1U; }

private:
    std::uint32_t next_var_ = 1;
    std::vector<Lit> inputs_;
    std::vector<AndGate> ands_;
};

class NetworkCopier {
public:
    NetworkCopier(const Aag& source,
                  AigBuilder& builder,
                  std::vector<Lit> group1_inputs,
                  std::vector<Lit> group2_inputs);

    Lit copy_lit(Lit source_lit, int group);

private:
    Lit copy_positive_var(std::uint32_t source_var, int group);

    const Aag& source_;
    AigBuilder& builder_;
    std::unordered_map<std::uint32_t, Lit> group_maps_[2];
};

} // namespace gene_multi_aig
