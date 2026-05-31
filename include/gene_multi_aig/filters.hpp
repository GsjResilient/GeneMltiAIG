#pragma once

#include "gene_multi_aig/aag.hpp"

#include <functional>
#include <string>
#include <vector>

namespace gene_multi_aig {

enum class GateKind {
    Generic,
    And,
    Xor
};

struct GatePattern {
    GateKind kind = GateKind::Generic;
    Lit root = 0;
    std::vector<Lit> fanins;
};

struct FilterContext {
    const Aag& source;
    AigBuilder& builder;
    NetworkCopier& copier;
};

using FilterBuildFn = std::function<Lit(FilterContext&, const GatePattern&)>;

struct FilterDefinition {
    std::string name;
    GateKind kind = GateKind::Generic;
    FilterBuildFn build_left;
    FilterBuildFn build_right;
};

class FilterRegistry {
public:
    void add(FilterDefinition definition);
    const FilterDefinition& select(GateKind kind) const;

private:
    std::vector<FilterDefinition> definitions_;
};

GatePattern classify_gate_for_filter(const Aag& source, Lit root);
FilterRegistry make_default_filter_registry();

} // namespace gene_multi_aig
