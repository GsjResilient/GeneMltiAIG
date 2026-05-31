#include "gene_multi_aig/filters.hpp"

#include <sstream>

namespace gene_multi_aig {
namespace {

Lit and3(AigBuilder& builder, Lit a, Lit b, Lit c) {
    return builder.add_and(builder.add_and(a, b), c);
}

Lit and4(AigBuilder& builder, Lit a, Lit b, Lit c, Lit d) {
    return builder.add_and(builder.add_and(a, b), builder.add_and(c, d));
}

Lit generic_left(FilterContext& ctx, const GatePattern& pattern) {
    const Lit x = ctx.copier.copy_lit(pattern.root, 1);
    const Lit y = ctx.copier.copy_lit(pattern.root, 2);
    const Lit xy = ctx.builder.add_and(x, y);
    return ctx.builder.add_and(x, invert_lit(xy));
}

Lit generic_right(FilterContext& ctx, const GatePattern& pattern) {
    const Lit x = ctx.copier.copy_lit(pattern.root, 1);
    const Lit y = ctx.copier.copy_lit(pattern.root, 2);
    const Lit x_or_y = ctx.builder.add_or(x, y);
    return ctx.builder.add_and(x_or_y, invert_lit(y));
}

Lit and_embedded_left(FilterContext& ctx, const GatePattern& pattern) {
    if (pattern.fanins.size() != 2U) {
        return generic_left(ctx, pattern);
    }

    const Lit x0 = ctx.copier.copy_lit(pattern.fanins[0], 1);
    const Lit x1 = ctx.copier.copy_lit(pattern.fanins[1], 1);
    const Lit y0 = ctx.copier.copy_lit(pattern.fanins[0], 2);
    const Lit y1 = ctx.copier.copy_lit(pattern.fanins[1], 2);

    const Lit term0 = and3(ctx.builder, x0, x1, invert_lit(y0));
    const Lit term1 = and3(ctx.builder, x1, x0, invert_lit(y1));
    return ctx.builder.add_or(term0, term1);
}

Lit and_embedded_right(FilterContext& ctx, const GatePattern& pattern) {
    if (pattern.fanins.size() != 2U) {
        return generic_right(ctx, pattern);
    }

    const Lit x0 = ctx.copier.copy_lit(pattern.fanins[0], 1);
    const Lit x1 = ctx.copier.copy_lit(pattern.fanins[1], 1);
    const Lit y0 = ctx.copier.copy_lit(pattern.fanins[0], 2);
    const Lit y1 = ctx.copier.copy_lit(pattern.fanins[1], 2);

    const Lit term0 = and3(ctx.builder, x1, invert_lit(y1), x0);
    const Lit term1 = and3(ctx.builder, invert_lit(y0), x0, x1);
    return ctx.builder.add_or(term0, term1);
}

Lit xor_embedded_left(FilterContext& ctx, const GatePattern& pattern) {
    if (pattern.fanins.size() != 2U) {
        return generic_left(ctx, pattern);
    }

    const Lit x0 = ctx.copier.copy_lit(pattern.fanins[0], 1);
    const Lit x1 = ctx.copier.copy_lit(pattern.fanins[1], 1);
    const Lit y0 = ctx.copier.copy_lit(pattern.fanins[0], 2);
    const Lit y1 = ctx.copier.copy_lit(pattern.fanins[1], 2);

    const Lit m0 = and4(ctx.builder, x0, invert_lit(x1), y0, y1);
    const Lit m1 = and4(ctx.builder, x0, invert_lit(x1), invert_lit(y0), invert_lit(y1));
    const Lit m2 = and4(ctx.builder, invert_lit(x0), x1, y0, y1);
    const Lit m3 = and4(ctx.builder, invert_lit(x0), x1, invert_lit(y0), invert_lit(y1));

    return ctx.builder.add_xor(ctx.builder.add_xor(m0, m1),
                               ctx.builder.add_xor(m2, m3));
}

Lit xor_embedded_right(FilterContext& ctx, const GatePattern& pattern) {
    if (pattern.fanins.size() != 2U) {
        return generic_right(ctx, pattern);
    }

    const Lit x0 = ctx.copier.copy_lit(pattern.fanins[0], 1);
    const Lit x1 = ctx.copier.copy_lit(pattern.fanins[1], 1);
    const Lit y0 = ctx.copier.copy_lit(pattern.fanins[0], 2);
    const Lit y1 = ctx.copier.copy_lit(pattern.fanins[1], 2);

    const Lit m0 = and4(ctx.builder, x1, invert_lit(x0), y1, y0);
    const Lit m1 = and4(ctx.builder, x1, invert_lit(x0), invert_lit(y1), invert_lit(y0));
    const Lit m2 = and4(ctx.builder, invert_lit(x1), x0, y1, y0);
    const Lit m3 = and4(ctx.builder, invert_lit(x1), x0, invert_lit(y1), invert_lit(y0));

    return ctx.builder.add_xor(ctx.builder.add_xor(m0, m2),
                               ctx.builder.add_xor(m1, m3));
}

std::string gate_kind_name(GateKind kind) {
    switch (kind) {
    case GateKind::Generic:
        return "generic";
    case GateKind::And:
        return "and";
    case GateKind::Xor:
        return "xor";
    }
    return "unknown";
}

} // namespace

void FilterRegistry::add(FilterDefinition definition) {
    definitions_.push_back(std::move(definition));
}

const FilterDefinition& FilterRegistry::select(GateKind kind) const {
    for (const FilterDefinition& definition : definitions_) {
        if (definition.kind == kind) {
            return definition;
        }
    }
    for (const FilterDefinition& definition : definitions_) {
        if (definition.kind == GateKind::Generic) {
            return definition;
        }
    }

    std::ostringstream oss;
    oss << "No filter registered for gate kind " << gate_kind_name(kind);
    throw AigError(oss.str());
}

GatePattern classify_gate_for_filter(const Aag& source, Lit root) {
    GatePattern pattern;
    pattern.root = root;

    if (const auto xor_match = match_xor(source, root)) {
        pattern.kind = GateKind::Xor;
        pattern.fanins = {xor_match->lhs, xor_match->rhs};
        return pattern;
    }

    if (const AndGate* gate = source.positive_and(root)) {
        pattern.kind = GateKind::And;
        pattern.fanins = {gate->rhs0, gate->rhs1};
        return pattern;
    }

    pattern.kind = GateKind::Generic;
    return pattern;
}

FilterRegistry make_default_filter_registry() {
    FilterRegistry registry;
    registry.add(FilterDefinition{
        "generic-no-xor",
        GateKind::Generic,
        generic_left,
        generic_right,
    });
    registry.add(FilterDefinition{
        "and-embedded-no-xor",
        GateKind::And,
        and_embedded_left,
        and_embedded_right,
    });
    registry.add(FilterDefinition{
        "xor-embedded",
        GateKind::Xor,
        xor_embedded_left,
        xor_embedded_right,
    });
    return registry;
}

} // namespace gene_multi_aig
