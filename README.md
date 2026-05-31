# GeneMltiAIG

`GeneMltiAIG` builds an expanded EC miter from an input AIGER circuit that has exactly one equivalence pair: the two roots compared by a single XOR miter output.

The current tool is `ec_expand`. It uses the external AIGER command `aigtoaig` to convert binary `.aig` files to ASCII `.aag`, performs the construction on `.aag`, then converts the result back to `.aig`. The converter is invoked with explicit `-a` and `-b` mode flags.

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Usage

```bash
ec_expand --in input.aig --out expanded.aig --aigtoaig aigtoaig
```

Useful options:

```bash
ec_expand --in input.aig --out expanded.aig --keep-temp
ec_expand --in input.aig --out expanded.aig --tmp-dir /path/to/tmp
ec_expand --in input.aig --out expanded.aig --lhs-lit 42 --rhs-lit 84
```

`--lhs-lit` and `--rhs-lit` are a fallback when the tool cannot recognize the input PO as an XOR miter. Normally the input must have exactly one PO, and that PO must be a common AIG encoding of `lhs XOR rhs`.

## Construction

For an input miter root pair `(a, b)`, the tool creates two independent copies of the original PI set:

```text
x^1, x^2
```

It then builds:

```text
left  = P(a(x^1), a(x^2))
right = P2(b(x^1), b(x^2))
out   = left XOR right
```

All default filters implement the same Boolean function:

```text
x & ~y
```

The left and right filter implementations use different structures. Filters are registered by gate kind, so new filters can be added in `src/filters.cpp` without changing the main flow.

Default filters:

- `generic`: no-XOR fallback for literals that are not recognized as a handled gate kind.
- `and-embedded`: no-XOR embedded filter for positive AND roots. It avoids materializing the removed AND root as a standalone copied node.
- `xor-embedded`: XOR-containing embedded filter for recognized XOR roots.

## Validation Hooks

The tool calls two validation functions:

- `validateInputSingleEqPair(...)`
- `validateOutputSingleEqPair(...)`

They are currently placeholders in `src/validation.cpp` and return success. They are intentionally isolated so fastLEC can be added later without changing the construction pipeline.

## Limitations

- Sequential AIGs are not supported. Inputs with latches are rejected.
- The parser supports standard ASCII AAG produced by `aigtoaig`.
- The built-in validation hooks do not prove uniqueness of the equivalence pair yet.
- If the input XOR miter uses an unusual AIG encoding, pass `--lhs-lit` and `--rhs-lit`.
