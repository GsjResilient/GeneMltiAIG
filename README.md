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
ec_expand --in input.aig --out expanded.aig --sub-out-dir /path/to/sub_aigs
ec_expand --in input.aig --out expanded.aig --lhs-lit 42 --rhs-lit 84
ec_expand --in input.aig --out expanded.aig --fastlec /path/to/fastLEC
```

`--lhs-lit` and `--rhs-lit` are a fallback when the tool cannot recognize the input PO as an XOR miter. Normally the input must have exactly one PO, and that PO must be a common AIG encoding of `lhs XOR rhs`.

Normal expansion preserves root literal polarity. If both miter roots are positive AND literals, `ec_expand` uses the embedded four-input replacement. If either root is inverted or otherwise not a positive AND literal, it uses a generic fallback that copies the literal with its polarity intact.

`--sub-out-dir` generates sub-AIG files after the expanded output is built. The tool records the internal variables created by the final `left XOR right`, enumerates every assignment to those variables, and writes one `.aig` per assignment. Since AIGER has no CNF clause section, each sub-AIG encodes the unit constraints by replacing the single output with `original_out AND assigned_lit_0 AND assigned_lit_1 ...`.

`--fastlec` is required. `ec_expand` runs `fastLEC` in `SAT_sweeping` mode before expansion, writes its output to a temporary log under `--tmp-dir`, and parses the log to require exactly one potential-equivalence class. If the input is not unique, the tool reports the fastLEC log path and stops. After constructing the output, it runs the same check on the generated AIGER and reports whether the output is unique; that post-check is informational and does not stop the run. `--keep-temp` keeps the logs for inspection.

## Batch Expansion

Use `scripts/batch_expand.sh` to expand a list of AIG files from one input directory:

```bash
scripts/batch_expand.sh \
  --input-dir ./cases \
  --list ./aig_list.txt \
  --out-dir ./expanded \
  --aigtoaig ../aiger/aigtoaig \
  --fastlec ../fastLEC/build/bin/fastLEC
```

The list file contains one AIG filename or relative path per line. Empty lines and lines starting with `#` are ignored. For `test_15_TOP32_72.aig`, the expanded output is written as `test_15_TOP32_72_expanded.aig`; relative subdirectories in the list are preserved under `--out-dir`.

Pass `--sub-out-dir ./sub_aigs` to also generate sub-AIG files for each expanded output.

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

Input validation calls `fastLEC` and stops the run unless the input has exactly one potential-equivalence class. Output validation also calls `fastLEC`, but only reports whether the generated AIGER has exactly one potential-equivalence class.

## Limitations

- Sequential AIGs are not supported. Inputs with latches are rejected.
- The parser supports standard ASCII AAG produced by `aigtoaig`.
- Uniqueness checks require a working `fastLEC` executable passed with `--fastlec`.
- If the input XOR miter uses an unusual AIG encoding, pass `--lhs-lit` and `--rhs-lit`.
