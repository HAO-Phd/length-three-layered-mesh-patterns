# Length-three layered mesh patterns

This repository contains the manuscript and reproducibility material for
**“The complete equidistribution classification of length-three mesh patterns
on layered permutations.”**

The paper introduces layer-constraint signatures for mesh patterns with
layered underlying permutation. At length three, the classification pipeline is

```text
644 signatures -> 523 pointwise classes -> 286 reverse-complement orbits
               -> 200 total-occurrence classes -> 164 support-refined classes
               -> 146 nonzero distribution classes + the zero class = 147 classes.
```

All verification scripts use exact integer or rational arithmetic; no random
sampling or floating-point comparison is used.

## Manuscript

- `main_revised.tex`: complete LaTeX source.
- `main_revised.pdf`: compiled manuscript.
- `equidistribution_classes_detailed.tex`: detailed class-membership table
  included by the main source.

## Verification and certificates

- `verify_classification.py`: primary exact audit and certificate generator.
- `independent_math_audit.py`: independent implementation of the central
  numerical checks.
- `pointwise_groups.json`: the 523 pointwise classes.
- `equidistribution_audit.json`: the 146 nonzero distribution classes and the
  local-orbit ordering used in Table 2.
- `coordinate_witnesses.json`: explicit coordinate witnesses for Table 2.
- `lower_bound_separation.json`: exact pairwise lower-bound separators.
- `REPRODUCIBILITY.md`: detailed certificate description and commands.

Run the committed-certificate audit with Python 3.10 or later:

```bash
python verify_classification.py
python independent_math_audit.py
```

To regenerate the JSON certificates:

```bash
python verify_classification.py --write-certificates
```

## Pattern-to-class lookup

`lookup_mesh_class.py` converts a raw length-three mesh pattern into its active
signature and Table 4 class, lists the other signatures in that class, and can
export every corresponding raw mesh diagram. For example:

```bash
python lookup_mesh_class.py --perm 132 --shade "0,0;1,1;1,2;3,1"
python lookup_mesh_class.py --class-id 19
python lookup_mesh_class.py --self-test
```

See `LOOKUP_GUIDE.md` for input conventions, JSON output, and raw-member export.

