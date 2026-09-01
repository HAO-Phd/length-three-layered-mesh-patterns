# Length-three layered mesh patterns

This repository contains the manuscript and auxiliary reproducibility material
for **The complete equidistribution classification of length-three mesh
patterns on layered permutations.**

The paper introduces layer-constraint signatures for mesh patterns with
layered underlying permutation. At length three, the mathematical reduction is

```text
644 signatures -> 523 pointwise classes -> 286 reverse-complement orbits
               -> 200 classes after the complete unique-support reduction
               -> 164 coordinate classes
               -> 146 nonzero distribution classes + the zero class = 147 classes.
```

All identities valid for arbitrary host size are proved symbolically in the
manuscript. The short supplement records only the finite coordinate incidence
data, the eighteen-bridge inventory, and the pairwise separation statement.
The programs are optional exact cross-checks of the finite arithmetic and
bookkeeping; they are not a substitute for the mathematical proofs.

## Manuscript

- `main_revised.tex` and `main_revised.pdf`: main manuscript source and compiled PDF.
- `supplementary_material.tex` and `supplementary_material.pdf`: short supplementary proof details.
- `coordinate_connections_detailed.tex`: the 43 oriented identities typeset as Supplementary Table S2.

## Optional classification data and checks

- `equidistribution_classes_detailed.tex`: human-readable 146-class membership ledger; retained for lookup and auditing but not typeset in the short supplement.
- `pointwise_groups.json`: the 523 exact pointwise classes.
- `equidistribution_audit.json`: the 146 nonzero distribution classes.
- `coordinate_witnesses.json`: explicit witnesses for the coordinate edges.
- `lower_bound_separation.json`: exact pairwise lower-bound separators.
- `verify_classification.py`: primary exact audit and certificate generator.
- `independent_math_audit.py`: separately implemented author-supplied cross-check; it is not external independent validation.
- `REPRODUCIBILITY.md`: detailed certificate description and commands.

Run the committed-data checks with Python 3.10 or later:

```bash
python verify_classification.py
python independent_math_audit.py
python lookup_mesh_class.py --self-test
```

To regenerate the JSON data:

```bash
python verify_classification.py --write-certificates
```

## Pattern-to-class lookup

`lookup_mesh_class.py` converts a raw length-three mesh pattern into its active
signature and distribution class, lists the other signatures in that class,
and can export every corresponding raw mesh diagram. For example:

```bash
python lookup_mesh_class.py --perm 132 --shade "0,0;1,1;1,2;3,1"
python lookup_mesh_class.py --class-id 19
python lookup_mesh_class.py --self-test
```

See `LOOKUP_GUIDE.md` for input conventions, JSON output, and raw-member export.
