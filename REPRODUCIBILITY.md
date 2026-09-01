# Reproducibility notes

The mathematical proof is contained in the main manuscript and its short
supplement. All identities valid for arbitrary host size are proved
symbolically. The programs provide optional exact checks of finite coefficient
collection, class bookkeeping, and lower-bound separation. No floating-point
comparison, probabilistic hashing, random sampling, or extrapolation from
bounded host sizes is used.

## Files

- `supplementary_material.tex` contains the final quotient count, the
  eighteen-bridge inventory, all 43 oriented coordinate edges, and the finite
  pairwise-separation proposition.
- `coordinate_connections_detailed.tex` is included as Supplementary Table S2.
- `equidistribution_classes_detailed.tex` is the full human-readable
  146-class ledger. It is deliberately not typeset in the short supplement.
- `pointwise_groups.json` records the 523 pointwise classes of the 644 signatures.
- `equidistribution_audit.json` records the 146 nonzero distribution classes.
- `coordinate_witnesses.json` records the source and target orbit, witness signatures, and coordinate map for every edge in Supplementary Table S2.
- `lower_bound_separation.json` records, for each of the
  \(\binom{146}{2}=10{,}585\) pairs, the first separating host size and both exact distribution polynomials.
- `verify_classification.py` reconstructs the signature set, checks the exact quotient counts, audits the coordinate witnesses, regenerates the lower-bound records, and compares the signature formula with the original mesh definition on small hosts.
- `independent_math_audit.py` is a separately implemented author-supplied cross-check of the central numerical claims; it is not external independent validation.
- `lookup_mesh_class.py` and `LOOKUP_GUIDE.md` provide deterministic pattern-to-class lookup.

## Commands

Regenerate the JSON data:

```text
python verify_classification.py --write-certificates
```

Verify the committed data without rewriting it:

```text
python verify_classification.py
python independent_math_audit.py
python lookup_mesh_class.py --self-test
```

The pointwise audit compares exact rational tensor series for every number of
host layers. Its bounded tensor checks supply regression tests and finite
witnesses between distinct canonical forms; they are not used to infer
all-layer equality.

Likewise, the bounded coordinate-witness loop is a regression test. The
uniform identities for arbitrary numbers of layers are established by the
symbolic support-insertion argument in Appendix A of the manuscript.
