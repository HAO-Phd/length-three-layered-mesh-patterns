# Reproducibility certificate

All audits use exact integer or rational arithmetic. No floating-point
comparison, probabilistic hashing, or random sampling is used.

## Files

- `verify_classification.py` reconstructs the signature set, verifies the
  arbitrary-length signature enumeration and the all-layer pointwise
  partition, checks the reverse-complement and
  unique-support counts, regenerates the local-orbit map, checks the explicit
  coordinate witnesses, regenerates all pairwise lower-bound witnesses, and
  compares the signature formula with the original mesh definition through
  host size six.
- `pointwise_groups.json` records the 523 pointwise classes of the 644
  signatures.
- `equidistribution_audit.json` records the 146 nonzero distribution classes.
  In each class, `local_orbits` gives the explicit meaning of the local orbit
  labels used in Table 2 of the paper.
- `coordinate_witnesses.json` records, for every coordinate edge in Table 2,
  the source and target local orbits, explicit witness signatures, the map
  \(U_i\), and the direction of the identity.
- `lower_bound_separation.json` records, for each of the
  \(\binom{146}{2}=10{,}585\) pairs of nonzero classes, the smallest separating
  host size and both exact distribution polynomials.
- `lookup_mesh_class.py` and `LOOKUP_GUIDE.md` provide the pattern-to-class
  lookup described with the final membership table.

## Commands

Regenerate the JSON certificates:

```text
python verify_classification.py --write-certificates
```

Verify the committed certificates without rewriting them:

```text
python verify_classification.py
```

The pointwise audit compares the exact rational tensor series for every number
of host layers. The additional tensor check through four layers supplies
finite witnesses between distinct pointwise classes; it is not used to infer
all-layer equality.
