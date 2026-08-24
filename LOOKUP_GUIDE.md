# Table 4 lookup guide

`lookup_mesh_class.py` is a deterministic reader for Table 4 of the paper. It
uses `equidistribution_audit.json`, the same 146-row data set from which the
table was produced. The script requires only Python 3.10 or later and the
standard library.

## Query a raw mesh pattern

Give the underlying permutation and the shaded cells. Separate cells with
semicolons:

```powershell
python lookup_mesh_class.py --perm 132 --shade "0,0;1,1;1,2;3,1"
```

The result reports the active signature, the Table 4 class, shaded inactive
cells that do not affect the statistic, every other signature in the same
class, and the number of raw mesh diagrams represented by the signature and
by the full class.

The worked example above has signature `(12;01;{0,1})` and belongs to Class
19. The cell `(1,2)` is inactive and is therefore ignored by the signature.

An empty shading may be queried by omitting `--shade`:

```powershell
python lookup_mesh_class.py --perm 321
```

The non-layered underlying patterns `231` and `312` are reported as members
of the zero class, which is not numbered in Table 4.

## Query a class directly

```powershell
python lookup_mesh_class.py --class-id 19
python lookup_mesh_class.py --class-id zero
```

This lists all signature members of the requested class and its total number
of raw mesh diagrams.

## Export all raw members

Use `--export-raw` with a `.csv`, `.jsonl`, or `.ndjson` destination:

```powershell
python lookup_mesh_class.py --class-id 19 --export-raw class_19.csv
```

Each exported row records the Table 4 class, underlying permutation, a
canonical 16-bit shading mask, the shaded-cell list, and the active signature.
The mask convention is that bit `4*x+y` represents cell `(x,y)`. A mask can
also be used as input:

```powershell
python lookup_mesh_class.py --perm 132 --mask 0x2061
```

Use the script's reported mask for a pattern rather than constructing it by
hand when possible.

## Machine-readable and AI-assisted use

Add `--json` to obtain stable machine-readable output:

```powershell
python lookup_mesh_class.py --perm 132 --shade "0,0;1,1;1,2;3,1" --json
```

An AI assistant may parse a natural-language description and call this
command, but the class assignment should be taken from the deterministic
script output rather than inferred by the language model itself.

## Consistency test

```powershell
python lookup_mesh_class.py --self-test
```

The self-test checks all `4*2^16=262144` raw diagrams with layered underlying
permutation, verifies that they realize exactly 644 signatures and all 146
Table 4 classes, adds the `2*2^16=131072` zero-class diagrams, and checks the
worked Class 19 example.


