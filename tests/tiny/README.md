# tiny — the only fixture that can be committed

24 records, 4 samples (8 haplotypes), a 6-marker map. Hand-written so the expected
arrays can be written out by hand in the tests rather than golden-filed, which
would just re-record whatever the code did on the day.

Every record is here to exercise something that has no other coverage:

| what | where |
|---|---|
| plain called sites | 100, 200, 400 |
| an explicit `.\|.` no-call with no structure → **code 4** | 300 |
| a 4-allele record → `allele_lengths` ragged, `n_alt` > 1 | 520 |
| `LV=1` child of 520; S1/S4 carry the deletion → **code 1** | 522 |
| `LV=2` grandchild; parent's own call missing → **code 2** | 523 |
| `CONFLICT` naming S3,S4 → **code 3** | 700 |
| **the precedence collision**: S3 is both CONFLICT-named *and* nested-not-applicable. Conflict must win | 720 |
| a 6-record missing run for S4 → `missing_runs.tsv`, and **no per-cell code** | 900–1040 |
| **deliberately wrong `INFO/AN`** (says 2, truth is 8) → the row-alignment canary must fire | 1300 |
| below the map span (`left=nan`) | 100–400 |
| above the map span (`right=nan`) | 1700, 1900 |
| a zero-dcM interior interval (500→700) → a real plateau, correct to leave | 500–700 |

The wrong-`AN` record at 1300 is the important one. `support(v)` is read from
`INFO/AN` because it is bit-identical to the non-missing count on real data — which
makes the comparison a free row-alignment check. A check nobody has watched fail is
not a check, so one record is built to make it fail.
