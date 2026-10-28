# Does their shipped code bias our comparison? MEASURED.

chr21, our `acc_ceiling` 0.80 against their `capacity` 0.9, paired, n = 44, one build.
⚠ The SHIPPED row is the number to publish; the counterfactual says how much of it is the quirk.

## Why an identical rule is not a neutral rule

| cells changed vs the ORIGINAL row (median) | ours | theirs | ratio |
|---|---|---|---|
| allele swapped, both cells called | 111,615 | 112,922 | 0.99 |
| uncalled -> called | 20,985 | 1,708 | 12.29 |
| called -> uncalled | 8,244 | 1,240 | 6.65 |
| called/uncalled status flipped | 31,894 | 2,930 | 10.89 |

- **Genuine allele swaps are within 1%** of each other, but we flip the called/uncalled status of **10.9x** more cells, for **44/44** individuals — our release supplies calls the source assembly left uncalled. `af_loss` scores that as allele change.

## af_loss

| -1 policy | ours | theirs | we are lower | sign p | paired z |
|---|---|---|---|---|---|
| `shipped` | 0.003580 | 0.002841 | **0/44** | 1.14e-13 | +25.28 |
| `called` | 0.003654 | 0.003251 | **1/44** | 5.12e-12 | +14.98 |
| `genotype_only` | 0.002737 | 0.002763 | **31/44** | 0.00956 | -4.58 |

- Dropping -1 and renormalising narrows the gap from 0.000740 to 0.000403 (**46%** of it) without changing the direction.
- ⚠⚠ **BUT HOLDING THE MISSINGNESS SKELETON FIXED REVERSES IT.** Scoring genotype change alone, WE are lower in 31 of 44. So `af_loss` AS SHIPPED CANNOT SETTLE THIS COMPARISON: across three defensible -1 policies the verdict runs from 'they win on all 44 at p = 1e-13' to 'we win 31/44 at p = 0.01', and nothing about either mechanism changes between them. Report the shipped figure, and report that it is policy-dependent.

| release uncalled count | range | distinct values | correlation with the master row |
|---|---|---|---|
| ours | 16,656–16,666 | 5 | -0.3190 |
| theirs | 7,113–54,955 | 44 | +0.9991 |
| the master row itself | 7,561–55,191 | 44 | — |

- **That is the whole mechanism.** Ours is a fixed-support path, theirs is an edited copy of the target's own row. `af_loss` charges the side whose support diverges from the master, which is structurally us.

## ld_loss

| p_AB denominator | ours | theirs | theirs/ours | we are lower | sign p | paired z |
|---|---|---|---|---|---|---|
| `shipped` | 0.005405 | 0.040671 | **7.53x** | **40/44** | 1.71e-08 | -19.65 |
| `consistent` | 0.002146 | 0.002800 | **1.30x** | **44/44** | 1.14e-13 | -31.01 |

| share of the whole metric carried by its largest 6 pairs | ours | theirs |
|---|---|---|
| `shipped` | 0.603 | 0.940 |
| `consistent` | 0.101 | 0.155 |

- ⚠⚠ **OUR LD ADVANTAGE IS REAL BUT IT IS 1.30x, NOT 7.53x.** Under the shipped mixed denominator the largest **6 pairs of ~7,500** carry 94% of THEIR total, so that figure is set by a handful of pairs with missing cells and is not a mean over 7,462 pairs in any useful sense. Under a consistent denominator the advantage is smaller and MORE robust: **44/44** rather than 40/44.
