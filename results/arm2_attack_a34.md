⚠⚠⚠ **[CORRECTION - 2026-10-28] THIS TABLE IS SUPERSEDED AND MUST NOT BE QUOTED.**
Its releases were drawn under an ARITHMETIC seed in `scripts/arm2_attack_a3a4.py`:
`SEED_BASE + i*10_000 + (ETA_SEED_OFFSET + ETA_GRID.index(eta))*1_000 + k`. Appending the
reported grid to `ETA_GRID` the same day pushed `index(eta)` to 10 and 11, where the two
terms ALIAS. **MEASURED on exactly this configuration: 225 colliding seeds, 450 of 1,750
cells = 25.71%, each `(target i, eta 378.5, draw k)` landing on the seed of
`(target i+1, eta 0.0, draw k)`.** The cells are therefore NOT independent and the
standard errors below are too narrow. See `docs/bugs.md`, occurrence SIX.
⚠⚠ **AND ITS COLUMNS ARE MISLABELLED INDEPENDENTLY OF THAT.** All ten targets ran at one
shared `eta_tau` list taken from HG00096's fit alone, while `c_quadratic` spreads 30.9274%
across them -- so nine of ten targets sat at a ceiling that is not the column's label.
✅ REPLACED BY a re-run at per-target `eta_tau` (`logs/B_acc_arm2/a34_ceil.md`), in which
the column IS the measured bound. This file is kept only until that lands.
⚠ The text below is preserved VERBATIM as the record of what was run.

methodology: docs/plan.md direction block, 2026-10-24 (git_rev 6979b5a4a1f096ff39a413271e4cd11c2e8bba16-dirty+ba903de77179)
A3 scores on 302,547 bi-allelic cohort sites of 332,521 chain positions.
computing the positive controls first -- a failed control voids the nulls...
  A3 unobfuscated link rate: 1.0000 (30/30)
  A4 planted-block detection: 1.0000 (30/30)

A3 power curve: 50 jobs x 4 untilted releases = 200 releases (30 workers)
  10/50 power
  20/50 power
  30/50 power
  40/50 power
  50/50 power
⚠⚠ ETA MODE: [0.0, 75.7, 151.4, 227.1, 302.8, 378.5, 454.19]. p_succ_bound is 1.0 to float64 at every one of these, so no bound comparison is printed; compare against the exact acc_ceiling in logs/B_eta_sweep/.


1,750 trials, each drawing 1 release under Q_p and 3 under Q_q (one per condition) = 7,000 releases
30 workers (0.60 GB/job, 173 GB usable) — CAPPED by --workers 30 (planned 46: CORE-bound: 48 cores less 2 reserved)

  87/1750 trials
  174/1750 trials
  261/1750 trials
  348/1750 trials
  435/1750 trials
  522/1750 trials
  609/1750 trials
  696/1750 trials
  783/1750 trials
  870/1750 trials
  957/1750 trials
  1044/1750 trials
  1131/1750 trials
  1218/1750 trials
  1305/1750 trials
  1392/1750 trials
  1479/1750 trials
  1566/1750 trials
  1653/1750 trials
  1740/1750 trials
methodology: docs/plan.md direction block, 2026-10-24 (git_rev 6979b5a4a1f096ff39a413271e4cd11c2e8bba16-dirty+ba903de77179)

ARM 2. 10 targets x 7 tau x 25 draws = 250 trials per (tau, condition) PER HYPOTHESIS -- matching A1's realised design of 250. chr21, scope chain, beta_share = 1.

### The pre-registered positive controls

| control | requirement | result | passed? |
|---|---|---|---|
| A3 (a) links the UNOBFUSCATED target | near-100% over 30 (target, decoy) pairs | **1.0000** (30/30) | **yes, but see below** |
| A4 finds a block copied VERBATIM from the target | the argmax block is the planted one, over 30 pairs | **1.0000** (30/30) | **yes** |

⚠⚠⚠ **A3's PRE-REGISTERED CONTROL CANNOT FAIL, SO ITS 1.0000 IS NOT A POWER STATEMENT.** Their score sums `-log P_HWE` only where genotypes AGREE, mismatches contribute zero, and the weights are non-negative, so `score(p, q)` is a SUBSET SUM of `score(p, p)` and `Lambda > 0` is algebraically forced. The control verifies the site axis, the genotype encoding and the weights -- the plumbing -- and nothing about power. A gate that cannot fail is decoration, so the power curve below replaces it as the actual power statement.
✅ **A4's control is NOT vacuous and it passed:** a 4,096-position block copied verbatim from the target was the argmax block in 30 of 30 pairs, against an untilted background. A verbatim leak of that size IS detectable by A4.

### A3's power against a KNOWN leak — a labelled DIAGNOSTIC, not an attack outcome

An untilted (`tau` = 0) release with a declared fraction `f` of its positions overwritten by the target VERBATIM. Balanced accuracy against `f`, so each column is A3's power against a leak of exactly that size. Subsets are prefixes of one random permutation, so the curve is monotone by construction.

| condition | f=0 | f=0.0001 | f=0.001 | f=0.01 | f=0.05 | f=0.1 | f=0.25 | f=0.5 | f=1 |
|---|---|---|---|---|---|---|---|---|---|
| D1_extremal | 0.5000 | 0.5000 | 0.5000 | 0.9900 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| D2_nearest | 0.5700 | 0.5700 | 0.5900 | 0.8000 | 0.9900 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| D3_random | 0.5200 | 0.5200 | 0.5300 | 0.8300 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

- **D1_extremal:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5000.

- **D2_nearest:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5700.

- **D3_random:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5200.

### A3 — rarity-weighted genotype linkage, re-posed as 1-of-2

Metrics are rows; `tau` is the column. Balanced accuracy, ties split.

| metric | 0.0 | 75.7 | 151.4 | 227.1 | 302.8 | 378.5 | 454.19 |
|---|---|---|---|---|---|---|---|
| D1_extremal balanced accuracy * | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| D2_nearest balanced accuracy * | 0.5000 | 0.5020 | 0.5300 | 0.5220 | 0.5480 | 0.5660 | 0.5380 |
| D3_random balanced accuracy * | 0.4940 | 0.4820 | 0.5280 | 0.5400 | 0.5360 | 0.5980 | 0.5620 |
| `p_succ_bound` * | 0.5000 | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) |

### A4 — segment-restricted likelihood ratio, maximised over blocks

⚠ EVERY declared block size is reported. The MAX over blocks is an optimistically biased estimate by construction -- maximising over `B` noisy accuracies -- so the Bonferroni-corrected comparison below is the only claimable one, and the raw max is shown beside it to make the size of that bias visible.

| metric | 0.0 | 75.7 | 151.4 | 227.1 | 302.8 | 378.5 | 454.19 |
|---|---|---|---|---|---|---|---|
| block 512: D1_extremal MAX over blocks | 0.5000 | 0.5050 | 0.5130 | 0.5080 | 0.5310 | 0.5110 | 0.5010 |
| block 512: D2_nearest MAX over blocks | 0.5000 | 0.5640 | 0.5620 | 0.5610 | 0.5780 | 0.5760 | 0.5720 |
| block 512: D3_random MAX over blocks | 0.5000 | 0.5570 | 0.5590 | 0.5540 | 0.5650 | 0.5990 | 0.5560 |
| block 4096: D1_extremal MAX over blocks | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| block 4096: D2_nearest MAX over blocks | 0.5000 | 0.5460 | 0.5420 | 0.5430 | 0.5480 | 0.5530 | 0.5590 |
| block 4096: D3_random MAX over blocks | 0.5000 | 0.5340 | 0.5420 | 0.5460 | 0.5540 | 0.5640 | 0.5630 |
| block 32768: D1_extremal MAX over blocks | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| block 32768: D2_nearest MAX over blocks | 0.5000 | 0.5420 | 0.5380 | 0.5660 | 0.5410 | 0.5520 | 0.5700 |
| block 32768: D3_random MAX over blocks | 0.5000 | 0.5230 | 0.5520 | 0.5460 | 0.5260 | 0.5460 | 0.5320 |
| `p_succ_bound` * | 0.5000 | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) |

**THE BOUND CHECK.** ⚠ SKIPPED where the point is an `eta_tau` (p_succ_bound is 1.0 to float64 there and cannot fail; compare against the exact `acc_ceiling` in logs/B_eta_sweep/). Otherwise one-sided on the excess over `p_succ_bound`, Bonferroni-corrected across the cells actually examined: A3 over 21 cells (z_crit = 2.82), A4 over 15,603 block-cells (z_crit = 4.51).


no condition exceeded `p_succ_bound` significantly, in either attack

✅ Both positive controls pass, so the nulls are lower bounds on distinguishability rather than broken instruments. ⚠ They remain LOWER bounds: the upper bound is `acc_ceiling`, which is proof-derived.
⚠ Minimum detectable advantage at 250 trials per arm per cell (25 draws x 10 targets) is 0.0556 UNCORRECTED, so tau = 0.05 and 0.10 do not resolve below the bound in either attack; those rows are reported, not claimed, and the PROOF carries them. A4's corrected threshold is 2.74x wider still.
⚠ D2's decoy pool is the other prepared targets, not the 3,202-sample panel, so it is a LOWER BOUND on decoy difficulty and flatters us.
⚠ The three conditions share each release and must not be pooled.
⚠ A3 scores on 302,547 bi-allelic cohort sites of 332,521 chain positions; the rest carry weight the attack cannot see, which flatters us.

JSON -> logs/B_acc_arm2/a34_acc.json
