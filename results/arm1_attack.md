# ARM 1 — A1 AND A3 ON THE `acc_ceiling` GRID

chr21 ONLY. K = 86 per-target leave-one-out. 44 targets x 7 labels x 10 draws per hypothesis per condition.
Weighting `uniform`; chr21 `beta_share` per target (~0.01502, NOT 1.0). Code revision `f76d3db`.
Methodology: `docs/plan.md` direction block, 2026-10-24; spec `docs/memory.md` 2026-10-27 (B-ATK-arm1-spec).

⚠⚠ POSING: FIXED MECHANISM — one mechanism (p's panel, chain, beta, eta), two input paths. That is what Theorem 1 bounds. The PROCEDURE posing (`Q_q` from q's own leave-one-out panel) is TWO mechanisms, is NOT bounded, and is not run.
⚠⚠ `acc_ceiling` IS A ROW, PER (target, label, condition) — the directory label is the D1_extremal ceiling ALONE. Reading a D2/D3 accuracy against the label would test it against a bound ~0.12 too loose.
⚠ The pooled accuracy is an AVERAGE over targets, so the bound it is tested against is the MEAN of the per-pair ceilings. The max is printed beside it; testing against the max would be a looser test that could hide an excess.
⚠ Pinsker is an INEQUALITY: every ceiling here is an OVER-estimate, so the true one is closer to 0.5 and an excess is conservative.
⚠ The 0.50 column is `eta` = 0, where `Q_p = Q_q = R_D` identically: the true accuracy is EXACTLY 0.5 for every attacker. It is a HARNESS CHECK.


## A1 — `lrt_pair`, the Neyman-Pearson OPTIMAL distinguisher

### D1_extremal

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| balanced accuracy* | 0.500000 | 0.534091 | 0.584091 | 0.613636 | 0.670455 | 0.689773 | 0.706818 |
| SE | 0.016855 | 0.016776 | 0.016615 | 0.016388 | 0.015842 | 0.015593 | 0.015339 |
| TPR (P_p[Lambda>0]) | 0.5000 | 0.5000 | 0.5818 | 0.6409 | 0.6795 | 0.6864 | 0.7205 |
| TNR (P_q[Lambda<0]) | 0.5000 | 0.5682 | 0.5864 | 0.5864 | 0.6614 | 0.6932 | 0.6932 |
| `acc_ceiling` THIS pair, mean — THE BOUND* | 0.500000 | 0.550000 | 0.600000 | 0.650000 | 0.700000 | 0.750000 | 0.800000 |
| `acc_ceiling` THIS pair, max (conservative) | 0.500000 | 0.550000 | 0.600000 | 0.650000 | 0.700000 | 0.750000 | 0.800000 |
| excess over the mean | +0.000000 | -0.015909 | -0.015909 | -0.036364 | -0.029545 | -0.060227 | -0.093182 |
| z | +0.00 | -0.95 | -0.96 | -2.22 | -1.86 | -3.86 | -6.07 |
| verdict | within | within | within | within | within | within | within |

- decoy at 0.80: `constructed`; `d_logz` median +261.243

### D2_nearest

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| balanced accuracy* | 0.500000 | 0.543182 | 0.528409 | 0.570455 | 0.584091 | 0.618182 | 0.636364 |
| SE | 0.016855 | 0.016788 | 0.016826 | 0.016687 | 0.016615 | 0.016377 | 0.016215 |
| TPR (P_p[Lambda>0]) | 0.5000 | 0.5318 | 0.5364 | 0.5682 | 0.5818 | 0.6227 | 0.6318 |
| TNR (P_q[Lambda<0]) | 0.5000 | 0.5545 | 0.5205 | 0.5727 | 0.5864 | 0.6136 | 0.6409 |
| `acc_ceiling` THIS pair, mean — THE BOUND* | 0.500000 | 0.531602 | 0.563228 | 0.594878 | 0.626552 | 0.658250 | 0.689971 |
| `acc_ceiling` THIS pair, max (conservative) | 0.500000 | 0.537388 | 0.574803 | 0.612247 | 0.649718 | 0.687217 | 0.724745 |
| excess over the mean | +0.000000 | +0.011580 | -0.034819 | -0.024424 | -0.042461 | -0.040068 | -0.053607 |
| z | +0.00 | +0.69 | -2.07 | -1.46 | -2.56 | -2.45 | -3.31 |
| verdict | within | within | within | within | within | within | within |

- decoy at 0.80: `HG00621`, concordance 0.7294 median; `d_logz` median -4.086

### D3_random

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| balanced accuracy* | 0.500000 | 0.515909 | 0.554545 | 0.562500 | 0.577273 | 0.622727 | 0.659091 |
| SE | 0.016855 | 0.016817 | 0.016754 | 0.016713 | 0.016613 | 0.016336 | 0.015978 |
| TPR (P_p[Lambda>0]) | 0.5000 | 0.4864 | 0.5568 | 0.5795 | 0.6114 | 0.6318 | 0.6545 |
| TNR (P_q[Lambda<0]) | 0.5000 | 0.5455 | 0.5523 | 0.5455 | 0.5432 | 0.6136 | 0.6636 |
| `acc_ceiling` THIS pair, mean — THE BOUND* | 0.500000 | 0.533224 | 0.566469 | 0.599734 | 0.633020 | 0.666327 | 0.699654 |
| `acc_ceiling` THIS pair, max (conservative) | 0.500000 | 0.537479 | 0.574982 | 0.612507 | 0.650056 | 0.687626 | 0.725219 |
| excess over the mean | +0.000000 | -0.017315 | -0.011923 | -0.037234 | -0.055747 | -0.043600 | -0.040563 |
| z | +0.00 | -1.03 | -0.71 | -2.23 | -3.36 | -2.67 | -2.54 |
| verdict | within | within | within | within | within | within | within |

- decoy at 0.80: `NA21309`, concordance 0.7082 median; `d_logz` median -2.811


## A3 — PanMixer's rarity-weighted gap score, as a 1-of-2 discriminator

### D1_extremal

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| balanced accuracy* | 0.500000 | 0.500000 | 0.500000 | 0.500000 | 0.500000 | 0.500000 | 0.500000 |
| SE | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| TPR (P_p[Lambda>0]) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| TNR (P_q[Lambda<0]) | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| `acc_ceiling` THIS pair, mean — THE BOUND* | 0.500000 | 0.550000 | 0.600000 | 0.650000 | 0.700000 | 0.750000 | 0.800000 |
| `acc_ceiling` THIS pair, max (conservative) | 0.500000 | 0.550000 | 0.600000 | 0.650000 | 0.700000 | 0.750000 | 0.800000 |
| excess over the mean | +0.000000 | -0.050000 | -0.100000 | -0.150000 | -0.200000 | -0.250000 | -0.300000 |
| z | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| verdict | DEGENERATE | DEGENERATE | DEGENERATE | DEGENERATE | DEGENERATE | DEGENERATE | DEGENERATE |

- decoy at 0.80: `constructed`; `d_logz` median +261.243

### D2_nearest

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| balanced accuracy* | 0.501136 | 0.521591 | 0.481818 | 0.517045 | 0.528409 | 0.531818 | 0.565909 |
| SE | 0.016474 | 0.016579 | 0.016328 | 0.016481 | 0.016321 | 0.016414 | 0.016380 |
| TPR (P_p[Lambda>0]) | 0.6068 | 0.6091 | 0.6045 | 0.6205 | 0.6500 | 0.6409 | 0.6636 |
| TNR (P_q[Lambda<0]) | 0.3955 | 0.4341 | 0.3591 | 0.4136 | 0.4068 | 0.4227 | 0.4682 |
| `acc_ceiling` THIS pair, mean — THE BOUND* | 0.500000 | 0.531602 | 0.563228 | 0.594878 | 0.626552 | 0.658250 | 0.689971 |
| `acc_ceiling` THIS pair, max (conservative) | 0.500000 | 0.537388 | 0.574803 | 0.612247 | 0.649718 | 0.687217 | 0.724745 |
| excess over the mean | +0.001136 | -0.010011 | -0.081410 | -0.077833 | -0.098143 | -0.126431 | -0.124062 |
| z | +0.07 | -0.60 | -4.99 | -4.72 | -6.01 | -7.70 | -7.57 |
| verdict | within | within | within | within | within | within | within |

- decoy at 0.80: `HG00621`, concordance 0.7294 median; `d_logz` median -4.086

### D3_random

| metric | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|---|
| balanced accuracy* | 0.495455 | 0.496591 | 0.509091 | 0.522727 | 0.538636 | 0.545455 | 0.565909 |
| SE | 0.016762 | 0.016732 | 0.016775 | 0.016745 | 0.016748 | 0.016740 | 0.016579 |
| TPR (P_p[Lambda>0]) | 0.4432 | 0.4364 | 0.4614 | 0.4705 | 0.4977 | 0.5091 | 0.5045 |
| TNR (P_q[Lambda<0]) | 0.5477 | 0.5568 | 0.5568 | 0.5750 | 0.5795 | 0.5818 | 0.6273 |
| `acc_ceiling` THIS pair, mean — THE BOUND* | 0.500000 | 0.533224 | 0.566469 | 0.599734 | 0.633020 | 0.666327 | 0.699654 |
| `acc_ceiling` THIS pair, max (conservative) | 0.500000 | 0.537479 | 0.574982 | 0.612507 | 0.650056 | 0.687626 | 0.725219 |
| excess over the mean | -0.004545 | -0.036633 | -0.057378 | -0.077007 | -0.094384 | -0.120872 | -0.133745 |
| z | -0.27 | -2.19 | -3.42 | -4.60 | -5.64 | -7.22 | -8.07 |
| verdict | within | within | within | within | within | within | within |

- decoy at 0.80: `NA21309`, concordance 0.7082 median; `d_logz` median -2.811


## SYNTHESIS

- **0 cell(s) exceed their own pair `acc_ceiling` at z > 1.96** and **7 are DEGENERATE** (min(TPR, TNR) < 0.05).
  - DEGENERATE A3 D1_extremal @ 0.50: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
  - DEGENERATE A3 D1_extremal @ 0.55: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
  - DEGENERATE A3 D1_extremal @ 0.60: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
  - DEGENERATE A3 D1_extremal @ 0.65: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
  - DEGENERATE A3 D1_extremal @ 0.70: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
  - DEGENERATE A3 D1_extremal @ 0.75: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
  - DEGENERATE A3 D1_extremal @ 0.80: TPR 1.0000 / TNR 0.0000 — a CONSTANT ANSWER, not a null. It must not be read as the mechanism holding.
- ⚠ An attacker that fails is not evidence the bound holds (standing constraint 4). A1 is the Neyman-Pearson optimum for its pair, so ITS null is the strongest available null — but A3's is not, and a DEGENERATE A3 row is no evidence at all.
- ⚠ D2/D3 decoys are cohort members IN the K = 86 panel while the target is HELD OUT, so the decoy is exactly reproducible by the mosaic and the target is not. That asymmetry is legitimate and inside the bounds, but it drives TPR and TNR apart; read those rows, not the balanced figure alone.
