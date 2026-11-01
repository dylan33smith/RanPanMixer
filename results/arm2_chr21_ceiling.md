methodology: chr21-scope exact KL for all 4 decoy conditions, arm 2, weights 'uniform', beta_share 1.0, K = 88
⚠ `eta` READ from logs/B_acc_arm2/arm2_chr21_kl_verified.json (refined block), not re-solved.
⚠ `acc_ceiling` is Pinsker on the exact KL -- an UPPER BOUND on TV, not TV. This closes the COLUMN-LABEL envelope, not the TV ranking of the conditions.
  HG00096: D1_extremal=constructed, D2_nearest=HG00097, D3_random=NA19056, D4_farthest=NA19350
  HG00097: D1_extremal=constructed, D2_nearest=HG00096, D3_random=HG02389, D4_farthest=NA20294
  HG01323: D1_extremal=constructed, D2_nearest=HG04211, D3_random=NA20294, D4_farthest=NA20294
  HG01468: D1_extremal=constructed, D2_nearest=NA19056, D3_random=NA19350, D4_farthest=NA20294
  HG02389: D1_extremal=constructed, D2_nearest=NA19056, D3_random=HG00096, D4_farthest=NA19350
  HG03757: D1_extremal=constructed, D2_nearest=NA19056, D3_random=HG00097, D4_farthest=NA20294
  HG04211: D1_extremal=constructed, D2_nearest=HG01323, D3_random=NA19056, D4_farthest=NA20294
  NA19056: D1_extremal=constructed, D2_nearest=HG02389, D3_random=NA20294, D4_farthest=NA20294
  NA19350: D1_extremal=constructed, D2_nearest=NA19056, D3_random=HG01468, D4_farthest=NA20294
  NA20294: D1_extremal=constructed, D2_nearest=HG02389, D3_random=NA19350, D4_farthest=NA19350

10 targets x 6 labels x 4 conditions = 240 exact KL evaluations (60 jobs, 30 workers)
  10/60 jobs  (59s)
  20/60 jobs  (59s)
  30/60 jobs  (60s)
  40/60 jobs  (117s)
  50/60 jobs  (118s)
  60/60 jobs  (118s)

✅ instrument check: D1 reproduces the verified `kl_exact` on all 60 cells, worst |err| = 6.82e-13 (atol 1e-07)

Exact chr21 `kl_exact` by decoy condition. Metrics are ROWS, `acc_ceiling` labels are COLUMNS.

| metric | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|
| D1_extremal `kl_exact` median [n=10] | 0.020000 | 0.080000 | 0.180000 | 0.320000 | 0.500000 | 0.720000 |
| D2_nearest `kl_exact` median [n=10] | 0.006621 | 0.026504 | 0.059680 | 0.106178 | 0.166028 | 0.239260 |
| D3_random `kl_exact` median [n=10] | 0.007445 | 0.029795 | 0.067074 | 0.119305 | 0.186512 | 0.268718 |
| D4_farthest `kl_exact` median [n=10] | 0.007811 | 0.031257 | 0.070361 | 0.125146 | 0.195633 | 0.281843 |
| D1_extremal `acc_ceiling` median [n=10] | 0.5500 | 0.6000 | 0.6500 | 0.7000 | 0.7500 | 0.8000 |
| D2_nearest `acc_ceiling` median [n=10] | 0.5288 | 0.5576 | 0.5864 | 0.6152 | 0.6441 | 0.6729 |
| D3_random `acc_ceiling` median [n=10] | 0.5305 | 0.5610 | 0.5916 | 0.6221 | 0.6527 | 0.6833 |
| D4_farthest `acc_ceiling` median [n=10] | 0.5312 | 0.5625 | 0.5938 | 0.6251 | 0.6564 | 0.6877 |
| worst real/D1 `kl_exact` ratio [max over targets] * | 0.476656 | 0.477005 | 0.477354 | 0.477704 | 0.478053 | 0.478403 |

✅ D1_extremal carries the HIGHEST exact `kl_exact` in all 180 (target, label, real condition) cells on chr21, arm 2, so every chr21 column label bounds all three real decoys.

most-divergent REAL condition, by cell: D2_nearest 12/60, D3_random 18/60, D4_farthest 30/60
⚠ This ranks the three real decoys by exact KL, which is an upper-bound input and NOT the TV that the sampled attack ranks them on.

JSON -> logs/B_acc_arm2/arm2_chr21_ceiling.json
elapsed 118s
