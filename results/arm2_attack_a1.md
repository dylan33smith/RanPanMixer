methodology: docs/plan.md direction block, 2026-10-24 (git_rev 6979b5a4a1f096ff39a413271e4cd11c2e8bba16-dirty+844e5210ca16)

⚠⚠ ETA MODE: driving the tilt by eta [0.0, 75.7, 151.4, 227.1, 302.8, 378.5, 454.19] directly. The equivalent tau is 1.0 to float64 at every point, so p_succ_bound carries NO information here -- compare against the exact acc_ceiling in logs/B_eta_sweep/ instead.

10 targets x 25 draws x 7 tau = 1,750 trials x 4 release(s) = 7,000 releases
250 trials per (tau, condition) per ARM -> minimum detectable advantage delta = 0.0556 (BALANCED 0.879/sqrt(n))
30 workers (0.94 GB/job, 174 GB usable) — CAPPED by --workers 30 (planned 46: CORE-bound: 48 cores less 2 reserved)

  175/1750 trials
  350/1750 trials
  525/1750 trials
  700/1750 trials
  875/1750 trials
  1050/1750 trials
  1225/1750 trials
  1400/1750 trials
  1575/1750 trials
  1750/1750 trials
Accuracy of the optimal distinguisher. Metrics are rows; `tau` is the column.

| metric | 0.0 | 75.7 | 151.4 | 227.1 | 302.8 | 378.5 | 454.19 |
|---|---|---|---|---|---|---|---|
| D1_extremal balanced accuracy * | 0.5000 | 0.5140 | 0.5660 | 0.6280 | 0.6540 | 0.7260 | 0.7080 |
| D2_nearest balanced accuracy * | 0.5000 | 0.5120 | 0.5340 | 0.5520 | 0.5840 | 0.5900 | 0.6260 |
| D3_random balanced accuracy * | 0.5000 | 0.5580 | 0.5500 | 0.5380 | 0.5860 | 0.6080 | 0.6080 |
| D1_extremal one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5080 | 0.6080 | 0.5880 | 0.6480 | 0.7440 | 0.6880 |
| D2_nearest one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5080 | 0.5640 | 0.5080 | 0.5640 | 0.6200 | 0.5960 |
| D3_random one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5200 | 0.5400 | 0.5560 | 0.5560 | 0.5840 | 0.6160 |
| equivalent `tau` | 0 | 1 - 1e-65.45 | 1 - 1e-131.2 | 1 - 1e-197 | 1 - 1e-262.7 | 1 - 1e-328.5 | 1 - 1e-394.2 |
| `epsilon_tau` | 0 | 151.4 | 302.8 | 454.2 | 605.6 | 757 | 908.38 |
| `p_succ_bound` * | ~1.0 | ~1.0 | ~1.0 | ~1.0 | ~1.0 | ~1.0 | ~1.0 |
| bound comparison meaningful? | **no (vacuous)** | **no (vacuous)** | **no (vacuous)** | **no (vacuous)** | **no (vacuous)** | **no (vacuous)** | **no (vacuous)** |

⚠ Bound check is against the EXACT `acc_ceiling` supplied for [0.0, 75.7, 151.4, 227.1, 302.8, 378.5, 454.19], NOT `p_succ_bound` (which is 1.0 to float64 at every one of these eta and so cannot fail).

⚠ 1 cell(s) exceed `p_succ_bound` by POINT ESTIMATE but NOT significantly (one-sided 5%, Bonferroni over 21 cells, z_crit = 2.82). Reported, not claimed:
    D3_random at tau=75.7: 0.5580 vs bound 0.5500 (+0.0080, z = +0.36)

no condition exceeded the exact acc_ceiling at any eta

⚠ delta = 0.0556, so a cell resolves only where the exact acc_ceiling exceeds 0.5556. Below that the CEILING is the stronger statement and the attack adds nothing.
⚠ D2's decoy pool is the other prepared targets, not the 3,202-sample panel, so it is a LOWER BOUND on decoy difficulty and flatters us.
⚠ The three conditions share each release and must not be pooled.

JSON -> logs/B_acc_arm2/a1_acc.json
