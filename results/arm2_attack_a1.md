methodology: docs/plan.md direction block, 2026-10-24 (git_rev 6882ad55cac661e4c270070ac1ed725e0e99e3da-dirty+141712599732)

⚠⚠ CEILING MODE: each of the 10 targets runs at ITS OWN `eta_tau`, from logs/B_acc_arm2/arm2_chr21_kl_verified.json. The column is `acc_ceiling` and IS the bound -- measured at the solved eta and refined once, not fitted.
      0.50: eta 0.00 .. 0.00 (median 0.00)
      0.55: eta 75.35 .. 86.24 (median 76.41)
      0.60: eta 150.69 .. 172.45 (median 152.81)
      0.65: eta 226.03 .. 258.66 (median 229.21)
      0.70: eta 301.38 .. 344.84 (median 305.59)
      0.75: eta 376.72 .. 431.01 (median 381.98)
      0.80: eta 452.06 .. 517.17 (median 458.36)
⚠ `p_succ_bound` is 1.0 to float64 at every one of these eta, so it is NOT the bound here and is printed as vacuous.

10 targets x 100 draws x 7 tau = 7,000 trials x 5 release(s) = 35,000 releases
1000 trials per (tau, condition) per ARM -> minimum detectable advantage delta = 0.0278 (BALANCED 0.879/sqrt(n))
44 workers (0.94 GB/job, 180 GB usable) — CAPPED by --workers 44 (planned 46: CORE-bound: 48 cores less 2 reserved)

  700/7000 trials
  1400/7000 trials
  2100/7000 trials
  2800/7000 trials
  3500/7000 trials
  4200/7000 trials
  4900/7000 trials
  5600/7000 trials
  6300/7000 trials
  7000/7000 trials
Accuracy of the optimal distinguisher. Metrics are rows; `acc_ceiling` is the column.

| metric | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 |
|---|---|---|---|---|---|---|---|
| D1_extremal balanced accuracy * | 0.5000 [0.4781, 0.5219] | 0.5530 [0.5312, 0.5748] | 0.5895 [0.5679, 0.6111] | 0.5905 [0.5690, 0.6120] | 0.6555 [0.6347, 0.6763] | 0.6905 [0.6702, 0.7108] | 0.7190 [0.6993, 0.7387] |
| D2_nearest balanced accuracy * | 0.5000 [0.4781, 0.5219] | 0.5485 [0.5267, 0.5703] | 0.5360 [0.5141, 0.5579] | 0.5705 [0.5488, 0.5922] | 0.6200 [0.5987, 0.6413] | 0.6205 [0.5992, 0.6418] | 0.6290 [0.6078, 0.6502] |
| D3_random balanced accuracy * | 0.5000 [0.4781, 0.5219] | 0.5220 [0.5001, 0.5439] | 0.5440 [0.5222, 0.5658] | 0.5830 [0.5614, 0.6046] | 0.5935 [0.5720, 0.6150] | 0.6265 [0.6053, 0.6477] | 0.6425 [0.6215, 0.6635] |
| D4_farthest balanced accuracy * | 0.5000 [0.4781, 0.5219] | 0.5315 [0.5096, 0.5534] | 0.5300 [0.5081, 0.5519] | 0.5730 [0.5513, 0.5947] | 0.5735 [0.5518, 0.5952] | 0.6235 [0.6023, 0.6447] | 0.6505 [0.6296, 0.6714] |
| D1_extremal one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5450 | 0.5860 | 0.5960 | 0.6460 | 0.6800 | 0.7140 |
| D2_nearest one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5610 | 0.5320 | 0.5620 | 0.6160 | 0.6310 | 0.6430 |
| D3_random one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5310 | 0.5420 | 0.5880 | 0.5860 | 0.6230 | 0.6380 |
| D4_farthest one-sided P_p(L>0) [diagnostic, NOT bounded] | 0.5000 | 0.5410 | 0.5460 | 0.5770 | 0.5720 | 0.6450 | 0.6430 |
| `eta_tau` median [per target] | 0.00 | 76.41 | 152.81 | 229.21 | 305.59 | 381.98 | 458.36 |
| `eta_tau` min [per target] | 0.00 | 75.35 | 150.69 | 226.03 | 301.38 | 376.72 | 452.06 |
| `eta_tau` max [per target] | 0.00 | 86.24 | 172.45 | 258.66 | 344.84 | 431.01 | 517.17 |
| `epsilon_tau` median [per target] | 0.00 | 152.82 | 305.62 | 458.41 | 611.19 | 763.96 | 916.72 |
| equivalent `tau` | 0 | 1 - 1e-65 or smaller | 1 - 1e-65 or smaller | 1 - 1e-65 or smaller | 1 - 1e-65 or smaller | 1 - 1e-65 or smaller | 1 - 1e-65 or smaller |
| `p_succ_bound` * | ~1.0 | ~1.0 | ~1.0 | ~1.0 | ~1.0 | ~1.0 | ~1.0 |
| `acc_ceiling` [THE BOUND, measured per target] * | 0.5000 | 0.5500 | 0.6000 | 0.6500 | 0.7000 | 0.7500 | 0.8000 |

⚠ Bound check is against the COLUMN LABEL, which in this mode IS the exact `acc_ceiling`: each target ran at its own `eta_tau`, solved from its own `c` and verified by re-measurement. No --acc-ceiling is involved and none is accepted.

⚠ 1 cell(s) exceed their own measured `acc_ceiling` by POINT ESTIMATE but NOT significantly (one-sided 5%, Bonferroni over 28 cells, z_crit = 2.91), each with a 95% bar REACHING BELOW the line. Reported, not claimed:
    D1_extremal at acc_ceiling=0.55: 0.5530 [0.5312, 0.5748] vs bound 0.5500 (+0.0030, z = +0.27)

⚠ MULTIPLICITY. 28 cells were examined, so at a one-sided 5% level about 1.4 of them are expected above the line BY CHANCE even if the mechanism is exact. A bare count of 95% exceedances is not evidence; the Bonferroni-corrected z is.

no condition exceeded its own measured acc_ceiling at any ceiling

⚠ delta = 0.0278, so a cell resolves only where the exact acc_ceiling exceeds 0.5278. Below that the CEILING is the stronger statement and the attack adds nothing.
⚠ EVERY decoy pool is the other nine prepared targets, not the 3,202-sample panel, and it cuts OPPOSITE ways for the two real conditions: widening it would find a NEARER neighbour for D2, making that pair less distinguishable and the attacker weaker, and a FARTHER genome for D4. So D4 as run is a LOWER BOUND on the real worst case and is the limitation to report; the narrow D2 is conservative, not flattering.
⚠ D2_nearest is the LEAST distinguishable real pair and so the WEAKEST test of the bound here, not the hardest: the optimal distinguisher's advantage is TV(Q_p, Q_q), which grows with target-decoy discordance. 'Nearest is hardest' is a LINKAGE intuition and belongs to Arm 1.
⚠ The four conditions share each release and must not be pooled.

JSON -> logs/B_acc_arm2/a1_ceil_n1000_d4.json
