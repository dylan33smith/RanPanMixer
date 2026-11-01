methodology: docs/plan.md direction block, 2026-10-24 (git_rev 6882ad55cac661e4c270070ac1ed725e0e99e3da-dirty+7e60fa4b3995)
A3 scores on 302,547 bi-allelic cohort sites of 332,521 chain positions.
computing the positive controls first -- a failed control voids the nulls...
  A3 unobfuscated link rate: 1.0000 (40/40)
  A4 planted-block detection: 1.0000 (40/40)

A3 power curve: 50 jobs x 4 untilted releases = 200 releases (44 workers)
  10/50 power
  20/50 power
  30/50 power
  40/50 power
  50/50 power
⚠⚠ CEILING MODE: 10 targets, each at ITS OWN `eta_tau`, from logs/B_acc_arm2/arm2_chr21_kl_verified.json. The column is `acc_ceiling` and IS the bound.


7,000 trials, each drawing 1 release under Q_p and 4 under Q_q (one per condition) = 35,000 releases
44 workers (0.60 GB/job, 180 GB usable) — CAPPED by --workers 44 (planned 46: CORE-bound: 48 cores less 2 reserved)

  350/7000 trials
  700/7000 trials
  1050/7000 trials
  1400/7000 trials
  1750/7000 trials
  2100/7000 trials
  2450/7000 trials
  2800/7000 trials
  3150/7000 trials
  3500/7000 trials
  3850/7000 trials
  4200/7000 trials
  4550/7000 trials
  4900/7000 trials
  5250/7000 trials
  5600/7000 trials
  5950/7000 trials
  6300/7000 trials
  6650/7000 trials
  7000/7000 trials
methodology: docs/plan.md direction block, 2026-10-24 (git_rev 6882ad55cac661e4c270070ac1ed725e0e99e3da-dirty+7e60fa4b3995)

ARM 2. 10 targets x 7 acc_ceiling x 100 draws = 1000 trials per (acc_ceiling, condition) PER HYPOTHESIS -- matching A1's realised design of 250. chr21, scope chain, beta_share = 1.

⚠⚠ CEILING MODE: each of the 10 targets runs at ITS OWN `eta_tau`, solved from its own `c_quadratic`. The column is `acc_ceiling` and IS the bound -- measured at the solved eta and refined once, not fitted. ⚠ Pinsker is an inequality, so the TRUE ceiling is closer to 0.5 than every column label, and the looseness runs in the direction that flatters us.
      0.50: eta 0.00 .. 0.00 (median 0.00)
      0.55: eta 75.35 .. 86.24 (median 76.41)
      0.60: eta 150.69 .. 172.45 (median 152.81)
      0.65: eta 226.03 .. 258.66 (median 229.21)
      0.70: eta 301.38 .. 344.84 (median 305.59)
      0.75: eta 376.72 .. 431.01 (median 381.98)
      0.80: eta 452.06 .. 517.17 (median 458.36)

### The pre-registered positive controls

| control | requirement | result | passed? |
|---|---|---|---|
| A3 (a) links the UNOBFUSCATED target | near-100% over 40 (target, decoy) pairs | **1.0000** (40/40) | **yes, but see below** |
| A4 finds a block copied VERBATIM from the target | the argmax block is the planted one, over 40 pairs | **1.0000** (40/40) | **yes** |

⚠⚠⚠ **A3's PRE-REGISTERED CONTROL CANNOT FAIL, SO ITS 1.0000 IS NOT A POWER STATEMENT.** Their score sums `-log P_HWE` only where genotypes AGREE, mismatches contribute zero, and the weights are non-negative, so `score(p, q)` is a SUBSET SUM of `score(p, p)` and `Lambda > 0` is algebraically forced. The control verifies the site axis, the genotype encoding and the weights -- the plumbing -- and nothing about power. A gate that cannot fail is decoration, so the power curve below replaces it as the actual power statement.
✅ **A4's control is NOT vacuous and it passed:** a 4,096-position block copied verbatim from the target was the argmax block in 40 of 40 pairs, against an untilted background. A verbatim leak of that size IS detectable by A4.

### A3's power against a KNOWN leak — a labelled DIAGNOSTIC, not an attack outcome

An untilted (`tau` = 0) release with a declared fraction `f` of its positions overwritten by the target VERBATIM. Balanced accuracy against `f`, so each column is A3's power against a leak of exactly that size. Subsets are prefixes of one random permutation, so the curve is monotone by construction.

| condition | f=0 | f=0.0001 | f=0.001 | f=0.01 | f=0.05 | f=0.1 | f=0.25 | f=0.5 | f=1 |
|---|---|---|---|---|---|---|---|---|---|
| D1_extremal | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.9900 [0.9706, 1.0000] | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT |
| D2_nearest | 0.5700 [0.4775, 0.6625] | 0.5700 [0.4775, 0.6625] | 0.5900 [0.4982, 0.6818] | 0.8000 [0.7241, 0.8759] | 0.9900 [0.9706, 1.0000] | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT |
| D3_random | 0.5300 [0.4322, 0.6278] | 0.5300 [0.4322, 0.6278] | 0.5400 [0.4424, 0.6376] | 0.8700 [0.8041, 0.9359] | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT |
| D4_farthest | 0.5600 [0.4640, 0.6560] | 0.5600 [0.4640, 0.6560] | 0.6100 [0.5160, 0.7040] | 0.8200 [0.7456, 0.8944] | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT | 1.0000 EXACT |

- **D1_extremal:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5000 DEGENERATE against a chance line of 0.5000, and carries NO bar, because it carries no information (SE 0.0000, n = 50 per arm).

- **D2_nearest:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5700 [0.4775, 0.6625] against a chance line of 0.5000, whose 95% bar reaches below it (SE 0.0472, n = 50 per arm).

- **D3_random:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5300 [0.4322, 0.6278] against a chance line of 0.5000, whose 95% bar reaches below it (SE 0.0499, n = 50 per arm).

- **D4_farthest:** A3 first reaches 0.75 accuracy at **f = 0.01**; at f = 0 it sits at 0.5600 [0.4640, 0.6560] against a chance line of 0.5000, whose 95% bar reaches below it (SE 0.0490, n = 50 per arm).

### A3 — rarity-weighted genotype linkage, re-posed as 1-of-2

Metrics are rows; `acc_ceiling` is the column. Balanced accuracy, ties split.

⚠ Each cell is `attacker_accuracy` [`acc_ci95`] -- the point estimate and its ERROR BAR, 1.96 x SE either side, clipped to [0, 1]. The bar is what makes an excess over the line checkable rather than asserted; where a cell sits above the line, read whether its lower limit falls below it. ⚠ Two cells carry no bar, because SE is exactly 0 and they are opposites: `DEGENERATE` is a CONSTANT ANSWER at 0.5000 carrying no information, and `EXACT` is PERFECT SEPARATION at 1.0000 carrying all of it.

| metric | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 |
|---|---|---|---|---|---|---|---|
| D1_extremal balanced accuracy * | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE |
| D2_nearest balanced accuracy * | 0.5105 [0.4889, 0.5321] | 0.5060 [0.4845, 0.5275] | 0.5355 [0.5138, 0.5572] | 0.5310 [0.5093, 0.5527] | 0.5505 [0.5289, 0.5721] | 0.5755 [0.5541, 0.5969] | 0.5785 [0.5571, 0.5999] |
| D3_random balanced accuracy * | 0.5220 [0.5001, 0.5439] | 0.5210 [0.4992, 0.5428] | 0.5265 [0.5046, 0.5484] | 0.5145 [0.4926, 0.5364] | 0.5370 [0.5152, 0.5588] | 0.5710 [0.5494, 0.5926] | 0.5610 [0.5393, 0.5827] |
| D4_farthest balanced accuracy * | 0.5005 [0.4788, 0.5222] | 0.5305 [0.5090, 0.5520] | 0.5275 [0.5059, 0.5491] | 0.5405 [0.5189, 0.5621] | 0.5525 [0.5309, 0.5741] | 0.5610 [0.5395, 0.5825] | 0.5680 [0.5466, 0.5894] |
| `p_succ_bound` * | 0.5000 (exact, eta = 0) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) |
| `acc_ceiling` [THE BOUND, measured per target] * | 0.5000 | 0.5500 | 0.6000 | 0.6500 | 0.7000 | 0.7500 | 0.8000 |

### A4 — segment-restricted likelihood ratio, maximised over blocks

⚠ EVERY declared block size is reported. The MAX over blocks is an optimistically biased estimate by construction -- maximising over `B` noisy accuracies -- so the Bonferroni-corrected comparison below is the only claimable one, and the raw max is shown beside it to make the size of that bias visible.

⚠ Two intervals per cell, because the statistic is a MAXIMUM. `acc_ci95` (square brackets) is 1.96 x SE at the argmax block and is ANTI-CONSERVATIVE for a max -- it ignores the selection over B blocks. `(simultaneous)` is 4.57 x SE, Bonferroni over all 20,804 block-cells, and is the interval the bound check below uses. The gap between them IS the selection penalty.

| metric | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 |
|---|---|---|---|---|---|---|---|
| block 512: D1_extremal MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5002 [0.4929, 0.5076] (0.4830, 0.5175) | 0.5020 [0.4953, 0.5087] (0.4864, 0.5176) | 0.5035 [0.4965, 0.5105] (0.4871, 0.5199) | 0.5020 [0.4945, 0.5095] (0.4844, 0.5196) | 0.5210 [0.5005, 0.5415] (0.4731, 0.5689) | 0.5012 [0.4994, 0.5031] (0.4970, 0.5055) |
| block 512: D2_nearest MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5377 [0.5159, 0.5596] (0.4868, 0.5887) | 0.5310 [0.5091, 0.5529] (0.4800, 0.5820) | 0.5292 [0.5075, 0.5510] (0.4784, 0.5801) | 0.5350 [0.5132, 0.5568] (0.4841, 0.5859) | 0.5360 [0.5143, 0.5577] (0.4855, 0.5865) | 0.5360 [0.5142, 0.5578] (0.4851, 0.5869) |
| block 512: D3_random MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5360 [0.5142, 0.5578] (0.4852, 0.5868) | 0.5327 [0.5110, 0.5545] (0.4820, 0.5835) | 0.5270 [0.5052, 0.5488] (0.4760, 0.5780) | 0.5400 [0.5182, 0.5618] (0.4890, 0.5910) | 0.5457 [0.5239, 0.5676] (0.4948, 0.5967) | 0.5453 [0.5235, 0.5670] (0.4945, 0.5960) |
| block 512: D4_farthest MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5337 [0.5120, 0.5555] (0.4829, 0.5846) | 0.5270 [0.5061, 0.5479] (0.4782, 0.5758) | 0.5312 [0.5106, 0.5519] (0.4832, 0.5793) | 0.5310 [0.5093, 0.5527] (0.4804, 0.5816) | 0.5373 [0.5156, 0.5589] (0.4868, 0.5877) | 0.5423 [0.5204, 0.5641] (0.4913, 0.5932) |
| block 4096: D1_extremal MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE |
| block 4096: D2_nearest MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5260 [0.5046, 0.5474] (0.4761, 0.5759) | 0.5240 [0.5028, 0.5452] (0.4745, 0.5735) | 0.5272 [0.5054, 0.5491] (0.4763, 0.5782) | 0.5270 [0.5057, 0.5483] (0.4773, 0.5767) | 0.5370 [0.5152, 0.5588] (0.4862, 0.5878) | 0.5500 [0.5282, 0.5718] (0.4992, 0.6008) |
| block 4096: D3_random MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5210 [0.4992, 0.5428] (0.4701, 0.5719) | 0.5232 [0.5014, 0.5451] (0.4722, 0.5743) | 0.5265 [0.5049, 0.5481] (0.4762, 0.5768) | 0.5305 [0.5089, 0.5521] (0.4801, 0.5809) | 0.5460 [0.5244, 0.5676] (0.4955, 0.5965) | 0.5310 [0.5092, 0.5528] (0.4801, 0.5819) |
| block 4096: D4_farthest MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5220 [0.5011, 0.5429] (0.4732, 0.5708) | 0.5250 [0.5047, 0.5453] (0.4776, 0.5724) | 0.5238 [0.5061, 0.5414] (0.4826, 0.5649) | 0.5300 [0.5086, 0.5514] (0.4800, 0.5800) | 0.5430 [0.5212, 0.5648] (0.4922, 0.5938) | 0.5360 [0.5146, 0.5574] (0.4860, 0.5860) |
| block 32768: D1_extremal MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE | 0.5000 DEGENERATE |
| block 32768: D2_nearest MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5168 [0.4953, 0.5382] (0.4668, 0.5667) | 0.5305 [0.5094, 0.5516] (0.4812, 0.5798) | 0.5265 [0.5053, 0.5477] (0.4770, 0.5760) | 0.5298 [0.5086, 0.5509] (0.4804, 0.5791) | 0.5268 [0.5058, 0.5477] (0.4780, 0.5755) | 0.5457 [0.5244, 0.5671] (0.4959, 0.5956) |
| block 32768: D3_random MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5160 [0.4942, 0.5378] (0.4652, 0.5668) | 0.5270 [0.5052, 0.5488] (0.4760, 0.5780) | 0.5222 [0.5010, 0.5435] (0.4728, 0.5717) | 0.5295 [0.5078, 0.5512] (0.4788, 0.5802) | 0.5360 [0.5149, 0.5571] (0.4869, 0.5851) | 0.5340 [0.5123, 0.5557] (0.4834, 0.5846) |
| block 32768: D4_farthest MAX over blocks | 0.5000 [0.4781, 0.5219] (0.4489, 0.5511) | 0.5185 [0.5029, 0.5341] (0.4820, 0.5550) | 0.5100 [0.4942, 0.5258] (0.4731, 0.5469) | 0.5138 [0.4978, 0.5297] (0.4766, 0.5509) | 0.5165 [0.5006, 0.5324] (0.4795, 0.5535) | 0.5292 [0.5133, 0.5452] (0.4920, 0.5665) | 0.5288 [0.5131, 0.5444] (0.4922, 0.5653) |
| `p_succ_bound` * | 0.5000 (exact, eta = 0) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) | ~1.0 (vacuous) |
| `acc_ceiling` [THE BOUND, measured per target] * | 0.5000 | 0.5500 | 0.6000 | 0.6500 | 0.7000 | 0.7500 | 0.8000 |

**THE BOUND CHECK.** One-sided on the excess over `acc_ceiling`, Bonferroni-corrected across the cells actually examined: A3 over 28 cells (z_crit = 2.91), A4 over 20,804 block-cells (z_crit = 4.57).

⚠⚠ THE LINE IS THE COLUMN LABEL, not `p_succ_bound`. In this mode each target ran at its own `eta_tau` solved to land on that label, so the label IS the measured `acc_ceiling` and IS the bound. Handing the label to `p_succ_bound` instead would test against (1 + label)/2 -- 0.7500 at the 0.50 column up to 0.9000 at 0.80 -- which is the defect fixed on 2026-10-29.

⚠ 2 cell(s) exceed `acc_ceiling` by POINT ESTIMATE, each with a 95% bar REACHING BELOW the line. Reported, not claimed:
    A3 D2_nearest at acc_ceiling=0.5: 0.5105 [0.4889, 0.5321] vs bound 0.5000 (+0.0105, z = +0.95)
    A3 D4_farthest at acc_ceiling=0.5: 0.5005 [0.4788, 0.5222] vs bound 0.5000 (+0.0005, z = +0.05)

⚠⚠ 1 cell(s) exceed `acc_ceiling` and their 95% bar does NOT reach below the line. Not significant after correction, so not claimable -- but stated here rather than absorbed, because this is the case an uncorrected reader would call a violation:
    A3 D3_random at acc_ceiling=0.5: 0.5220 [0.5001, 0.5439] vs bound 0.5000 (+0.0220, z = +1.97)

⚠ MULTIPLICITY. 112 cells were examined. Of these, 28 are A3 SINGLE estimates, where a one-sided 5% level puts about 1.4 above the line BY CHANCE even if the mechanism is exact. The other 84 are A4 MAXIMA over hundreds of blocks each (20,804 block-cells in total), whose chance of landing above the line is far higher than 5% and is NOT quoted here -- that is precisely what the z_crit of 4.57 corrects for. A bare count of 95% exceedances is not evidence in either attack; the corrected z is.

no condition exceeded `acc_ceiling` significantly, in either attack

✅ Both positive controls pass, so the nulls are lower bounds on distinguishability rather than broken instruments. ⚠ They remain LOWER bounds: the upper bound is `acc_ceiling`, which is proof-derived.
⚠ Minimum detectable advantage at 1000 trials per arm per cell (100 draws x 10 targets) is 0.0278 UNCORRECTED. A cell resolves only where the bound exceeds 0.5278, so the 0.50 column(s) do NOT resolve below it in either attack; those rows are reported, not claimed, and the PROOF carries them. A4's corrected threshold is 2.78x wider still.
⚠ EVERY decoy pool is the other nine prepared targets, not the 3,202-sample panel, and it cuts OPPOSITE ways for the two real conditions: widening it would find a NEARER neighbour for D2, making that pair less distinguishable and the attacker weaker, and a FARTHER genome for D4. So D4 as run is a LOWER BOUND on the real worst case and is the limitation to report; the narrow D2 is conservative, not flattering.
⚠ D2_nearest is the LEAST distinguishable real pair and so the WEAKEST test of the bound here, not the hardest. 'Nearest is hardest' is a LINKAGE intuition and belongs to Arm 1.
⚠ The 4 conditions share each release and must not be pooled.
⚠ A3 scores on 302,547 bi-allelic cohort sites of 332,521 chain positions; the rest carry weight the attack cannot see, which flatters us.

JSON -> logs/B_acc_arm2/a34_ceil_n1000_d4.json
