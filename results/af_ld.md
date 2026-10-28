ARM 1 IN-GRAPH UTILITY -- EXPLORATORY, NOT PRE-REGISTERED (decision H).
⚠ `af_loss`/`ld_loss` are UTILITY metrics in PanMixer's framing, NOT privacy.
chr21. build: OURS -- vg 1.68 GRCh38 cohort44, 340,849 records
Ours: 44 targets x 7 acc_ceiling, median over targets.
Theirs: 44 subjects x 8 capacity, median over subjects.
⚠ ld_loss window is 5 POSITION UNITS, not 5 kb (`ld_loss.py:107`).
⚠ acc_ceiling 0.50 is eta = 0 -- NO protection; its level is not a privacy result.

**OURS** (columns: `acc_ceiling`; n = 44 targets per cell)

| variable | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 |
|---|---|---|---|---|---|---|---|
| `af_loss` [all] | 0.003600 | 0.003604 | 0.003577 | 0.003592 | 0.003591 | 0.003596 | 0.003580 |
| `af_loss` [complex] | 0.004283 | 0.004276 | 0.004300 | 0.004269 | 0.004290 | 0.004299 | 0.004302 |
| `af_loss` [snp_only] | 0.003336 | 0.003315 | 0.003310 | 0.003328 | 0.003311 | 0.003321 | 0.003305 |
| `af_loss` [snp_only_maf_0.01] | 0.000968 | 0.000978 | 0.000962 | 0.000973 | 0.000976 | 0.000974 | 0.000967 |
| `af_loss` [snp_only_maf_0.05] | 0.002171 | 0.002167 | 0.002191 | 0.002183 | 0.002169 | 0.002190 | 0.002192 |
| `af_loss` [snp_only_maf_0.1] | 0.003624 | 0.003720 | 0.003698 | 0.003699 | 0.003671 | 0.003696 | 0.003709 |
| `af_loss` [snp_only_maf_0.5] | 0.006964 | 0.006898 | 0.006960 | 0.006929 | 0.006931 | 0.006943 | 0.006927 |
| `ld_loss` per pair | 0.005426 | 0.005571 | 0.005392 | 0.005310 | 0.005436 | 0.005291 | 0.005405 |
| `wd_af` = 2x [all] -- AF [their Fig 4A at `eps_private`: 0.006] | 0.007200 | 0.007209 | 0.007154 | 0.007185 | 0.007182 | 0.007192 | 0.007161 |
| `wd_af_snp_only` = 2x [snp_only] -- AF (SNPs only) [their Fig 4A at `eps_private`: 0.006] | 0.006673 | 0.006630 | 0.006620 | 0.006656 | 0.006622 | 0.006643 | 0.006609 |
| `wd_af_complex` = 2x [complex] -- AF (complex) | 0.008566 | 0.008552 | 0.008599 | 0.008538 | 0.008580 | 0.008597 | 0.008605 |
| `wd_af_snps_0_01` = 2x [snp_only_maf_0.01] -- AF (MAF < 0.01) | 0.001936 | 0.001956 | 0.001923 | 0.001947 | 0.001951 | 0.001948 | 0.001935 |
| `wd_af_snps_0_05` = 2x [snp_only_maf_0.05] -- AF (MAF in [0.01, 0.05)) [their Fig 4A at `eps_private`: 0.005] | 0.004342 | 0.004334 | 0.004383 | 0.004366 | 0.004337 | 0.004381 | 0.004384 |
| `wd_af_snps_0_1` = 2x [snp_only_maf_0.1] -- AF (MAF in [0.05, 0.1)) | 0.007248 | 0.007440 | 0.007396 | 0.007397 | 0.007342 | 0.007392 | 0.007418 |
| `wd_af_snps_0_5` = 2x [snp_only_maf_0.5] -- AF (MAF in [0.1, 0.5]) | 0.013929 | 0.013795 | 0.013920 | 0.013858 | 0.013861 | 0.013887 | 0.013855 |

**THEIRS** (columns: `capacity`; n = 44 per cell)

| variable | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| `af_loss` [all] | 0.000216 | 0.000898 | 0.001730 | 0.002779 | 0.002840 | 0.002841 | 0.002841 | 0.002841 |
| `af_loss` [complex] | 0.000570 | 0.001393 | 0.001864 | 0.002591 | 0.002643 | 0.002644 | 0.002644 | 0.002644 |
| `af_loss` [snp_only] | 0.000029 | 0.000670 | 0.001644 | 0.002940 | 0.002991 | 0.002991 | 0.002991 | 0.002991 |
| `af_loss` [snp_only_maf_0.01] | 0.000083 | 0.000231 | 0.000453 | 0.000670 | 0.000691 | 0.000691 | 0.000691 | 0.000691 |
| `af_loss` [snp_only_maf_0.05] | 0.000001 | 0.000504 | 0.001058 | 0.001713 | 0.001776 | 0.001776 | 0.001776 | 0.001776 |
| `af_loss` [snp_only_maf_0.1] | 0.000001 | 0.000749 | 0.001923 | 0.003415 | 0.003452 | 0.003453 | 0.003453 | 0.003453 |
| `af_loss` [snp_only_maf_0.5] | 0.000000 | 0.001222 | 0.003441 | 0.006502 | 0.006625 | 0.006625 | 0.006625 | 0.006625 |
| `ld_loss` per pair | 0.000077 | 0.037462 | 0.039025 | 0.040408 | 0.040671 | 0.040671 | 0.040671 | 0.040671 |
| `wd_af` = 2x [all] -- AF [their Fig 4A at `eps_private`: 0.006] | 0.000431 | 0.001797 | 0.003460 | 0.005558 | 0.005680 | 0.005682 | 0.005682 | 0.005682 |
| `wd_af_snp_only` = 2x [snp_only] -- AF (SNPs only) [their Fig 4A at `eps_private`: 0.006] | 0.000059 | 0.001340 | 0.003289 | 0.005881 | 0.005982 | 0.005982 | 0.005982 | 0.005982 |
| `wd_af_complex` = 2x [complex] -- AF (complex) | 0.001140 | 0.002786 | 0.003729 | 0.005183 | 0.005285 | 0.005288 | 0.005288 | 0.005288 |
| `wd_af_snps_0_01` = 2x [snp_only_maf_0.01] -- AF (MAF < 0.01) | 0.000165 | 0.000462 | 0.000906 | 0.001341 | 0.001381 | 0.001382 | 0.001382 | 0.001382 |
| `wd_af_snps_0_05` = 2x [snp_only_maf_0.05] -- AF (MAF in [0.01, 0.05)) [their Fig 4A at `eps_private`: 0.005] | 0.000001 | 0.001009 | 0.002116 | 0.003426 | 0.003552 | 0.003552 | 0.003552 | 0.003552 |
| `wd_af_snps_0_1` = 2x [snp_only_maf_0.1] -- AF (MAF in [0.05, 0.1)) | 0.000002 | 0.001498 | 0.003847 | 0.006829 | 0.006905 | 0.006907 | 0.006907 | 0.006907 |
| `wd_af_snps_0_5` = 2x [snp_only_maf_0.5] -- AF (MAF in [0.1, 0.5]) | 0.000001 | 0.002445 | 0.006882 | 0.013004 | 0.013250 | 0.013250 | 0.013250 | 0.013250 |

af_loss[all] per-target direction, acc_ceiling 0.5 -> 0.8: 25/44 down, 19/44 up, median -8.352e-06
  -> NO systematic direction, unlike target_fidelity and e_q_u (44/44).

**PAIRED BY INDIVIDUAL, ONE BUILD** -- our `acc_ceiling` 0.8 against their `capacity` 0.9, same 44 people, same `pangenome.npy`.
⚠⚠ CROSS-OPERATING-POINT: their `capacity` is not a privacy level in our sense and `acc_ceiling` does not transfer to their mechanism, so this pairs the two tools at THEIR OWN endpoints, not at matched privacy (standing constraint 5). It answers 'at each tool's most-perturbing setting, which disturbs the graph less', nothing more.

| metric | ours (median) | theirs (median) | we win | sign p | paired z |
|---|---|---|---|---|---|
| `af_loss` [all] | 0.003580 | 0.002841 | **0/44** | 1.14e-13 | +25.28 |
| `ld_loss` per pair | 0.005405 | 0.040671 | **40/44** | 1.71e-08 | -19.65 |
⚠ 'we win' counts individuals where OUR loss is lower. n is the paired count after dropping exact ties.

**ON THEIR PUBLISHED SCALE (`wd_af`), ours against their `eps_private` point.** ⚠ Cross-operating-point, and the published figures are from THEIR vg-1.36 build -- see the THEIR_FIG4A note in this script. A ratio below 1.00 means WE distort less.

| panel | ours (median, worst column) | theirs published | ratio |
|---|---|---|---|
| AF | 0.007209 | 0.006 | **1.20x** |
| AF (SNPs only) | 0.006673 | 0.006 | **1.11x** |
| AF (MAF in [0.01, 0.05)) | 0.004384 | 0.005 | **0.88x** |
⚠ 'worst column' is the largest median across the grid, so this is the least flattering reading of our own sweep.

af_loss[all] grid range 0.003577..0.003604 (span 2.76e-05, 0.77% of the floor); theirs 0.000216..0.002841
ld_loss grid range 0.005291..0.005571; theirs 0.000077..0.040671
