PanMixer's OWN attacks, run on OUR ARM 1 releases. 44 targets x 7 `tau` = 308 releases, chr21, K = 86 release on their K = 88 axis.
Metrics are ROWS, `tau` is the COLUMN. Median over 44 targets unless stated.
⚠ PRIMARY subset `snp_only` -- the paper's 'SNPs overlapping with 1000 Genomes'.
⚠ Truth axis `/home/ds85/pm_ourbuild/starting_data` (340,849 sites), asserted against the manifest.
⚠⚠ THIS ROOT WAS RUN FOR `beagle` ONLY. Every row of a non-run evaluator reads `n/a (NOT RUN on this root)` -- that is an absence of measurement, NOT a null result, and must never be read as one.
⚠ Provenance `logs/B_arm1_attacks.md`; published values `/data/ds85/RanPanMixer/runs/beagle_attack/published/PROVENANCE.md`.

| metric | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 | PanMixer PUBLISHED |
|---|---|---|---|---|---|---|---|---|
| **BEAGLE RECONSTRUCTION — paper-faithful (target EXCLUDED from panel)** | | | | | | | | |
| (exclude-self runs not yet present) | n/a | n/a | n/a | n/a | n/a | n/a | n/a | |
| **BEAGLE RECONSTRUCTION — their CODE's protocol (no exclusion)** | | | | | | | | |
| `beagle_accuracy` − `major_allele_accuracy` | -0.046384 | -0.044739 | -0.045508 | -0.046010 | -0.045470 | -0.044975 | -0.044885 | their code passes no `excludesamples` |
| `reconstruction_accuracy` | 0.232482 | 0.228238 | 0.232528 | 0.227460 | 0.227124 | 0.225460 | 0.227292 | no number published |
| `obfuscated_accuracy` [the attacker's floor] | 0.826582 | 0.828554 | 0.830391 | 0.828874 | 0.828950 | 0.829078 | 0.829854 | no number published |
| **WHERE WE SIT ON THEIR FIG S2 x-AXIS** | | | | | | | | |
| cells altered, `snp_only` [their '# of Obfuscated Variants'] | 81,486 | 82,197 | 81,510 | 82,376 | 81,664 | 81,588 | 81,304 | their curve spans 1,657 to 232,613 |
| sites altered, `snp_only` | 68,254 | 68,608 | 68,383 | 68,775 | 68,496 | 68,343 | 68,050 |  |
| their DIGITIZED curve at our cell count [interpolated] | -0.149148 | -0.149592 | -0.149163 | -0.149704 | -0.149259 | -0.149211 | -0.149035 | **OUR digitization of their Fig S2** |
| **GAP SCORE — rarity-weighted GENOTYPE linkage** | | | | | | | | |
| `genotype_g_to_gstar` | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | no box statistic printed |
| targets LINKED (Gap Score >= 0) of 44 * | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | **55 of 100 Can Link** across their eps grid |
| **DIPLOID GAP SCORE — the HAPLOTYPE version** | | | | | | | | |
| `haplotype_g_to_gstar` | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | ablation, no numbers published |
| targets LINKED (>= 0) of 44 * | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) |  |
| **MEMBERSHIP INFERENCE** | | | | | | | | |
| `hap_original_rank` of 44 * | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | **no number published**, either document |
| targets at rank 1 of 44 * | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | 'the red line falls within the bulk' |
| `hap_bin_match_frac_vs_original` | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | n/a (NOT RUN on this root) | no number published |

**THE PANEL-MEMBERSHIP CONTROL.** 39 of 44 targets are IN the 1000G reference panel and 5 are not (HG01123, HG02109, HG02486, HG02559, NA21309). Their code's protocol lets Beagle copy an in-panel target's own true row out of the reference; the out-of-panel targets cannot be copied at all, so the gap between the two groups is how much of the 'reconstruction' comes from the PANEL rather than from our release.

| group | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 |
|---|---|---|---|---|---|---|---|
| `reconstruction_accuracy`, IN panel (n=39) | 0.230998 | 0.231867 | 0.236556 | 0.229170 | 0.227609 | 0.225431 | 0.227547 |
| `reconstruction_accuracy`, OUT of panel (n=5) | 0.234503 | 0.216547 | 0.225243 | 0.218863 | 0.222125 | 0.232373 | 0.225228 |

⚠ `sampled_gap_score` is NOT runnable here: it reads a PRE-stacker object from `starting_data/chr21/subjects/`, which is absent from the shipped tree. THREE of their four privacy evaluators ran, never four.
⚠ `wged` and `pmi` are NOT Beagle statistics in their paper. 'WGED' appears once, as *Weighted Graph Edit Distance*, a conceptual foil; PMI is the mechanism's own privacy metric. Their `tools/beagle/beagle_stats.py`, which computes both, is DEAD CODE that cannot run as shipped. The published statistic is an accuracy DIFFERENCE.
⚠ NOT a same-build comparison (decision F): their published numbers come from their own vg-1.36 pipeline and ours from our vg-1.68 build. Dylan's instruction for this run was to compare anyway and name the gap, which is what this line does.
⚠ Gap Score uses `>=` per Eq 6 and Fig 3A. The published Results text says '>0'. ⚠ The PDF's TEXT LAYER returns '>' for the Fig 3A glyph, which renders as '>=' -- anyone grepping that file for the rule gets the wrong answer.
⚠ `eps_private` is 0.002 published and **0.001** in the preprint, with the definition flipping from 'largest' to 'smallest'. No comparison here depends on it.

JSON -> logs/B_arm1_beagle_acc/beagle_acc.json
