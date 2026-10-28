READ MAPPING ON THE RELEASED GRAPH. ARM 1, chr21, target HG00438, 5 external donors, median over donors. Metrics are ROWS, the knob is the COLUMN. vg 1.68.0, single-end, `vg autoindex --workflow giraffe`.
⚠ Each metric is a percentage of ALIGNED reads, per `external/PanMixer/plot_util/plot_helpers.py:15-19` — not of reads submitted. ⚠ UTILITY, not privacy. ⚠ Levels are NOT comparable to their published figures (vg 1.68 vs 1.36); the DELTA is.

| metric | baseline | 0 | 250 | 500 | 1262 | 2821 | 5642 |
|---|---|---|---|---|---|---|---|
| `perfect` % | 76.451 | 76.452 | 76.453 | 76.456 | 76.448 | 76.453 | 76.455 |
| `gapless` % | 95.947 | 95.953 | 95.952 | 95.954 | 95.950 | 95.950 | 95.953 |
| `mapq60` % | 79.191 | 79.192 | 79.190 | 79.128 | 79.189 | 79.191 | 79.190 |
| reads aligned (median) | 11,185,151 | 11,185,214 | 11,185,267 | 11,185,242 | 11,185,178 | 11,185,159 | 11,185,156 |

**THE DELTA AGAINST OUR OWN BASELINE** — the only comparable quantity.

| metric | 0 | 250 | 500 | 1262 | 2821 | 5642 | worst | PanMixer's own delta |
|---|---|---|---|---|---|---|---|---|
| `perfect` delta, pts | +0.001 | +0.002 | +0.004 | -0.003 | +0.002 | +0.004 | **-0.003** | -0.01 |
| `gapless` delta, pts | +0.006 | +0.005 | +0.007 | +0.003 | +0.003 | +0.006 | **+0.003** | -0.01 |
| `mapq60` delta, pts | +0.000 | -0.001 | -0.063 | -0.002 | -0.001 | -0.001 | **-0.063** | +0.00 |

worst degradation anywhere on the grid: **-0.063 points**; largest movement in either direction: **-0.063**.
⚠ A positive delta is the obfuscated graph aligning BETTER than the baseline, which at this magnitude is not a gain.
⚠ THE ONE MOVEMENT ABOVE NOISE IS GRAPH-LEVEL, NOT DONOR-LEVEL, and the distinction matters. At `eta_tau` = 500 ALL FIVE donors drop together (-0.038 to -0.063 pts each), so that graph really does align slightly worse for everyone — it is not one noisy donor. But it is NOT monotone in the knob (1262, 2821 and 5642 return to -0.001/-0.002), so it is draw-to-draw variation in WHICH variants the released row happens to carry, not an effect of `eta_tau`. One release per point is the limit here: with n = 1 graph per column this spread cannot be separated from sampling, and quoting -0.063 as the knob's cost would overstate it.
⚠ One target (HG00438) substituted into a 44-member graph, so the released row is 1 of 44: this measures the cost of obfuscating ONE member, exactly as their Fig 5A does, and must not be read as the cost of obfuscating all.
⚠ `e_differ`-scale context: the substituted row differs from the original in ~137,000-139,000 of 681,648 cells at every eta, so the graph really is being perturbed; a flat result is NOT a null edit.

JSON -> logs/B_readmap_eta/readmap.json
