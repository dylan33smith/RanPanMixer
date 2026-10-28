# Results — the `acc_ceiling` 0.50–0.80 regime, chr21, Arm 1 against PanMixer

Every table here is **emitted by a script, not typed**, and carries its own provenance
line. This directory is the curated set the paper cites; the raw run logs it was taken
from live under `logs/`, which is gitignored on every branch and is therefore not in git
on any of them. That is why these are copies rather than a merge of another branch.

⚠ **THIS IS A SNAPSHOT.** Re-running a script rewrites its file under `logs/` and does
NOT update the copy here. Each file states the revision it came from; if a number here
disagrees with a fresh run, the fresh run wins and this directory needs re-copying.

## The regime

Arm 1, **chr21 only**, K = 86 per-target leave-one-out, 44 HPRC individuals in the same
order as PanMixer, reported on an `acc_ceiling` axis from 0.50 to 0.80. Each target
carries its OWN fitted `eta_tau` (the per-target `c` spreads 35.7%, so one shared `eta`
could not place every individual on the same ceiling) and its OWN chr21 `beta_share`
(~0.0150249, NOT 1.0).

⚠ `acc_ceiling` is a **Pinsker bound**, so the true ceiling is closer to 0.5 than every
figure here. ⚠ And it is **per pair**: a directory or column labelled 0.80 is the
`D1_extremal` ceiling alone, and the same cell's `D2_nearest` ceiling is 0.678293.
⚠ `acc_ceiling` 0.50 is `eta` = 0 — no protection, zero target information — so its
level is a floor, never a privacy result.

## Files

| file | what it is | emitted by |
|---|---|---|
| `headtohead.md` | the paired alteration-rate comparison, ours at 0.80 against PanMixer at saturation | `scripts/arm1_headtohead.py` |
| `matched_fidelity.md` | the capacity at which a privacy comparison becomes legitimate under constraint 5 | `scripts/matched_fidelity_capacity.py` |
| `pm_draw_noise.md` | PanMixer's own draw noise, 5 seeds per cell — the analysis that changed a conclusion | `scripts/pm_draw_noise.py` |
| `af_ld.md` | `af_loss` and `ld_loss` for BOTH tools on ONE build, their evaluators unmodified | `scripts/arm1_af_ld_results.py` |
| `af_ld_sensitivity.md` | **read this before quoting `af_ld.md`** — both its headlines are policy-dependent | `scripts/arm1_afld_sensitivity.py` |
| `readmap.md` + `.json` | read mapping on the released graph, PanMixer's Fig 5A quantity | `scripts/readmap_results.py` |
| `summary.html` | the figure page: frontier chart and 44-individual paired dot plot, drawn to scale | generated from `summary_figures.json` |
| `summary_figures.json` | every figure on that page, so no number on it is retyped by hand | consolidated |
| `kl_fit_verified.json` | the per-target fitted `eta` at each ceiling label, and its measured error | `scripts/arm1_chr21_kl.py --verify` |
| `arm1_release_record.json` | the release record for all 308 cells (44 targets x 7 ceilings) | `scripts/arm1_release.py --ceilings` |
| `arm1_exact_endpoint.json` | the exact `E[differ]` and `E_Q[u]` endpoint, zero Monte-Carlo error | same run, rides along |
| `arm1_attack.md` + `.json` | **ARM 1: A1 and A3 at power on this grid** | `scripts/arm1_attack.py` |
| `arm2_attack_a1.md` + `.json` | **ARM 2: A1 on this grid**, chr21 at share 1.0, K = 88, nothing held out | `scripts/arm2_attack.py` |

## The headline numbers, each with the caveat that bounds it

- **33 of 44 individuals are more faithful under our release** at `acc_ceiling` 0.80
  against PanMixer at saturation — sign test p = 0.00126, paired z = -5.25.
  ⚠ Ours is an EXACT expectation and theirs a 5-draw mean; residual draw noise on their
  side is ~29% of the effect. ⚠ On ONE draw this was 27/44 with the two tests
  disagreeing; `pm_draw_noise.md` is why it moved, and the median barely changed — the
  dispersion fell.
- **For 11 of 44, PanMixer is more faithful at every setting it offers.** A quarter of
  the cohort with no matched point. The same 11 are exactly those the paired test says we
  lose, by two independent routes.
- **A matched-fidelity capacity exists for 33 of 44**, median 0.2524. ⚠ A crossing is not
  a privacy result; it is the capacity at which one becomes legitimate.
- **The grid is solved to 7.224e-08** worst `|acc_ceiling - label|` over 44 targets x 6
  points after one refinement.
- **`af_loss`: they distort allele frequencies less than we do** — but
  `af_ld_sensitivity.md` shows that verdict runs from 0/44 against us to **31/44 for us**
  across three defensible policies for uncalled cells, with nothing about either
  mechanism changing. The cause is structural: our release is a fixed-support path while
  theirs edits a copy of the target's own row.
- **`ld_loss`: we distort linkage less**, and the honest size is **1.30x, not the 7.53x**
  the shipped code reports — whose largest 6 pairs of ~7,500 carry 94% of their total.
- **The optimal attacker never beats its own pair ceiling.** A1 is the Neyman-Pearson
  optimum for its pair, so no attack in any literature can beat it there. Over 44 targets
  x 7 labels x 3 conditions x 10 draws (n = 440 per arm per row), **0 of 21 pooled rows
  exceed**, and the margin WIDENS: 0.95 SE below the bound at 0.55, **6.07 SE below at
  0.80**, where it reaches 0.707 against a permitted 0.800. ⚠ Pinsker is conservative so
  the true margin is smaller than it looks; the readable thing is the direction.
  ⚠ `acc_ceiling` is PER PAIR here and printed as a row -- the 0.80 label is D1's ceiling
  alone, and D2's at that label is 0.689971 (mean over targets).
  ⚠ A3's `D1_extremal` row is DEGENERATE at all 7 labels (TPR 1.0000 / TNR 0.0000),
  exactly as predicted in writing before the run, and A3 sits z -4.99..-8.07 below the
  bound on D2/D3. A3 is a weak attacker; only A1's null is evidence.
- **ARM 2's optimal attacker lands where ARM 1's did, and no condition significantly
  exceeds its ceiling in either arm.** At the 0.80 label A1 reaches **0.7080** in Arm 2
  against **0.706818** in Arm 1 -- within 0.0012, against the same 0.800 bound.
  ⚠⚠ THIS IS NOT A REPRODUCTION. The arms release DIFFERENT objects: Arm 1 is chr21
  carrying the genome-normalised `beta_share` 0.0150249 under K = 86 leave-one-out; Arm 2
  is chr21 at share **1.0** under the full K = 88 with nothing held out. What agrees is
  where the Neyman-Pearson optimum lands against a bound the two were matched on by
  `acc_ceiling`, which is what matching on `acc_ceiling` rather than `eta_tau` buys.
  ⚠ Arm 2 is n = 250 per cell per arm against Arm 1's 440, and its `delta` is 0.0556 --
  a cell resolves only where the ceiling exceeds 0.5556, so the 0.55 and 0.60 columns add
  nothing over the ceiling itself. One point-estimate exceedance (D3_random at 0.55,
  0.5580 vs 0.5500, z = +0.36) sits inside exactly that band.
  ⚠⚠ **A SCOPE SPLIT THAT MUST BE STAMPED, NEVER POOLED.** Arm 2's ATTACK evidence here is
  chr21 at share 1.0 -- the only scope its driver implements -- while Arm 2's UTILITY
  numbers are genome-wide. Two different released objects; they may not share a table.
- **Read mapping is flat**: worst degradation anywhere is -0.063 points, against
  PanMixer's own published -0.01. ⚠ n = 1 graph per column; absolute levels are not
  comparable to their published figures (vg 1.68 against 1.36), only deltas.

## What is NOT here

- **PanMixer's own privacy evaluators** (`gap_score`, `diploid_gap_score`,
  `MIA_privacy`) on this grid -- DEFERRED by decision, 2026-10-28. ⚠ And an audit
  established what a fresh run could conclude: all three are defeated by our `eta` = 0
  release too, so they measure the release FORM (a complete mosaic of other people's
  haplotypes) and NOT the tilt. They are needed for the head-to-head parity claim --
  theirs is identified at rank 1 of 44 at `capacity` 0.01/0.05/0.1 where ours never is --
  but never as evidence that the tilt protects anyone.
- **Arm 2's A3/A4, and Beagle's reconstruction attack** on this grid -- both RUNNING as
  of this snapshot, not yet here. Arm 2's A3/A4 is on the same seven etas as its A1 above;
  Beagle is 308 rows on OUR 340,849 build (PanMixer's own protocol pass).
  ⚠⚠ **BEAGLE'S EXCLUDE-SELF PASS IS DROPPED BY DECISION (2026-10-28), and that costs a
  control.** PanMixer's published protocol passes no `excludesamples=`, and **39 of our 44
  targets ARE in Beagle's 1000G reference panel** -- so for those 39 the panel holds the
  target's own true row and Beagle can copy from it. The exclude-self pass is what
  separates what the RELEASE leaks from what the PANEL leaks; without it a reconstruction
  figure cannot be attributed to the release alone. The 44 per-target panels exist on
  disk, so it is a ~2.2 h addition whenever wanted. ⚠ The head-to-head against PanMixer is
  unaffected: their published figures use the same no-exclusion protocol, so the
  cross-tool comparison is like-for-like. It is our own column that loses a control.
- **Read mapping** on this grid. The existing numbers are on the old genome-wide `eta`
  grid and are deliberately NOT included here.
- **PanMixer's privacy evaluators** (`gap_score`, `diploid_gap_score`, `MIA_privacy`) on
  this grid. They ran on the earlier genome-wide `eta` grid, where they were defeated by
  our `eta` = 0 release too — which makes them a statement about release FORM rather than
  about the tilt.
- **Arm 2**, which carries genome-wide alone and is not part of this comparison.
- **Anything at matched privacy between the two tools.** `capacity` is a utility-loss
  budget and bounds nothing about distinguishability; `acc_ceiling` does not transfer to
  their mechanism. Every cross-tool row here pairs the two at their OWN endpoints, and
  putting both on one build removed the build confound while leaving that one untouched.
