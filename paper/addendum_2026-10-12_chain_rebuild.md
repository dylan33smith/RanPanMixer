# Addendum 2 to `paper/prereg_2026-10-10.md` — the chain rebuild's effect on two quoted figures

**Status: SIGNED 2026-10-12.** Authorised by Dylan Smith on 2026-10-12, selecting
"Dated addendum, as with decision K" when asked how to handle the stale figures.
`git_rev` a925fd4 — the code state this addendum attests to.

⚠ **WHY THIS IS A SEPARATE FILE AND THE PRE-REGISTRATION IS UNTOUCHED.** Standing
constraint 6 freezes that file: "If any of them must change, it is a NEW dated
pre-registration file and this one stays in the repository — superseded, never edited."
Editing a signed document after the fact is precisely what a signature is supposed to
make impossible, so nothing in `prereg_2026-10-10.md` has been altered.

⚠ **IT IS DELIBERATELY NOT NAMED `prereg_*.md`.** `src/attacks.py:require_signed_prereg`
globs that pattern (`src/attacks.py:218`) and returns the first signed, clean match, so a
second file matching it could become the document an experiment stamps itself with. An
addendum must never be able to impersonate the pre-registration.

⚠⚠ **THIS ADDENDUM CHANGES NO PRE-REGISTERED FIELD.** `tau`, the utility FUNCTION, the
endpoints, the attacker set, the attack database and the decoy design are all untouched,
and no outcome was inspected before writing it. It corrects two DESCRIPTIVE numbers that
the signed file quotes in passing, and which were measured on a chain that no longer
exists. Unlike addendum 1, nothing here was forced by a crash and nothing is a decision —
it is a measurement correction.

---

## 1. What changed

Nothing in the design. The **chain itself** was rebuilt between the signature and now, by
decision C (the declared-rate fill, which replaced the coordinate cut) and decision K
(untraversable zero-distance blocks excluded, addendum 1). chr21's chain went from
**T = 305,886** positions to **T = 332,521**.

Two figures in the signed file were measured on the old chain and are now wrong.

## 2. The figures

**Section 7, R1 (`prereg_2026-10-10.md:94-95`)** reads:

> R1 observes that `beta_v` is not normalised over readable positions, so ~18.5% of the
> budget is inert and `u` cannot exceed ~0.815.

⚠ **Both numbers are superseded.** Measured on the live chain, 2026-10-12, HG00096,
chr21, `beta_scope` = `chain`, **ARM 2**:

| quantity | signed file (T = 305,886) | **live (T = 332,521)** |
|---|---|---|
| share of `beta_v` the target cannot call | ~18.5% | **30.749104%** |
| `u_ceiling` | ~0.815 (0.815455) | **0.692509** |
| positions the target calls | 253,846 of 305,886 = 0.829871 | 258,208 of 332,521 = **0.776516** |

**The inert share very nearly doubled.**

## 3. Why it is a build change and not a miscalculation

The rebuild added **26,635** chain positions, and only **4,362** of them are called for
HG00096 — the added positions are predominantly in map-silent and interior-plateau
regions, which is exactly where a single sample's calls are sparsest. That drops the
called fraction from 0.829871 to 0.776516, and the inverse-support weighting pushes the
weighted ceiling lower still than the raw count, because rare variants are both
up-weighted and more often uncalled.

⚠ **The decisive check:** restricting the live arrays to the old coordinate span
reproduces **T = 305,886 and `u_ceiling` = 0.815455 exactly.** The arithmetic is
unchanged; only the axis moved. Had the two disagreed under the old span, this would be a
bug report rather than an addendum.

This was found because a new experiment's first row failed to reproduce a recorded number,
and it was *diagnosable* only because the record carried its **denominator** (253,846 of
305,886) rather than the ratio alone.

## 4. What it affects

**(a) R1's magnitude, not its disposition.** The signed file *declines* R1's fix
(`:92-104`) on deployment-realism grounds — normalising over `readable` would make the
normalising set a function of the target's own VCF, which is the standing-constraint-1
violation the scope choice exists to avoid. **That decision stands and is unchanged.**
What changes is the size of the thing being accepted: the accepted cost is ~30.7% of the
budget, not ~18.5%.

**(b) `u_ceiling` as a denominator.** Anything divided by 0.815455 on this build
overstates by a factor of 1.178. Every live utility figure must be read against
**0.692509** for HG00096 chr21 at `beta_share` = 1.

**(c) It is per target, and Arm 1 has 44 of them.** `u_ceiling` is a readout of one
target's missingness pattern, so 0.692509 belongs to HG00096 alone. Arm 1's 44
leave-one-out targets each carry their own ceiling, none of them measured yet, and
constraint 8 forbids pooling them.

**(d) A newly identified inert class, logged here for completeness.** Separately from the
target's missingness, **4,448 live chr21 positions (1.3377%) are monomorphic among called
donors** and carry **3.520479%** of chr21's `beta_v`. Every allowed state gives the same
allele there, so the potential is a positive constant that cancels in the normalised
forward and backward weights: that weight provably cannot move the release. Only
**0.446875%** of the budget is *new* inert mass beyond the target-uncalled share already
counted above. Unlike R1's `readable` fix this test reads only the donor matrix and is
target-independent, so reclaiming it would not engage constraint 1 — but it is far too
small to matter and is **not** adopted here.

## 5. What it does NOT affect — stated explicitly

**No privacy claim moves.** `delta_u`, the global range of the utility and the only
property the bound depends on, is **1 by construction** on every release on record,
because `sum_v beta_v` = the declared share and the shares sum to 1. `u_ceiling` is a
utility DENOMINATOR, never an input to the calibration. So `tau`, `eta_tau`,
`epsilon_tau`, `tv_bound` and `p_succ_bound` are all unchanged, and the Arm 2 attack
table and `acc_ceiling` figures stand as computed.

⚠⚠ **AND `u_ceiling` MUST NOT BECOME A CALIBRATION INPUT.** Dividing `eta_tau` by
`u_ceiling` would raise the exponent — since `u_ceiling` < 1 — and because `u_ceiling` is
a direct readout of the private input, that is local-sensitivity calibration, which
Nissim–Raskhodnikova–Smith show is not private at any `epsilon`. Ruled inadmissible
2026-10-12. The safe form, dividing by the supremum of `u_ceiling` over all admissible
inputs, gains precisely nothing: that supremum is exactly **1**, attained by any path
called at every chain position.

## 6. Signature

- [x] The two superseded figures are named with their old denominators, not silently replaced.
- [x] The cause is demonstrated, by reproducing both old values from the live arrays on the old span.
- [x] No pre-registered field is altered; R1's disposition is unchanged.
- [x] It is explicit that no privacy claim moves, and that `u_ceiling` may not enter the calibration.
- [x] The file is not named `prereg_*.md` and cannot shadow the pre-registration.

**Dylan Smith, 2026-10-12.**
