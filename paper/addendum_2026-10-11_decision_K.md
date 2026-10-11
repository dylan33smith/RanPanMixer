# Addendum 1 to `paper/prereg_2026-10-10.md` — decision K

**Status: SIGNED 2026-10-11.** Authorised by Dylan Smith on 2026-10-11 ("Decision K is
okay. Just make a dated addendum.").
`git_rev` 2fc39ebea1c6baaf8f4e3f5d568fd139f641addf — the code state this addendum attests to.

⚠ **WHY THIS IS A SEPARATE FILE AND THE PRE-REGISTRATION IS UNTOUCHED.** Standing
constraint 6 freezes that file: "If any of them must change, it is a NEW dated
pre-registration file and this one stays in the repository — superseded, never edited."
Editing a signed document after the fact is precisely what a signature is supposed to
make impossible, so nothing in `prereg_2026-10-10.md` has been altered. This file
records what changed after it was signed, and the main pre-registration remains in
force for everything not named here.

⚠ **IT IS DELIBERATELY NOT NAMED `prereg_*.md`.** `src/attacks.py:require_signed_prereg`
globs that pattern and returns the first signed, clean match, so a second file matching
it could become the document an experiment stamps itself with if the original were ever
renamed or unsigned. An addendum must never be able to impersonate the pre-registration.

---

## 1. What changed

**Decision K (Dylan, 2026-10-10, AFTER the signature).** A maximal run of chain
positions with zero genetic distance between neighbours is a block inside which the
copying model cannot switch donors, because the transition is exactly the identity
there. If no single cohort donor has a called allele at every position in the block, no
path traverses it, and the whole block now leaves the chain.

## 2. Why it was forced, not chosen

The chain as the pre-registration described it **could not be released from.** The
forward pass raised `zero total mass` on the rebuilt chr21 at position 10,543,014. This
was not a preference between two workable designs; one of them does not run.

⚠ That matters for how this addendum should be read. The fields constraint 6 exists to
protect — `tau`, the utility FUNCTION, the endpoints, the attacker set, the attack
database, the decoy design — are untouched, and none of them was inspected against an
outcome before this change. The change was forced by a crash, which is as far from
outcome-driven analysis choice as a change can be.

## 3. What it affects in the signed file

**(a) Section 7, item 6** reads: "**3.27% more** sits in interior zero-cM plateaus and
is emitted as verbatim donor runs." ⚠ **Now only partly true.** Most plateaus are still
emitted as verbatim donor runs; **20 blocks of them, 172,142 chain positions, are
EXCLUDED instead.** Of 52,870 zero-distance runs genome-wide only those 20 are
untraversable; 14 of the 20 overlap a centromere.

**(b) `beta_v`'s normalising set.** Section 2 pins Arm 2's `beta_scope` to the value
"chain". The NAME is unchanged and the rule is unchanged — normalise over the chain —
but **the chain's CONTENTS changed**, so the same pre-registered field now covers a
different set of positions:

| | before decision K | after |
|---|---|---|
| genome-wide `T` | 22,306,704 | **22,134,562** |
| mean `beta_v` per position | 4.4830e-08 | **4.5178e-08** (+0.778%) |
| chr9's `beta_share` | 0.046553 | 0.043559 (**−6.43%**) |
| chr20's `beta_share` | 0.022990 | 0.022340 (−2.83%) |

Every chromosome's share moved, because decision J allocates `beta_share` by variant
count and 19 of the 22 chromosomes lost positions.

## 4. What it does NOT affect

Unchanged and still in force exactly as signed: the `tau` grid (0, 0.05, 0.10, 0.25,
0.50, 0.90, 0.99) and each point's declared role; `w_v` = 1/support(v); the primary
endpoint and the marginals route it must use, with the derivative route still forbidden
at genome scope; the attacker set A1–A4 and their mandatory positive controls; the arm
guard; the attack database and its V1/V2 membership variants; the three decoy
conditions and seed 20261010; **`delta` = 0.056 and the balanced 25-draw design**; and
every other item in section 7.

⚠ No attack outcome had been inspected when decision K was made, and none has been
inspected as of this addendum's date.

## 5. Reporting obligation this creates

Any table quoting a chain-position denominator, a `beta_share`, or a `mean beta_v`
**must state that it is post-decision-K**, because a pre-decision-K figure differs by
0.772% and the two are not comparable. `docs/plan.md`, `docs/data.md` and
`docs/memory.md` carry in-place corrections at every figure that moved, including a
marker on the memory table's total row.

⚠ **And a consequence for ARM 1, recorded before it is reached.** Leave-one-out removes
both of the target's haplotypes from the donor panel, which lowers support(v) and can
EMPTY a block's donor intersection that was non-empty for Arm 2. **Arm 1's chain is
therefore per-target** — 44 chains, 44 `T` values, 44 `beta_share` allocations, and the
untraversable-block set recomputed per target. Reusing Arm 2's mask there would
reintroduce the dead end decision K removed.

---

**Signed:** ✅ **Dylan Smith, 2026-10-11.** `git_rev` 2fc39ebea1c6baaf8f4e3f5d568fd139f641addf.
