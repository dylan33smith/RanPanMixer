# memory.md — the linear ledger

**DO NOT read this file whole.** `grep` it.

```
grep -n "^## 20" docs/memory.md          # list entries
grep -n "INCORRECT\|CORRECTION" docs/memory.md   # what we got wrong
```

`docs/plan.md` exists so a new session does not have to come here.

## Rules

1. **Append only.** Never delete, never overwrite. New entries go directly
   BELOW the APPEND marker, so everything after the marker reads
   **newest-first**. (The rule previously read "newest at the bottom", which
   never matched the file; the file is the source of truth for history.)
2. **In-place correction.** Prepend `[INCORRECT] - ` to the wrong line, preserving
   its text verbatim; insert `[CORRECTION - YYYY-MM-DD]: ` directly below it. The
   wrong version stays — it is what makes the reasoning legible later.
3. **Provenance or it didn't happen.**
4. Metric names come from `docs/terms.md`.

---

## UNEDITED ARCHIVE — prior history

There is none. The repository's first commit (`30b3a1a`, 2026-08-29) contained two
PDFs and nothing else: the proposal and the documentation framework. No prior
notes, logs or results exist to import. This banner records the absence
deliberately, so that a later reader does not go looking for a history that was
never written.

---

## 2026-08-29 — A-FIX-docs-system: adopted the six-file documentation system

**Goal.** Stand up the documentation structure of
`archive_docs/The_Six_File_Lab_Record.pdf` on a greenfield repository, with the
content derived from `paper/Private_Genome_Path_Release.pdf`.

**Method.** Read both source PDFs. Wrote `docs/terms.md` and `docs/data.md` first
(the framework's advice: they carry the anti-drift load and cannot be generated
from existing material). Moved the two PDFs out of `docs/` so that directory holds
exactly the six files — the proposal to `paper/`, the framework to `archive_docs/`,
both with `git mv` so history follows. Wrote the verifier before any code exists.

**Result.** Six files plus `tests/test_docs_contract.py`. `CLAUDE.md` is 137 lines
against a 150-line budget and carries no findings. Every path in `docs/data.md`
except the documentation set is `PLANNED`, which is the honest state: this project
has no code, no graph, no cohort and no target genome.

**Deviations from the framework, and why.**
- The framework's Stage A / Stage B measurement split (selection vs
  characterisation) does not apply to a release mechanism. It is replaced in
  `docs/terms.md` by the denominator trap this project actually has: utility over
  ALL blocks versus over SWITCHED blocks only. Same failure mode, different axis.
- The reference verifier asserts that `memory.md` contains at least one
  `[INCORRECT]` line. On a new project that fails for the wrong reason, so that
  check skips until the first retraction exists and fails loudly thereafter.
- Two checks in the reference implementation had no analogue here (a pinned
  regression baseline, and a contested-values registry). They are replaced by a
  check that recomputes the `tau` calibration table from its own definitions —
  the closest equivalent for a project whose primary artifact is a proof.

**Source:** `paper/Private_Genome_Path_Release.pdf`,
`archive_docs/The_Six_File_Lab_Record.pdf`.

### Decision — `viterbi_release` retired before use

MAP/Viterbi decoding of the tilted HMM is not a draw from `release_distribution`
and carries no privacy guarantee, while producing output indistinguishable in
appearance from a valid release. It is the natural thing to reach for in an HMM
codebase. Retired in `docs/terms.md` so that the name resolves to its refutation
rather than to nothing. **This is a decision from the proposal's text, not a
measurement.**

### Finding — two typesetting defects in the proposal PDF

Verified by rendering page 3 at 300 dpi, not by text extraction (which mangles
superscripts and would not be evidence).

- **Theorem 1**, the pointwise likelihood-ratio bound, is typeset as
  `e^{-eps} <= Q_p(y)/Q_q(y) ^{eps}` — the upper half of the inequality has lost
  its `<= e`, leaving `eps` floating as a superscript beside the fraction. The
  intended statement is `e^{-eps} <= Q_p(y)/Q_q(y) <= e^{eps}`.
- **The proof's first line** has the same defect: `1 <= exp(eta*u(p,y))^{eta}`
  where `1 <= exp(eta*u(p,y)) <= e^{eta}` is intended.

Both are the same LaTeX slip losing `\le e` before a superscript. **The
mathematics is unaffected:** the chain still yields `TV <= tanh(eta_tau) = tau`,
and `epsilon_tau = 2*arctanh(tau)` is correct as printed elsewhere.

⚠ **Why this is recorded rather than ignored.** A one-sided bound transcribed
literally into an assertion would pass on every input, because the missing half is
the half that can fail. `docs/terms.md` states both halves explicitly and the
`log_z_p` bracket `[0, eta_tau]` is written as a two-sided assertion for the same
reason. Fix the manuscript when `paper/` is next edited.

**What this does NOT establish:** nothing about the correctness of the mechanism.
It is a defect in the document, found by reading it, before any code was written.

---

<!-- APPEND NEW ENTRIES BELOW THIS LINE -->

## 2026-10-10 (later) — B-THY-diploid-gate: the falsifier now covers the path that runs

**Hypothesis.** Every falsification clause in this project was HAPLOID and
single-chromosome, while the diploid path was the only one ever run and decision E
makes Arm 1 diploid across 22 autosomes under one genome-wide budget. If that gap
mattered, injecting a real calibration error into the diploid code path would pass
the existing suite. Chunk 1 of the review-gated plan.

**Method.** Added the diploid and multi-segment reference primitives `exact.py` did
not have (`Factor`, `composite_release`, `composite_utility`,
`assert_log_z_bracket`) plus a two-segment toy fixture (`toy.toy_segments`,
`toy.segment_shares`). Then `tests/test_diploid_gate.py`, 64 clauses in five groups:
diploid on the 96-output toy (product support 9,216); two segments haploid,
**EXHAUSTIVE** over all 32 composite inputs and all 1,024 ordered pairs; diploid x
two segments (four chains, the Arm 1 shape); three mutations; and the clauses tying
the reference to `ffbs`. Then mutation-tested by editing the real modules, running
the gate, and reverting — each revert verified byte-identical with `cmp`.

**The algebra the gate encodes.** A release is a product of independent chains, one
per (haplotype, chromosome). Each factor is tilted at `eta_tau / ploidy` with a
`beta` summing to that chromosome's SHARE `w_c`, so `u` stays in [0,1] and `log Z`
in [0, `eta_tau`] however many chains there are. **One invariant catches both
composition errors and is therefore asserted rather than documented:
`sum over all factors of beta.sum() == ploidy`.**

### Which clauses fire under which injected error

Clauses are rows; the injected error is the column. ✓ = fired, — = did not.

| clause | unmutated | M1 full `eta` per haplotype | M2 `beta` per chromosome | M3 `epsilon_tau` for `eta_tau` |
|---|---|---|---|---|
| budget invariant (`sum beta == ploidy`) * | — | ✓ | **✓** | — |
| composite `u` in [0,1] * | — | ✓ | **✓** | — |
| `log_z` bracket on the TOTAL * | — | ✓ at `tau` = 0.9 only | ✓ | — |
| per-segment `log_z_c <= eta_tau * w_c` | — | n/a | ✓ | — |
| exact TV <= `tau` * | — | — | — | — |
| `ratio_spread <= exp(epsilon_tau)` * | — | ✓ | — | **✓** |
| the sharp pair SATURATES the bound * | — | ✓ | — | ✓ |
| **`ffbs` diploid `log_z_p` == exact composite** * | — | **✓ 5 of 6 `tau`** | n/a | **✓ 5 of 6 `tau`** |
| **`ffbs` diploid marginals == exact single chain** * | — | **✓** | n/a | **✓** |

**Provenance.** `tests/test_diploid_gate.py` · toy K=4 T=7 (96 output paths) and
segments K=3 T=3 / K=3 T=2 with declared shares (0.7, 0.3) · `tau` grid
0 / 0.05 / 0.10 / 0.25 / 0.50 / 0.90 · constructed extremal pairs plus a seeded
sample, seed 20261010, 400 pairs · M1 and M3 injected into `src/mechanism/ffbs.py`,
M2 into `src/mechanism/exact.py` · 64 clauses, 166 tests in the suite, all passing.

**COLUMNS** · *unmutated* fires nothing, which is the only column where that is the
right answer. · *M1* is the SUM aggregator in disguise — the defect that was accepted
at every `tau` before a diploid bracket existed. · *M2* is the genome-wide error:
normalising `beta` to 1 per chromosome, which turns one release into `n_chrom`
releases of one person. · *M3* is the one-character `epsilon_tau`/`eta_tau` slip.

**ROWS** · *budget invariant* and *composite u* are the only clauses that see M2 at
all, and both were verified to have power by disabling each in turn — disabling
either breaks the M2 test, disabling both also breaks M1 at four of five `tau`. ·
*exact TV fires on nothing*, which reproduces the project's existing finding at
diploid scale: TV reaches only a fraction of `tau`, so it cannot detect a mis-scaled
exponent and must never be the headline clause. · *the two `ffbs` rows are the ones
that caught M1 and M3*, and that is the finding below.

**SYNTHESIS — AND THE THING I ALMOST SHIPPED.** The first version of this gate
checked `exact.py` against itself and nothing else. It passed all 57 clauses, and
**it would have caught neither M1 nor M3**, because both live in `ffbs.py` and
nothing in the gate read `ffbs.py`. A reference implementation agreeing with its own
gate is two unverified things agreeing. The clauses that give this gate its power are
the two that tie the reference to the mechanism, and the stronger of them needs no
draws at all: **`log_z_p` is deterministic**, so comparing it between `ffbs.release`
and `exact.composite_release` pins the exponent, the per-haplotype factorisation and
the baseline subtraction to 1e-12 with zero Monte-Carlo error. The sampled-marginal
clause is the complement — it checks that the DRAW follows the factorisation, not
merely that the normaliser does, and it is what a shared RNG stream between
haplotypes would fail.

⚠ **Three honest limits on what this gate covers.**
1. **At `tau` = 0 every mutation is undetectable**, because `eta_tau` = 0 makes the
   halved and the full exponent identical. That is correct, not a hole, and the
   mutation clauses exclude `tau` = 0 rather than claiming coverage there.
2. **The `log_z` bracket alone catches M1 only at `tau` = 0.9.** Below that the
   doubled exponent stays inside [0, `eta_tau`] and the bracket is silent — the same
   shape as the recorded haploid result that a doubled exponent passes clause (b)
   everywhere and clause (c) below `tau` = 0.5. The layered detectors are load-bearing.
3. **M2 cannot yet be injected into the mechanism**, because the genome-wide share
   does not exist in `utility.beta_v` — its assertion still hard-demands `sum == 1`.
   So M2 is caught at the REFERENCE level only, and the gate is forward-looking: it
   will bite when chunk 2 implements the share. That ordering is deliberate and
   matches how the haploid toy gate was built before the sampler it checks.

**One real fix found on the way.** `toy_segments`'s own fixture check rejected my
first segment B: the invariant "some donor is another donor minus a cell" was written
as "some donor is DONOR 0 minus a cell", and B3 is B2 minus a cell. The narrow version
would pass a fixture whose only collapsing pair sat elsewhere, so the collapse that
makes TV a statement about OUTPUT paths rather than state paths would go untested.
Generalised to any ordered pair. ⚠ **`toy_cohort` has the same narrow check** and is
left alone for now — it happens to be satisfied by donor 0 on that fixture, so
widening it is a no-op today, but it is the same latent weakness.

## 2026-10-10 — CHUNK 0: the known-wrong facts landed, and one of them was a real measurement

**Goal.** Clear the stale and wrong facts out of the authority documents before any Phase B
code reads them. Chunk 0 of the review-gated plan; no new capability, only corrections.

**What was corrected, and the two that mattered.**

**1. The float32 headroom pair, wrong in five live places.** `1.61e-06` and `13.5x` were
computed over the **340,849-record site axis** instead of the chain's T = 305,886. The
denominator is the chain — `eta_tau * mean(beta_v) = eta_tau / T` — so at `tau` = 0.5 it is
0.549306 / 305,886 = **1.80e-06**, i.e. **15.1x** float32's 1.19e-07, not 13.5x. The old
pair understated the tilt's headroom by 10%. Fixed at `docs/plan.md`, `docs/terms.md`,
`tests/test_mechanism.py` and `scripts/results_table.py`; corrected in place in this file.
⚠ **The conclusion is stronger than the old text claimed, not weaker:** under decision B
`beta_v` sums to 1 genome-wide, which drops the per-position tilt BELOW float32 eps
entirely, so float64 is not a 15x margin — it is the only thing preventing a silently
untilted draw. In `docs/bugs.md`.

**2. The `21 TB` dense-transition figure, wrong by 1000x in three places.** `340,824 x 88 x
88 x 8 B` is **21.1 GB**. It stood in `src/mechanism/ffbs.py`, `src/mechanism/hmm.py` and
this file, with `hmm.py` printing its own arithmetic one line below the wrong unit. It
survived because the conclusion it supports — never materialise `A` — is correct either
way. **A number supporting a correct conclusion is never checked**, and that is the general
lesson, not the unit slip.

**3. `15.03%` -> `14.8077%`.** Same statistic, same denominator (cells differing from the
target, over 507,692 target-known cells, HG00096 diploid chr21), but 15.03% is 76,283 from
ONE draw at seed 1000 and 14.8077% is 75,177.401 computed exactly from the posterior
marginals with zero Monte-Carlo error. The sampled draw ran **1,106 cells high**, ~1.5%
relative. Corrected on the live board; a `[CORRECTION]` block at the head of the 2026-10-03
entry covers the seven copies there, whose conclusions all survive the substitution. ⚠ The
`capacity ~= 0.17` crossing needed **no** change — the capacity sweep had already used the
exact 14.81%, so only one plan.md line was ever out of step.

**4. `K_states` and `donor_switches` stamped by arm.** `K = 88` is Arm 2; `K = 86` is Arm 1
under per-target leave-one-out. Also removed a stale "whether duplicates OCCUR on this
cohort is unmeasured — measure it before quoting K again": it WAS measured 2026-09-30, all
88 haplotypes are distinct over the chain and the closest pair differs at 32,741 positions
(10.7%). The gap is a missing assertion, not a wrong number.

### ⚠ AND ONE "UNMEASURED" THAT WAS A REAL MEASUREMENT WAITING TO BE TAKEN

Yesterday's pass wrote, in two places, that the split of chr21's 34,371 dropped-below
positions into map-silent and in-span-plateau was UNMEASURED. An independent verifier
flagged it as a false UNMEASURED; chasing it down found **two totals on record that did not
reconcile** (34,371 and 34,345). Rather than pick one, measured directly.

**Method.** `data/cohort/chr21/positions.npy` (340,849 records) and
`data/cohort/chr21/support.npy` against the extracted `external/genetic_maps/chr21.b38.gmap`
(span 10,326,676 - 46,680,243). Counted records by class and checked `support > 0` in each.

| class | positions | `support > 0` | declared rate reaches it? |
|---|---|---|---|
| below the map start (< 10,326,676) — map-silent | **26,044** | 26,043 | **yes** |
| in [10,326,676, 12,968,320] — in-span zero-cM plateau | **8,327** | 8,327 | **no** |
| sum | **34,371** | — | — |
| above the map end (> 46,680,243) — map-silent | **592** | — | **yes** |

**Provenance.** chr21/hprc-v1.0-pggb, our 340,849-record build · 305,886 chain positions ·
measured 2026-10-10 · the plateau interval contains exactly **one** distinct cM value among
all map markers (0.584144), confirming it is a genuine map plateau and not an artifact.

**Result, and it resolves the discrepancy.** The sum matches the recorded 34,371 exactly.
**26,018 is PanMixer's vg-1.36 340,824-record figure**, and pairing it with 8,327 gives
34,345 — which is why a second total existed. Ours is 26,044.

⚠⚠ **THE CONSEQUENCE IS A LIVE RESIDUAL WITH A SIZE: 8,327 positions, 2.7% of the chain.**
Decision C re-admits everything the cut dropped, but its declared rate reaches only the
map-silent class. The 8,327 in-span-plateau positions come back with `dcM = 0`, hence
`P(switch) = 0`, hence a **verbatim donor run** — and all 8,327 are usable (`support > 0`),
so none is filtered out downstream. **So decision C restores part of the p arm under the
declared rate and part of it under zero switching.** Whether that is acceptable is open
question 2 and is Dylan's; what changed today is that it is a known quantity rather than an
unknown. In Arm 1 the donors ARE the protected population, so a verbatim run there is a
privacy cost, not merely a utility one.

**The same split on the other 21 autosomes remains unmeasured** and is `B-DAT-map-coverage`
(chunk 3). Four more acrocentrics — chr13, chr14, chr15, chr22 — will have the same shape.

**Gate.** 102 tests pass. Every distinctive stale number greps clean in live text, and each
historical copy carries an `[INCORRECT]`/`[CORRECTION]` pair rather than an overwrite.

## 2026-10-09 — PHASE B OPENED: the two-arm paper, and the four measurements that sized it

**Goal.** Dylan set a new scope on 2026-10-08: the paper becomes two sections. The
first is a direct head-to-head against PanMixer, obfuscating donors from within the
same 44-person pangenome using our mechanism instead of theirs, with their own
evaluations so the numbers are comparable. The second is the external-target work of
Phase A carried to genome scale. This entry records the decisions, the measurements
taken to size them, and two defects the audit found in our own plan.

**The decisions (A-I) and the seven items still open live in `docs/plan.md`** under
"THE DECISIONS OF 2026-10-08" — not duplicated here. They were taken after a
seven-surface, 71-agent audit of what the change implies; 20 residual decisions
survived adversarial refutation and collapsed to five distinct forks plus five
classes nobody had asked, of which Dylan answered three.

### M1 — autosomal genetic-map coverage: the chromosome-ends problem is a chr21 artifact

**Method.** Extracted all 25 files from `external/genetic_maps/genetic_maps.b38.tar.gz`
(SHAPEIT4 b38; the archive holds all 22 autosomes plus chrX/PAR, so the 21 missing
maps are an EXTRACTION, not a download). For each autosome, parsed the three-column
marker table (`pos`, `chr`, `cM`), compared first and last marker positions against GRCh38 primary assembly
chromosome lengths, and summed intervals where consecutive markers carry identical cM.
Script in session scratchpad, not committed.

| quantity | autosomal total | share |
|---|---|---|
| total length | 2,875.0 Mb | — |
| bp before the first marker | 84.2 Mb | 2.93% |
| bp after the last marker | 2.0 Mb | 0.07% |
| **no cM at all (lead + trail)** | **86.2 Mb** | **3.00%** |
| interior zero-cM intervals (inside the span) | ~93.9 Mb | 3.27% |

- **98.5% of all lead bp is the five acrocentric p-arms**: chr13 18.45, chr14 19.03,
  chr15 19.81, chr21 10.33, chr22 15.29 Mb = 82.91 of 84.2 Mb.
- Every non-acrocentric autosome has a lead under 0.25 Mb and a trail under 0.25 Mb.
- Largest interior zero-cM extents: **chr9 25.61 Mb**, chr16 10.69 Mb.
- ⚠ **chr21 is the WORST autosome by fraction** — 10.33 of 46.71 Mb is 22.1% lead,
  27.8% to the end of its flat plateau. The chromosome we built everything on
  overstates this problem by roughly 8x.

**Reading.** The exclusion rule costs 3.00% genome-wide, not the ~10% chr21 implied,
and it is a UTILITY decision rather than a privacy one: `R_D` stays target-independent
under any map choice. This is what made decision C affordable. ⚠ The interior
zero-cM total is LARGER than the ends problem, sits INSIDE `chain_span`, and is where
the chain cannot recombine — so a donor is emitted verbatim there. Decision C
deliberately does NOT touch it; see open question 2.
⚠ **Position counts are UNMEASURED.** These are base pairs. How many chain positions
sit in each class, and how many survive `support(v) > 0`, is what `B-DAT-map-coverage`
must produce before the declared rate can be given a value.

### M2 — Arm 1's representability ceiling: a cohort member is reconstructible from the rest

**Method.** Loaded `data/cohort/chr21/haplotypes.npy` (88 x 340,849 int16, -1 missing),
restricted to `in_chain.npy` (T = 305,886). Verified from `haplotype_ids.npy` that
haplotypes 2i and 2i+1 belong to individual i. For each of the 88 haplotypes, removed
BOTH of that individual's haplotypes (K = 86, decision E) and asked whether any
remaining haplotype carries the same allele at the same position.

| metric | value |
|---|---|
| target's called alleles carried by >= 1 of the other 86 | **99.51%** (min 99.12, max 99.78) |
| target's non-ref alleles carried by >= 1 of the other 86 | **97.33%** (min 95.71, max 98.57, median 97.23) |
| mean called alleles per haplotype | 304,073 of 305,886 chain positions |
| mean non-ref alleles per haplotype | 54,712 |

**Reading.** Arm 1 is VIABLE, and its ceiling sits ABOVE Arm 2's (~90% of non-ref
calls have a graph record at all). ⚠ The two are therefore NOT comparable to each
other — see the `u_ceiling` and `target_fidelity` entries in `docs/terms.md`.

### M3 — and the mechanism protects the rarest alleles by construction

**Method.** Same arrays and leave-one-out as M2, stratified by cohort-wide allele
count over all 88 haplotypes and weighted by `-log f` with f = AC/88.

| stratum (AC over 88 haplotypes) | non-ref alleles | unrepresentable under LOO | share |
|---|---|---|---|
| AC = 1 (private singleton) | 126,065 | 126,065 | **100.0%** |
| AC = 2-4 | 231,577 | 4,068 | 1.8% |
| AC = 5-8 | 257,031 | 0 | 0.0% |
| AC > 8 | 4,199,977 | 0 | 0.0% |

Per haplotype: **1,479 of 54,712** non-ref alleles (2.70%) are carried by no other
individual, and they hold **11.03%** of that haplotype's total `-log f` information —
a **4.08x** over-representation.

**Reading, and the predicted Arm 1 headline.** Under leave-one-out the sampler
**cannot emit a cohort singleton at any `tau`**, because no donor state carries it.
That is protection by CONSTRUCTION, not by `tau`. PanMixer starts from a verbatim copy
of the target's own column (`stacker.py:53`) and releases whatever the knapsack did not
select, and its privacy score reads only anchor alleles, so the most identifying sites
are the ones most likely to survive. **We resample everything; they resample a selected
subset.**
⚠ **UNVERIFIED — the comparison half.** PanMixer's survival rate for this same allele
class has NOT been measured. The claim above is a prediction, not a result.

### M4 — float64 at genome scale: the hazard is the accumulator, not the potentials

**Method.** Arithmetic, not a run. `mean(beta_v) = 1/T`; chr21's T = 305,886 gives
3.269e-06. A genome-wide T is UNMEASURED; at 50-62x chr21 the mean falls to roughly
5e-08 to 6.5e-08. The float64 boundary is `exp(x) == 1.0` iff `x <= eps/2` =
1.110e-16.

- At `tau` = 0.5 the genome-wide exponent sits ~1.6e+08 ABOVE that boundary. **No
  precision problem**: T would need ~5e+15 positions to reach it.
- ⚠ **But float32 eps is 1.19e-07, and the genome-wide exponent falls BELOW it.** A
  dtype downgrade would make every tilt potential exactly 1.0 and silently turn the
  mechanism into an untilted `R_D` draw. At chr21 the same ratio is 15x — merely
  marginal. The dtype refusal in `ffbs.forward` stops being a preference and becomes
  the only thing standing between us and a mechanism that does nothing.
- ⚠ **The real cancellation hazard is the `log_c` ACCUMULATOR.** `log_z_p` is a
  difference of two sums whose magnitude grows with T. Materialised as a `(T,)` array
  and summed by `np.sum`'s pairwise algorithm, the error stays ~1.5e-10 at genome
  scale. Accumulated as a running Python float — the obvious fix when the alpha memory
  bites — it reaches **~9.9e-05 against a 1e-6 tolerance**, a 100x breach that passes
  every chr21 test. `math.fsum` or Kahan costs nothing and should go in before the
  driver.

### What PanMixer actually evaluated (gates how Arm 1 must be built)

From `paper/s41467-026-77591-0_reference.pdf` and the pinned `external/PanMixer`.

- **DIPLOID, not haploid:** *"Each individual genome inside the pangenome graph has
  two haplotype paths, which we obfuscate at the same time"*, with all downstream
  analysis at the genome level.
- **44 targets, ALL cohort members:** *"applying PanMixer 44 times for all 44
  individuals in the HPRC v1 pangenome reference graph created using PGGB"*.
- **22 autosomes.** chrX stripped by `remove_X.sh`; chrY never mentioned in either PDF.
- **`capacity` is applied PER CHROMOSOME** — against our genome-wide `tau` (decision
  B), which makes us strictly more conservative in scope. Say so rather than matching.
- The release object **starts as a verbatim copy of the target's own column**
  (`stacker.py:53`); only knapsack-selected blocks are overwritten.
- **50 log-spaced capacities**, 1.0 -> 0.000110.
- ⚠ **R = 1, and the SAME presampled draw is reused across all 50 capacities** —
  `stacker.py:57` loads one `new_haplotypes.npy` from the shared starting-data path
  inside every per-capacity run. Their published frontier is nested subsets of one
  sample, not 50 independent releases. This confirms the suspicion recorded in
  `docs/plan.md` and is the specific file:line that proves it.
- ⚠ The published seed is **claimed to be in the repository and is not there.**
- **Nothing in their suite tests two-target indistinguishability** — a tree-wide grep
  for `decoy|two_target|indistinguish|total variation` returns zero hits. The
  structural reason is in their own code: two inputs differing in an UNSELECTED block
  have disjoint output supports.
- Genome-wide figures: the attack and AF/LD. **chr21-only:** read mapping and the
  Beagle attack. And `starting_data/` contains **chr21 only**, so the genome-wide gap
  score cannot run as shipped — decision F's rebuild is a precondition, not a polish.

### Two defects the audit found in OUR plan, both corrected in place

1. **`docs/plan.md`'s Reporting contract carried a standing-constraint violation.** It
   read *"The head-to-head aligns the arms on MEASURED empirical attack success
   only."* Selecting the comparison point after inspecting attack outcomes is exactly
   what standing constraint 6 forbids — the plan held the violation it warns about
   everywhere else. **Corrected:** the arms align on `target_fidelity`, which is
   computable with no attack run; empirical attack success is the OUTCOME, never the
   alignment. Both full frontiers are still reported.
2. **`A-EVL-headtohead`'s exit gate demanded two things model C cannot have.** It
   required per-subject recomputation of *allele frequencies* and *block assignment*.
   Model C has **neither** — no allele-frequency table, no blocks — so those clauses
   could never have been met. Removed, with the reason recorded in the item.

### The standing constraints were renumbered, and the stale references hunted

Inserting the genome-wide `beta_v` rule as constraint 4 shifted the old constraints
4-8 to 5-9. **15 stale references** across `docs/plan.md`, `docs/terms.md`,
`src/mechanism/` and `scripts/` were corrected by line-targeted replacement with
per-line verification of the expected text. `docs/memory.md` was deliberately NOT
changed: its entries were accurate when written. The renumbering itself is recorded in
the `tau` entry of `docs/terms.md`.

### ⚠ OUTSTANDING — not done, do not assume otherwise

- ✅ **DONE 2026-10-10 (chunk 0).** The 15.03% differ-rate. `docs/plan.md` quoted it as
  ours; the exact artifact gives **14.8077%**. Corrected in place at `docs/plan.md`, and a
  `[CORRECTION]` block added at the head of the 2026-10-03 entry covering all seven copies
  there. The matched-capacity crossing (~0.17) needed NO change — it was already computed
  against the exact 14.81%, so only the one plan.md line was out of step.
- `docs/bugs.md` has no entry for the float32/accumulator hazards of M4.
- No code has been written for Phase B. The work order is in `docs/plan.md`.

## 2026-10-08 (later still) — the directly comparable number, in both units

**Dylan's question.** Our table shows ~74,000 changes. What is that comparable to in
PanMixer's table?

**Answer: PanMixer's `cells changed` row, and nothing else in it.** `moves selected` is
not comparable to anything of ours — we have no selection step — and the two must never
be put side by side. The comparable quantity is `target_fidelity`: cells where the
released allele differs from the target's own, which is computable on either output.

| | per CELL (variant x strand) | share | per VARIANT (>=1 strand) | share |
|---|---|---|---|---|
| **ours, model C**, `tau` = 0.5 | **74,988** of 506,788 | **14.80%** | **61,750** of 253,394 | **24.37%** |
| **ours, model B**, `tau` = 0.5 | 74,311 of 506,788 | 14.66% | 60,138 of 253,394 | 23.73% |
| PanMixer, `capacity` 0.01 | 7,838 of 650,720 | 1.20% | 5,360 of 336,677 | 1.59% |
| PanMixer, `capacity` 0.1 | 68,514 | 10.53% | 56,033 | 16.64% |
| PanMixer, `capacity` 0.2 | 107,207 | 16.48% | 86,069 | 25.56% |
| PanMixer, `capacity` 0.5+ (saturated) | 108,027 | 16.60% | 86,578 | 25.72% |

**Provenance.** OURS: HG00096 (external), chr21, `tau` = 0.5, seed 1000, diploid, common
axis of 304,664 positions. PANMIXER: HG00438 (a cohort member), chr21, seed 123, my
recount from `new_haplotypes.npy` against row 0 of `pangenome.npy`.

**SYNTHESIS.** **Our ~75,000 lands inside PanMixer's range of 7,838 to 108,027**, and the
matched-fidelity capacity is **~0.17 on the per-cell axis and ~0.19 per-variant** — close
enough that either is a fair statement of where the arms meet. That crossing point is the
whole value of this axis: it is the only place the two tools can be compared at all,
because `capacity` and `tau` are not comparable parameters.

⚠ **THE SAME MEASUREMENT, OPPOSITE MEANING, AND THIS IS THE POINT TO MAKE OUT LOUD.**
For PanMixer a changed cell is DELIBERATE — the release starts as the target's own
genotype column and a subset is rewritten, so changes are the product and more of them is
more obfuscation. For us a differing cell is a RESIDUAL — the release never contained the
target's data, it is a cohort mosaic, and a difference means the mosaic failed to
reproduce the target there. **14.80% of ours is reached at `tau` = 0**, so almost none of
it is caused by the mechanism. Same number, opposite sign of intent.

⚠ **THE DENOMINATORS ARE NOT EQUAL AND THE TARGETS ARE DIFFERENT PEOPLE.** Three reasons,
all measured: (1) our chain excludes **34,937 of their records (10.25%)**, the p arm and
the telomere; (2) we lose **104,080 cells (17.0%)** where our join finds no record for the
target, against about 4.5% for them, because their target is a cohort member read from
the same VCF; (3) HG00096 is external and HG00438 is in the cohort. **So compare the
SHARES, treat the raw counts as indicative, and do not quote a difference of counts as a
result.**


## 2026-10-08 (later) — selected moves vs actual changes, swept; and a clash with the recorded 30.9%

**Dylan's question.** `moves selected` counts items taken from the knapsack, but a
resampled block can come back carrying the target's own alleles, so many are effectively
no moves. How many moves actually CHANGE the target's haplotype?

**Method.** A move is a (block, strand) slot of `xsol.npy`. It is a no-op if no site in
that block changed on that strand, i.e. `(new_haplotypes != original) & (original != -1)`
is false throughout. Blocks from `simple_blocks_idx.npy`; original is row 0 of
`pangenome.npy` (HG00438). Comparison arm only; nothing in `src/` reads any of it.

| `capacity` | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| moves selected | 14,342 | 34,340 | 47,070 | 85,264 | 159,225 | 177,865 | 177,865 | 177,865 |
| **no-op moves** | 13,820 | 26,673 | 31,315 | 61,532 | 126,177 | **142,718** | 142,718 | 142,718 |
| **no-op share** | **96.4%** | 77.7% | **66.5%** | 72.2% | 79.2% | **80.2%** | 80.2% | 80.2% |
| moves that changed something | 522 | 7,667 | 15,755 | 23,732 | 33,048 | **35,147** | 35,147 | 35,147 |
| cells changed | 7,838 | 33,698 | 68,514 | 107,207 | 108,012 | **108,027** | 108,027 | 108,027 |

**Marginal cost of raising the knob:**

| step | extra moves | of which no-ops | extra cells changed |
|---|---|---|---|
| 0.01 -> 0.05 | +19,998 | 64% | +25,860 |
| 0.05 -> 0.1 | +12,730 | 36% | +34,816 |
| 0.1 -> 0.2 | +38,194 | 79% | +38,693 |
| **0.2 -> 0.3** | **+73,961** | **87%** | **+805** |
| **0.3 -> 0.5** | **+18,640** | **89%** | **+15** |
| 0.5 -> 0.9 | +0 | — | +0 |

**Provenance.** `xsol.npy` and `new_haplotypes.npy` under `/data/ds85/RanPanMixer/runs/capacity_sweep/`, one directory per capacity
· `external/PanMixer` @ `c182c38`, unmodified · HG00438, seed 123 · candidate list
100,757 blocks x 2 strands = 201,514 slots.

**SYNTHESIS.** Dylan's reading is right and the magnitude is the finding: **at
saturation, 80.2% of the 177,865 selected moves change nothing at all.** Only **35,147**
moves ever touch the genome, and they account for all 108,027 changed cells. The no-op
share is WORST at the extremes -- 96.4% at `capacity` 0.01, where the LP buys the
cheapest moves and the cheapest moves are the ones that resample a common allele back to
itself -- and the marginal moves at the top are almost pure no-ops: **the 92,601 moves
added between `capacity` 0.2 and 0.5 are 87-89% no-ops and buy 820 changed cells
between them.** So `moves selected` is close to meaningless as a measure of how much the
release was altered, and the two should never be quoted as if they tracked each other.

⚠⚠ **THIS CLASHES WITH THE 30.9% RECORDED ON 2026-09-24 AND I HAVE NOT RESOLVED IT.**
That entry reports "of 47,070 selected moves, 14,526 changed NOTHING (30.9%)" at
`capacity` 0.1; I measure **31,315 (66.5%)** at the same capacity and seed. Both cannot
be right. Candidate explanations, none verified: a different notion of "changed"
(perhaps restricted to the block's ANCHOR variants rather than all its sites); a
different block mapping; or different handling of blocks where the target's genotype is
missing throughout, which my definition counts as no-ops. **Do not quote either figure
as settled until this is closed.** What is NOT in doubt, because it does not depend on
the no-op definition, is the cells-changed row and the saturation it shows.

⚠ Also unresolved, and noted before: my recount gives 108,027 changed cells at the
ceiling against the stacker's logged 106,549, a 1.39% gap. The columns above use my
recount throughout so they are internally consistent.

⚠ **AND THE RIGHT DENOMINATOR FOR MY RECOUNT IS 650,720, NOT THEIRS.** Cells with a
called target allele across all 340,824 variants x 2 = **650,720**. So my 108,027 is
**16.60%**, while the stacker's 106,549 over its own 649,242 is **16.41%**. Both are
internally consistent; crossing them (my numerator over their denominator) gives a
spurious 16.64% and must not be quoted.

### How a move relates to a cell, and a partial explanation of the no-op gap

A **move** is one (block, strand): 100,757 x 2 = 201,514 candidates. A **cell** is one
(variant, strand): 681,648 slots, 650,720 with a called allele. Blocks hold median 1,
mean 3.38, max 618 variants, so **one move can change anywhere from 0 to 618 cells.**

At saturation: 177,865 moves selected, 142,718 no-ops, **35,147 changing moves**, which
carry **104,390** of the 108,027 changed cells at median 1 / mean 2.97 / max 262 cells
each. **68.7% of changing moves altered exactly one cell** and 6.0% altered ten or more.

⚠ **The remaining 3,553 cells (3.3%) changed in slots that were NOT selected**, plus 84
at sites outside any block id. That is possible because **PanMixer's blocks are not
disjoint** — `docs/memory.md` 2026-09-22 measured 9.2% of records sitting inside another
record's span and 15.4% of blocks containing such a record — so a write from a selected
block can land on a site whose own block id is unselected, and my per-block attribution
scores that move as a no-op.
**But 3.3% cannot close a 2.2x gap**, so this is a minor contributor and the
66.5%-versus-30.9% discrepancy stays unexplained.


## 2026-10-08 — churn of the mismatch set across tau, and a model-B release bug

**Dylan's question.** The exact `E[differ]` figures are NET. Is the tilt moving a
handful of positions, or reshuffling many and netting to a handful? Both models, same
seed.

### ⚠ FIRST, A BUG: model B could not produce a release at all

`from_arrays(..., blocks=...)` DROPPED the 4 units no donor can cover, which left
`unit_spans` indexing the full 305,886-position chain while `n_positions` was 304,664.
Two coordinate systems. `release()` refused outright (the constraint-1 shape guard
caught it), so model B was unusable for actual releases — and **every model-B fidelity
share reported on 2026-10-07 used the full chain's denominator while its numerator
covered only 304,664 positions.**

**Fixed:** `from_arrays` now REFUSES dead units and names the remedy, and
`hmm.dead_block_positions()` returns the mask to compose into `in_chain` — so the
exclusion lives on the axis every artifact shares, exactly as `support == 0` does.
Re-measured on the COMMON axis both models can emit (304,664 positions, 506,788 known
cells), the floors barely move and the gap survives: **C 14.8073%, B 14.5959%**. So
model B's advantage was not a denominator artifact.

### The churn, same seed, mean of 3 seeds

| `tau` | model C mismatches | net | churn | model B mismatches | net | churn |
|---|---|---|---|---|---|---|
| 0.0 | 74,760.0 | — | 0 | 73,175.7 | — | 0 |
| 0.05 | 74,760.0 | +0.0 | **0** | 73,175.7 | +0.0 | **0.7** |
| 0.10 | 74,760.0 | +0.0 | **0** | 73,175.7 | +0.0 | 0.7 |
| 0.25 | 74,756.0 | -4.0 | 5.3 | 73,175.7 | +0.0 | 0.7 |
| 0.50 | 74,756.0 | -4.0 | 5.3 | 73,175.7 | +0.0 | 1.3 |
| 0.90 | 74,755.3 | -4.7 | 6.7 | 73,185.7 | +10.0 | 19.3 |
| 0.99 | 74,747.3 | -12.7 | **40.7** | 73,199.0 | +23.3 | **84.7** |

`churn` = positions that moved in either direction (fixed + broken); `net` = broken −
fixed, which is what `E[differ]` reports.

⚠ **EVERY CELL IS A MEAN OVER 3 SEEDS OF ONE TARGET (HG00096), NOT OVER THE 10 TARGETS**
— which is why the figures are fractional, and the table header now says so rather than
only the provenance line. **And the per-seed spread is enormous:** model C's churn at
`tau` = 0.99 is **35 / 83 / 4** across seeds 1000/2000/3000, and model B's is
**72 / 3 / 179**. A factor of 20 to 60 between seeds, so the mean of 3 is a weak summary
of a highly variable quantity — churn depends on whether a flip lands in a long donor
segment or a short one, and segment length is itself the thing that varies. Treat the
churn means as an order of magnitude, not a measurement; the exact `E[differ]` figures
are the quotable ones.

**Provenance.** `scripts/churn_across_tau.py`, `logs/A_churn_across_tau.md` · chr21 ·
HG00096 diploid · seeds 1000/2000/3000, each `tau` column drawn from the SAME uniform
stream within a seed · common axis 304,664 positions.

**SYNTHESIS — the answer is "a bit of both, and still tiny".** The distribution DOES
shift more than the net suggests: at `tau` = 0.99 model C moves **40.7** positions to
net **−12.7**, and model B moves **84.7** to net **+23.3**. So roughly **3x more
positions move than the net change**, because the tilt both fixes and breaks matches.
But the absolute scale is unchanged: tens of positions out of 304,664. And below
`tau` = 0.25 model C's churn is **exactly zero** — the released path is bit-identical,
which is the thin-tilt finding seen position by position rather than in aggregate.

⚠ **The sampled nets are noise and the exact values are authoritative.** Exact
`E[differ]` at `tau` = 0.99 gives **−3.30** cells for C and **−4.72** for B; the
sampled −12.7 and +23.3 are 3-seed estimates against a paired sd near 50. Model B's
positive sign is sampling noise, not a reversal — the exact figure is negative, i.e.
the tilt helps. Quote the exact numbers.

**AND THE CHANGES CLUSTER, which is the structural finding.** They are not spread
across the chromosome: **12 of model C's 23 changed positions fall in one 400-variant
window, and 51 of model B's 64.** That is what the mechanism predicts — a flip changes
which donor is copied, and the donor persists for a segment, so consecutive positions
change together. The tilt does not nudge positions independently; it occasionally swaps
one donor segment for another.

### Whole-chromosome mismatch map, added 2026-10-08 after Dylan asked why the figure showed only 200 variants

⚠ **The first figure showed a 200-variant window — 0.07% of the chromosome — chosen for
maximum churn, so BOTH its coverage and its local density were unrepresentative.** Its
window density ran about 20% against the chromosome-wide 14.5%. Replaced by a binned
whole-chromosome map; `scripts/mismatch_map_chr21.py` emits the data, with the window
kept only as a labelled detail because 540 px cannot resolve 304,664 variants.

| model | mean density | min bin | max bin | max shift caused by `tau` | bins holding a flip at `tau` = 0.99 |
|---|---|---|---|---|---|
| C | 0.1451 | 0.0000 | **0.4071** | **0.0077** | **4 of 540** |
| B | 0.1430 | 0.0000 | **0.4051** | **0.0293** | **4 of 540** |

**Provenance.** `logs/A_mismatch_map_chr21.md` · 540 bins of 564 variants over the
304,664-position common axis · density PER CELL over 506,788 target-known cells · seed
1000, every `tau` from the same uniform stream.

**SYNTHESIS, and it is a better picture than the window gave.** Mismatch density has
**strong regional structure that has nothing to do with `tau`**: bins range from
**0.0% to 40.7%**, so some 564-variant stretches are reproduced almost exactly and
others are nearly half wrong. That is ancestry and local haplotype sharing between the
target and the cohort, visible as banding. Against that, **the largest shift `tau`
causes anywhere on the chromosome is 0.0077** (model C) — two orders of magnitude below
the regional variation, which is why the seven `tau` rows are indistinguishable by eye.
And the flips are confined to **4 of 540 bins** in both models. So the whole-chromosome
view says the same thing as the exact derivative, in a form that also shows what DOES
vary: region, not `tau`.


## 2026-10-07 (later) — why PanMixer's knob saturates, and whether it moves the same cells

**Dylan's question.** A knapsack with a bigger budget should buy more. Why does
`capacity` saturate, and are the same alleles changed at every setting?

**A. THE KNAPSACK RUNS OUT OF ITEMS, not out of budget.** From the optimizer's own
log lines:

| `capacity` | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| moves SELECTED | 14,342 | 34,340 | 47,070 | 85,264 | 159,225 | **177,865** | 177,865 | 177,865 |
| utility loss SPENT, of the budget | 1.0% | 5.0% | 10.0% | 20.0% | 30.0% | **33.1%** | 33.1% | 33.1% |
| cells changed (its stacker) | 7,838 | 33,682 | 68,335 | 105,737 | 106,534 | 106,549 | 106,549 | 106,549 |

**At `capacity` = 0.5 it spends 33.1% of an allowed 50% — because 177,865 moves is
EVERY candidate move that exists.** The constraint stops being the budget and becomes
the inventory. So 0.7 and 0.9 are byte-identical to 0.5: there is nothing left to buy.

**B. THE MARGINAL MOVES ARE ALMOST ALL NO-OPS**, which is the 2026-09-24 finding at
the top end rather than the middle:

| capacity step | extra moves | extra cells changed | cells per move |
|---|---|---|---|
| 0.05 -> 0.1 | +12,730 | +34,653 | 2.72 |
| 0.1 -> 0.2 | +38,194 | +37,402 | 0.98 |
| **0.2 -> 0.3** | **+73,961** | **+797** | **0.011** |
| **0.3 -> 0.5** | **+18,640** | **+15** | **0.0008** |

The last 92,601 moves bought 812 changed cells between them. So the knapsack's upper
range is spent almost entirely on moves that resample the target's own alleles back.

**C. THE CHANGED-CELL SETS ARE PERFECTLY NESTED.** Recomputed independently from each
run's `new_haplotypes.npy` against the original column:

| pair | lower is a SUBSET of higher | kept | dropped | new |
|---|---|---|---|---|
| 0.01 -> 0.05 | **True** | 7,838 | **0** | 25,860 |
| 0.05 -> 0.1 | **True** | 33,698 | **0** | 34,816 |
| 0.1 -> 0.2 | **True** | 68,514 | **0** | 38,693 |
| 0.2 -> 0.3 | **True** | 107,207 | **0** | 805 |
| 0.3 -> 0.5 | **True** | 108,012 | **0** | 15 |
| 0.5 -> 0.7 -> 0.9 | **True** | 108,027 | **0** | **0** |

**Not one cell is ever un-changed by raising the budget**, and 100% of the cells
changed at `capacity` 0.01 are still changed at 0.9. The knob is purely ADDITIVE: it
does not relocate edits, it only accretes them. That follows from the mechanism --
`capacity` is a budget on a fixed item list and the LP is monotone in it, with the same
seed fixing each block's resampled value.

**D. THE 16.41% IS OF ALLELE CELLS, NOT OF VARIANTS.** At the ceiling the release
differs from the original at **108,027 cells**, which fall on **86,578 distinct sites of
340,824 = 25.40%** of the variants. So "16.41% changed" and "a quarter of the variants
touched" are both true and describe different denominators; quoting one as the other
inflates or deflates by 1.5x.

**Provenance.** `logs/pmsweep.log` (sentinel 0) for moves and loss; the nesting and the
site count recomputed here from each run's `new_haplotypes.npy` against row 0 of
`pangenome.npy` (HG00438) · `external/PanMixer` @ `c182c38`, unmodified · seed 123.

⚠ **AN UNEXPLAINED 1.39% DISCREPANCY, recorded rather than smoothed.** My independent
recount gives **108,027** changed cells at the ceiling against the stacker's logged
**106,549** — a difference of 1,478. I count `(new != orig) & (orig != -1)` over the
whole site axis; their counter presumably restricts to some mask. The NESTING result is
computed entirely from my own consistent counts and is unaffected, but the absolute
figures should be quoted from one source or the other, never mixed.

### Two corrections to the entry above this one

⚠ **"6 units with no occupiable donor" was 6 on PANMIXER'S block set and is 4 on ours.**
The 6 came from a feasibility check run before our blocks existed, and then survived as a
HAND-TYPED literal in `scripts/compare_models_chr21.py` -- in an emitted table, which is
the one thing emitting tables is supposed to prevent. With our 8,773-block set it is **4
units, costing 1,222 chain positions**: 534, 257, 246 and 185 variants, at 33.3 Mb,
44.4 Mb, 45.2 Mb and 46.2 Mb. Fixed to compute from the data.

**And the reason only the giants fail is exponential.** The atomic requirement needs ONE
donor called at EVERY position of the unit. At about 4% missingness that is roughly
`0.96^n` per donor for a unit of `n` variants: at n = 20 about 44% of donors qualify, at
n = 100 about 1.7%, and at n = 534 about 2e-10 — so for the largest units it is only
clustered missingness that leaves any donor at all, and for four of them none does.

⚠ **Model C never had the start-to-start bug.** Asked directly: the fix applies to model
B only, and model C could not be affected, because a one-position unit has
`last == first`, so the inter-unit gap and start-to-start spacing are THE SAME NUMBER.
Verified: model C's `switch_mass` is bit-identical to model B's on identity blocks.


## 2026-10-07 — A-IMP-model-b: the LD-block chain, built; and PanMixer's knob, swept

**Goal (Dylan).** Build the model-B fallback so switching happens only between LD
blocks; measure how many variants each architecture changes across `tau`; and sweep
PanMixer's own knob for comparison.

### Model B is now a code path, and model C is its identity case

`hmm.from_arrays(..., blocks=...)` takes a block id per chain position; consecutive
positions sharing an id become ONE ATOMIC unit emitted entirely by one donor.
`blocks=None` is model C. **With identity blocks the two are asserted identical** --
`log_z_p`, `E[differ]`, the released path and `u_achieved` all agree -- so model B is
a regrouping of one mechanism rather than a second mechanism. That assertion is the
test that made the addition safe to make.

**The blocks are OURS.** `scripts/build_blocks_chr21.sh` runs PLINK v1.9 on our staged
30x panel with **the 10 test targets removed** (3,192 samples), giving **8,773** called
haploblocks, longest span 172.7 kb. ⚠ PanMixer's own block set was computed on all
3,202 samples, which CONTAINS our test targets, so its boundaries are weakly a function
of the target -- a standing-constraint-1 problem we avoid by rebuilding. ⚠ RESIDUAL,
the same shape as R4: a first-degree relative still carries half a target's haplotypes
and we do not exclude those. For a real external target the block panel must be
DECLARED IN ADVANCE and must not contain the target.

### Model C against model B, chr21, HG00096, EXACT

| metric | C (one variant) | B (LD blocks) |
|---|---|---|
| `T_blocks` chain units | 305,886 | **59,472** |
| unit size median / p90 / max | 1 / 1 / 1 | 1 / 4 / **1,212** |
| `expected_switches` under the prior | 7,021.0 | **3,394.1** |
| alpha storage | 215 MB | 42 MB |
| longest verbatim donor run | 1 variant | **1,212 variants** |
| units with no occupiable donor | 0 | **4**, dropped (1,222 positions) |
| `E[cells differing]` at `tau` = 0 | 75,177.401 | **73,970.207** |
| floor as a share of 507,692 known cells | 14.8077% | **14.5699%** |
| change from the floor at `tau` = 0.9 | -1.8227 cells | **-2.6114 cells** |
| slope, cells per unit `eta_tau` | -1.2381 | **-1.7738** |
| `sum beta_v` after regrouping | 1.000000000000 | 1.000000000000 |

**Provenance.** `scripts/compare_models_chr21.py`, `logs/A_models_B_vs_C_chr21.md` ·
EXACT from the posterior state marginals, zero Monte-Carlo error · diploid, half
`eta_tau` per haplotype · `chain_span` PROVISIONAL.

**COLUMNS** · *C* is the decided architecture. · *B* is the retained fallback, now
runnable in one argument.

**ROWS** · *`expected_switches` 3,394.1* against the **3,532** predicted independently
on 2026-09-28 -- a 3.9% agreement that is the best available check the implementation
matches the specified model, and the discrepancy that caught a real bug (below). ·
*floor* is BETTER under B, 14.57% against 14.81%: copying a donor across a whole LD
block preserves haplotype structure, so the mosaic resembles a real genome more
closely. · *slope* is **43% steeper** under B. · *1,212 variants verbatim* is
departure (b) at full scale.

**SYNTHESIS.** Model B helps, and it helps by **1.43x on the tilt's effect** -- not by
the 16x a contiguous utility gives. Both the floor and the slope improve, for the same
reason: fewer, longer donor runs mean the utility averages over 3,394 independent
segments instead of 7,021, and variance scales inversely with that count. **So model B
is a mild version of the contiguity lever**, reached by changing the chain rather than
the utility, and it buys about a third of what re-scoping the utility buys. It also
costs: 1,222 chain positions leave the release entirely because no donor is called
across their unit, and the release carries a real donor's haplotype verbatim for up to
1,212 variants, which `tau` does not bound.

⚠ **A BUG THE PREDICTION CAUGHT.** The first implementation spaced units start-to-start
rather than by the inter-unit gap, contradicting its own docstring and
`docs/plan.md`'s model-B spec. Every model-B number was wrong until it was fixed
(floor 14.78% -> 14.57%, slope -1.3355 -> -1.7738). **No test failed**; it surfaced only
because `expected_switches` had an independently predicted value to disagree with.
`docs/bugs.md` has the entry.

### PanMixer's knob, swept -- and it SATURATES

`capacity` is its utility-loss budget as a fraction 0-1 (`obfuscate.py:10`). Swept on
HG00438, chr21, seed 123, its own pinned checkout unmodified:

| `capacity` | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| alleles changed | 7,838 | 33,682 | 68,335 | 105,737 | 106,534 | 106,549 | 106,549 | 106,549 |
| share | 1.20% | 5.18% | 10.50% | 16.29% | **16.41%** | 16.41% | 16.41% | 16.41% |

**Provenance.** `scripts/run_panmixer_capacity_sweep.sh`, `logs/pmsweep.log`, sentinel
0 · `external/PanMixer` @ `c182c38`, unmodified · one run per capacity.

**SYNTHESIS.** **The knob saturates at 16.41% from capacity 0.3 upward** -- four
settings return the byte-identical count -- and is already within 0.12 points of its
ceiling at 0.2. Its entire reachable range is **1.20% to 16.41%**. The 2026-09-18 note
that "the knapsack saturates" is confirmed and sharpened: it saturates at **0.2-0.3**,
not somewhere above 0.5, and the ceiling is a property of how many blocks are
selectable at all, not of the budget.

⚠⚠ **CORRECTION to the 2026-10-03 entry.** That entry put the matched-fidelity
operating point "near capacity 0.45", interpolating from only two measured capacities
(0.1 and 0.5) through a region where the curve is in fact FLAT. With the full sweep our
model-C 14.81% is reached at **capacity ~= 0.17**. Stale copies hunted in
`docs/plan.md` and the status deck.

**And the frontier now exists.** Ours spans 14.8077% to 14.8070% across `tau` = 0 to
0.99 -- a span of 3.28 cells, i.e. a point, not a curve. PanMixer spans 1.20% to 16.41%
across its knob. **Our fixed point sits inside their range**, which is what makes a
matched-fidelity comparison possible at all, and it is also the sharpest statement yet
of R5: one tool has a usable dial on this axis and the other does not.


## 2026-10-03 (later) — the differ-rate at honest n, and a finding I nearly published from n=10

**Dylan caught it.** The differ-rate table reported 76,283 cells at `tau` = 0, 0.05, 0.1,
0.25 and 0.5 — identical to the digit — and he asked whether that was the same random
seed. **It was.** One seed, so one uniform stream; with a per-position tilt of 1.8e-06 the
inverse-CDF draw lands on the same state at essentially every position, so identity is the
EXPECTED outcome, not a finding. ⚠ **This is the SECOND time this error class has been
committed** — `docs/bugs.md`, "A headline finding measured under common random numbers" —
and the first time it was caught by an audit rather than by me.

**A. THE LEVEL, 40 independent seeds (3000-3039), HG00096, diploid, of 507,692
target-known cells:**

| `tau` | 0.0 | 0.5 | 0.9 |
|---|---|---|---|
| mean cells differing | 75,139.9 | 75,151.7 | 75,149.8 |
| sd between draws | 693.2 | 696.2 | 696.5 |
| **differ-rate** | **14.800%** | **14.803%** | **14.802%** |

**B. THE EFFECT, paired within seed — AND THE SIGN IS NOT STABLE.** Pairing (common random
numbers within a seed) is the right design for a small treatment effect, because it cancels
between-draw noise; the earlier error was reporting a paired result as an absolute LEVEL at
n = 1, not the pairing itself.

| paired `differ(tau) - differ(0)` | n = 10 (seeds 2000-2009) | n = 40 (seeds 3000-3039) |
|---|---|---|
| `tau` = 0.5 | **-7.0** (sd 21.1, t = -1.05) | **+11.8** (sd 51.4, t = +1.45) |
| `tau` = 0.9 | **-21.2** (sd 39.2, t = -1.71) | **+9.8** (sd 60.6, t = +1.02) |
| seeds whose draw moved, `tau` = 0.9 | 10/10 | 37/40 |
| state cells moved, `tau` = 0.9 | 457.2 | 542.3 |

⚠⚠ **AT n = 10 THE EFFECT LOOKED CORRECTLY SIGNED AND MONOTONIC IN `tau` (-0.1, -7.0,
-21.2) AND I WAS ABOUT TO PUBLISH THAT AS A FINDING.** At n = 40 the sign FLIPS to +11.8
and +9.8. Neither n reaches significance (|t| from 1.02 to 1.71). **It was noise with a
pattern, and the pattern was the thing I found compelling.** Only running a larger n
caught it — not the reasoning, which had already produced a tidy story about direction and
monotonicity.

**C. THE EXPECTED EFFECT, DERIVED RATHER THAN SAMPLED, and it settles the question.**
`E_Q[u]` moves **2.94e-06** across the whole `tau` grid (exact, three forward passes). A
flipped cell is worth about `mean(beta_v)` = 3.269e-06 of `u`, so:

    cells the tilt should flip  =  2.94e-06 / 3.269e-06  =  0.90 per haplotype
                                =  **1.80 cells of 507,692**  (3.5e-06 of them)

Against a paired sd of 51-61 cells, that needs **3,300 to 4,500 paired draws** for 2 sigma
(8-11 hours). So neither n = 10 nor n = 40 could ever have resolved it, and the sign
flipping between them is exactly what a 1.8-cell effect under 56-cell noise looks like.

⚠⚠ **[CORRECTION - 2026-10-03, within the hour, by independent audit]** Everything in
sections A and B above was measured by a script that **did not sample this mechanism**: it
tilted each haplotype at the FULL `eta_tau` instead of half of it, so every number in them
belongs to a mechanism with twice the intended tilt. And the whole sampling exercise was
unnecessary — `E[differ]` is a LINEAR FUNCTIONAL of the posterior state marginals and so is
EXACT. Both defects are in `docs/bugs.md`. **The exact answer, `scripts/exact_fidelity.py`,
validated against brute-force enumeration to 4.8e-15:**

| `tau` | 0.0 | 0.05 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 0.99 |
|---|---|---|---|---|---|---|---|---|
| `E[differ]` | 75177.401 | 75177.339 | 75177.277 | 75177.085 | 75176.721 | 75176.196 | 75175.578 | 75174.124 |
| shift vs `tau`=0 | 0 | **-0.062** | **-0.124** | **-0.316** | **-0.680** | **-1.205** | **-1.823** | **-3.277** |

**The tilt is REAL, correctly SIGNED, and EXACTLY LINEAR in `eta_tau` at -1.2381 cells per
unit** (range -1.2381 to -1.2381 across all seven points). Baseline **75,177.40 of 507,692
= 14.8077%**. So the sign never was unstable — the SAMPLING could not see a 1.8-cell effect
under 56-cell noise, and the n=10 / n=40 flip was noise around a true value of -0.680 and
-1.823. ⚠ My claim that "the sign is not stable" is therefore also withdrawn: it is the
estimator that was unstable, not the quantity.

**SYNTHESIS, and it corrects my own framing twice over.**

1. **Pairing is genuinely worth a lot**: it cuts the sd from ~695 cells (between-draw) to
   ~56 (paired), a **154x** reduction in draws required. That is real and it is the right
   design for any future effect measurement here.
2. **But it does not rescue the finding, because detectability was never the binding
   problem.** The effect is **1.8 cells of 507,692**, and we have known its size exactly
   since 2026-09-30 from the derivative of `log Z` — no sampling needed. My "6,552,960
   draws" headline invited exactly the wrong follow-up question ("can we get n down?")
   when the answer is that the magnitude, not the measurability, is what makes it
   negligible. The very first write-up said "the effect is exactly known and negligible",
   and every later framing drifted away from that.
3. **The differ-RATE Dylan asked for is 14.80%**, and it is flat to three decimal places
   from `tau` = 0 to 0.9. That is the number for the head-to-head axis, and the fact that
   it does not move with `tau` is the honest content of R5.

⚠ Cosmetic bug in the measurement script's output, recorded so the log is not misread: the
"seeds moved" column prints `/10` hard-coded regardless of n, so the n = 40 rows read
"27/10" and "37/10" where they mean 27/40 and 37/40. The counts themselves are right.


## 2026-10-03 — `target_fidelity` as the comparable axis: our differ-rate measured

[CORRECTION - 2026-10-10]: **every `15.03%` in this entry is a SINGLE-DRAW SAMPLING
ESTIMATE and is superseded for quoting by the exact `14.8077%`.** Both are the same
statistic on the same denominator — cells differing from the target, over 507,692
target-known cells, HG00096 diploid chr21 — but 15.03% is 76,283 from one draw at seed
1000, while 14.8077% is 75,177.401 computed exactly from the posterior marginals with
zero Monte-Carlo error (`logs/A_exact_fidelity_HG00096.md`, 2026-10-03; `docs/data.md`
records that the exact value supersedes every sampling estimate of this quantity). The
sampled draw ran **1,106 cells high**, about 1.5% relative.
**The entry's conclusions all survive the substitution**, which is why it is corrected
here rather than rewritten: ours 14.81% against theirs 10.50-11.70% on a matched region
is still about 1.3-1.4x; PanMixer's capacity dial still SPANS our fixed point; and the
flatness in `tau` is confirmed, not weakened — the exact movement across the whole grid
is **3.28 cells**, which no sampled table of this size could have resolved.
⚠ The affected lines below are the two measured tables (which are honest records of what
those draws gave, and are left verbatim) and four CLAIM lines that should be read at
14.8077%: the matched-comparison sentence, the synthesis, the three-way share table, and
the "superficially similar" sentence. The `capacity ~= 0.17` crossing quoted elsewhere was
already computed against the exact 14.81% and needed no change.

**Correction upstream.** The 2026-10-03 entry below this one first declared the
altered-cell count non-comparable. Corrected in place: it is `target_fidelity`, the one
axis that transfers between arms, and Dylan was right to push back.

**A. Our differ-rate barely moves with `tau`** — HG00096, chr21, seed 1000, diploid:

| `tau` | 0.0 | 0.05 | 0.10 | 0.25 | 0.50 | 0.90 |
|---|---|---|---|---|---|---|
| cells differing from the target | 76,283 | 76,283 | 76,283 | 76,283 | 76,283 | 76,269 |
| of 507,692 target-known cells | 15.03% | 15.03% | 15.03% | 15.03% | 15.03% | **15.02%** |

⚠ **15.03% AT `tau` = 0.** A pure prior draw, using zero target information, already
coincides with the target at **84.97%** of its called positions — because two humans share
most of their alleles. Our differ-rate is therefore a FLOOR set by cohort-to-target
mismatch, not a dial the mechanism turns. This is the `target_fidelity` free-lunch warning
measured at the cell level rather than argued.

**B. Target choice moves it by less than 2 points** — `tau` = 0.5, seed 1000:

| target | HG02389 (EAS) | HG00097 (EUR) | HG00096 (EUR) | HG04211 (SAS) | NA19350 (AFR) | NA20294 (AFR) |
|---|---|---|---|---|---|---|
| differ-rate | 14.83% | 14.90% | 15.03% | 15.09% | 16.43% | **16.63%** |

Range **1.80 points**, and ordered as population genetics predicts — the two AFR targets are
highest against a multi-continental cohort. **So comparing across different targets is NOT
what invalidates a cross-arm count:** the target effect is small beside the 4.5-point gap to
PanMixer.

**C. The three real obstacles, with the biggest one quantified.**

1. **Denominator mismatch — MEASURED.** **34,937 of PanMixer's 340,824 records (10.25%)**
   fall outside our `chain_span` (34,345 below, 592 above). That region carries 36.76%
   cohort missingness against 6.57% elsewhere, and missing genotypes are released VERBATIM,
   so it is low-change territory diluting their rate downward. Bounding it: if they change
   nothing there their matched-region rate is 10.50%/0.8975 = **11.70%**; if they change it
   at their average rate it stays **10.50%**. So the matched comparison is ours **15.03%**
   against theirs **10.50-11.70%** — we alter about **1.3-1.4x** as many cells.
2. **Different VCFs.** Ours is vg 1.68, theirs vg 1.36: 99.55% of sites agree by POS but
   only 92.7% by (POS, REF, ALT), and ALT ordering differs at multi-allelic sites, so
   allele INDICES are not comparable cell-by-cell there. The defensible common subset is the
   biallelic SNPs, where ordering is unambiguous: **274,458 of their records (80.5%)**, of
   which **249,750** sit inside our `chain_span`. That intersection is the site set a
   rigorous count should run on.
3. **Their obfuscated output was not kept.** The 68,335 figure is from `logs/pm_obf.log`,
   not from a saved artifact, so a matched-region recount needs their obfuscation re-run
   (~68 s). Cheap, but it must run on the SAME VCF for the count to mean anything, which is
   the condition `docs/plan.md` already states for the head-to-head.

**SYNTHESIS, and it answers the question as asked.** Yes, the two tools alter a similar
NUMBER of cells — ours 15.03%, theirs 10.50% raw and 10.50-11.70% on a matched region, so
within about 1.4x. But the two numbers are the same statistic of different operations, and
the difference shows up in how they respond to the knob. **PanMixer's capacity dial SPANS
our fixed point**: at capacity 0.1 it alters fewer cells than we do (10.50%), at capacity
0.5 it alters more (16.41%), while ours sits at 15.0% for every `tau` from 0 to 0.9. Theirs
is a budget of deliberate edits to a genome it starts from; ours is the residual mismatch of
a mosaic that never saw the target's data. **So the comparable axis is comparable, and what
it shows is that a matched-fidelity comparison is achievable — PanMixer at a capacity near
0.45 would sit at our differ-rate — which is exactly the alignment the head-to-head needs.**
That is a concrete, cheap experiment and it is now the obvious first step for
`[A-EVL-headtohead]`.


## 2026-10-03 — How many alleles each tool actually changes on chr21

**Question (Dylan).** How many variants does PanMixer actually change?

**PanMixer's figures are already measured** (2026-09-18 run, HG00438, chr21, seed 123)
and are restated here only so they sit beside ours. **Ours is new.**

| descriptive count | PanMixer, capacity 0.1 | PanMixer, capacity 0.5 | ours, `tau` = 0.5 |
|---|---|---|---|
| selected moves | 47,070 | 177,865 | n/a — no selection step exists |
| of those, moves changing NOTHING | 14,526 (30.9%) | not measured | n/a |
| allele cells altered vs the target | **68,335** | **106,549** | **76,283** |
| denominator | 650,541 | 649,242 | 507,692 target-known cells of 611,772 released |
| share | **10.50%** | **16.41%** | **15.03%** |
| cells released as the target's OWN data | 582,206 (89.50%) | 542,693 (83.59%) | **0** |
| sequence rewritten | 3.13 Mb (6.7% of chr21) | not measured | n/a — no edit, the path is built from the cohort |

**Provenance.** PanMixer: `logs/pm_obf.log`, sentinel 0, `external/PanMixer` @ `c182c38`,
HG00438 (a COHORT MEMBER), capacities 0.1 and 0.5, seed 123. Ours: HG00096 (EXTERNAL),
`tau` = 0.5, seed 1000, diploid, `indicator/w=inv_support/beta=chain`, `chain_span`
provisional, T = 305,886.

[INCORRECT] - ⚠⚠ **THESE COLUMNS ARE NOT A COMPARISON AND MUST NEVER BE QUOTED AS ONE.** Different targets, different denominators, and `capacity` and `tau` are formally incomparable (`docs/terms.md`). They are descriptive counts of two different operations, put side by side to answer one factual question. Nothing here says either tool is better.
[CORRECTION - 2026-10-03]: **OVER-STRONG, AND IT BLOCKED THE RIGHT COMPARISON** (Dylan).
The altered-cell row IS comparable, and it is already a named term: `target_fidelity`,
which `docs/terms.md` records as computable mechanism-agnostically on either arm's output.
It is the axis `docs/plan.md` `[A-EVL-headtohead]` lists as a MUST. What is incomparable is
the OPERATING POINT — reading a privacy claim off `capacity` 0.1 beside `tau` 0.5 — not the
quantity. The correct use is as the x-axis of a frontier, with the arms aligned on measured
attack success, exactly as the head-to-head already specifies. The three real obstacles are
below, and they are about denominators and VCFs, not about the quantity being meaningless.

**COLUMNS** · *PanMixer 0.1 / 0.5* — its two reported operating points. · *ours* — one
release at the mid-grid `tau`.

**ROWS** · *selected moves* — PanMixer picks a SUBSET of blocks to rewrite; we have no
such step. · *allele cells altered* — the counts land within a factor of 1.4 of each other,
which is a coincidence of two unrelated mechanisms and nothing more. · *cells released as
the target's own data* — **this is the only row that matters**, and it is the whole
difference between the two threat models.

**SYNTHESIS, and it is the row nobody would think to put in the table.** The altered-cell
counts are superficially similar — 10.50% and 16.41% against our 15.03% — and that
similarity is meaningless. PanMixer starts from the target's own genotype column and
rewrites a selected subset, so **89.50% of its release at capacity 0.1 is the target's
REAL DATA, verbatim.** We never release the target's data at all: every one of the 611,772
cells is drawn from the cohort, and the 431,409 that coincide with the target do so
because two humans share most of their alleles, not because anything was copied from them.

**That single structural fact is why the two guarantees differ by construction, and it is
worth having in one sentence for the writeup.** Because PanMixer releases unselected
blocks verbatim, two inputs differing in an unselected block produce outputs with DISJOINT
supports, so TV between their release distributions is exactly 1 and no bound of the kind
Theorem 1 states is available at any capacity. Our support `Y_D` is a function of `G` and
`D` alone, identical for every target, so every pair of inputs is bounded. The altered-cell
count is a measure of how much a tool EDITS; it is not a measure of what it LEAKS.

⚠ Also note what the 30.9% no-op figure means in this light: 14,526 of PanMixer's 47,070
selected moves resampled the target's own alleles back, so they are counted as privacy
work while leaving the genome untouched — and those blocks are then indistinguishable from
unselected ones in the output.


## 2026-10-01 — R5 knob inventory: every degree of freedom in the tilt, measured

**Question (Dylan).** Can `beta_v` be made larger by changing something else, keeping the
maths valid?

**Frame.** Detectability of the tilt in one draw is `gain/noise = eta_tau * sd(u)`, and
`Var(u) ~ Var(local match rate) / (independent segments the utility spans)`. So there are
exactly three factors: `eta_tau`, the local match variance, and the segment count. Every
knob below is one of those three.

| knob | effect on gain/noise | legitimate? | cost |
|---|---|---|---|
| `phi_v` functional form | **1.00x — no room at all** | yes | a 0/1 variable already MAXIMISES variance for any [0,1] variable of given mean |
| reward matches AND penalise mismatches | **1.00x — literally the same mechanism** | yes | none, and no benefit: measured max abs difference between the two release distributions is **2.776e-17** |
| negative `beta_v` summing to 1 | **< 1x — it COSTS** | yes but perverse | the range is `sum |beta_v|`, not `sum beta_v`; one weight at -0.3 gives range **2.0534**, so `eta_tau` must SHRINK to 0.2675 |
| release ONE haplotype instead of a pair | **sqrt(2) = 1.41x** | yes | half the release object |
| cohort heterozygosity weights | **2.22x** | yes — function of D alone | `E_Q[u]` 0.695 -> 0.520; a DIFFERENT utility |
| contiguous 1,000-position utility | **16.16x** | yes | target information only in the window |
| contiguous 100 + heterozygosity | **86.3x** | yes | both of the above, compounded |
| `tau` 0.5 -> 0.9 | **2.68x** | a weaker privacy claim | `eta_tau` enters linearly |
| fewer donor segments via `Ne*r` | **1.36x, and NON-MONOTONIC** | not target-dependent, but see cost | the release becomes a long verbatim copy of ONE real donor — departure (b), and `tau` protects the TARGET, not the donors |

**Provenance.** `Var_Q(u)` from `d2 log Z_p / d eta2`, three forward passes per row, no
sampling · HG00096 haplotype 0 · chr21 · `tau` = 0.5 unless stated · the equivalence and
negative-weight rows are exact on the toy.

**COLUMNS** · *effect* is on `eta_tau * sd(u)`, which is dimensionless and comparable
across utilities even though the utilities themselves are not (constraint 7). ·
*legitimate* means it does not make `R_D`, `Y_D` or `beta_v` a function of the target.

**ROWS** · the first three are the answer to the question as asked, and all three are
**no**. · the middle five are the real levers. · the last one is the surprise.

**SYNTHESIS, and three things worth keeping.**

**1. The two formulations of the utility are the same mechanism.** Rewarding matches
(`phi` valued 0 or 1, range 1, `eta_tau`) and rewarding matches while penalising mismatches
(`phi'` = 2`phi` - 1 in -1 or +1, range 2, `eta_tau`/2) give the IDENTICAL release
distribution, to 2.776e-17. Algebraically `(eta/2)(2u - 1) = eta*u - eta/2` and the
constant is absorbed by `Z_p`. The tilt sees only DIFFERENCES in `u` between outputs, so
`u`'s zero point and scale are both free and both cancel. Any reformulation that merely
rescales or shifts `u` is not a new mechanism.

**2. ⚠ CORRECTION TO A NUMBER QUOTED SINCE 2026-09-30.** The "3,276,480 draws" figure is
the HAPLOID one. The diploid release — the only kind we have ever run — has
`Var(u_dip) = V/2` and gain `(eta_tau/2) * V`, so `gain/noise = (eta_tau/sqrt2) * sqrt(V)`
= **7.8129e-04** and **6,552,960 draws**, twice what was reported. Releasing a single
haplotype therefore buys exactly **sqrt(2)**, not the 2x a first look at half `eta_tau`
suggests. Stale copies hunted in `docs/plan.md` and `docs/terms.md`.

**3. LENGTHENING THE SEGMENTS BARELY WORKS, AND SATURATES.** Dividing the recombination
rate by 10, 100 and 1,000 takes expected switches from 7,021 to 854 to 87.8 to 8.8, and
raises `Var` only 1.6x, 1.8x, 1.6x — it PEAKS around 100 switches and then falls. The
reason is instructive: with one donor spanning the chromosome, `Var(u)` becomes the
between-donor spread in GLOBAL match rate, and every donor matches a human target at about
the same rate, so that spread is tiny. **So the contiguity lever does not work by
lengthening segments.** It works because a utility over a small window measures a LOCAL
match rate, and local haplotype sharing genuinely varies between donors while the global
average does not. Shrinking the utility's footprint helps on both factors at once;
stretching the segments helps on only one and runs into a floor.

**The best legitimate stack** is a haploid release with a heterozygosity-weighted
contiguous 100-position utility at `tau` = 0.9: gain/noise **0.2463** against 7.8129e-04
as built, a factor of **315x**, which is **66 draws** instead of 6,552,960 — about
99,000x fewer. The two working levers compose almost exactly multiplicatively (16.16 x 2.22
= 35.85 predicted, 36.66 measured). ⚠ Every one of those configurations is a DIFFERENT
utility function and a different release object; this table says what is reachable, not
what we should pre-register.

⚠ **ONE TRAP NAMED.** Weighting `beta_v` by `p_v(1 - p_v)`, the per-position match variance,
would be the direct way to raise the diagonal — and `p_v` is a function of the TARGET, so it
is a standing-constraint-1 violation that voids Theorem 1. Cohort heterozygosity is the
legitimate proxy and is what was measured above; it is a function of D alone.


## 2026-09-30 (later still) — R5 REOPENED: contiguity is the lever, and there is an exact ceiling

**Why this entry exists.** Two hours earlier the same day, R5's option list was narrowed to
"nothing works" on the strength of a table that concentrated `beta_v` on the top-n positions
by `w_v`. **Those positions are SCATTERED across the chromosome**, and that confound was
never separated. It is the whole answer.

**Measured, sampling-free, `Var_Q(u)` from `d2 log Z_p / d eta2` at `tau` = 0.5, HG00096
haplotype 0, chr21. Same n, same `sum beta_v` = 1, only the LAYOUT differs.**

| n positions | `Var_Q(u)` scattered | `Var_Q(u)` contiguous | ratio |
|---|---|---|---|
| 100,000 | 1.1404e-05 | 1.5149e-05 | 1.3x |
| 10,000 | 1.1140e-05 | 1.6422e-04 | **14.7x** |
| 1,000 | 1.1132e-05 | 1.0565e-03 | **94.9x** |
| 100 | 6.6763e-05 | 8.3509e-03 | **125.1x** |

**Provenance.** `src/mechanism/ffbs.py` forward passes only, no sampling; three passes per
row. Contiguous windows are centred on the chain; scattered sets are the top n by
`w_v = 1/support(v)`.

**COLUMNS** · *scattered* is what the earlier table measured. · *contiguous* is a single
window of the same size. · *ratio* is the finding.

**ROWS** · *100,000* — a third of the chromosome, so both layouts span most of the segments
and barely differ. · *10,000 down to 100* — scattered `Var` is PINNED near 1.1e-05 however
few positions it covers, while contiguous `Var` rises roughly as 1/n.

**SYNTHESIS, and it corrects this project twice over.** `Var_R(u)` is governed by **how many
independent mosaic SEGMENTS the utility spans**, not by `sum beta_v^2`. Scattered positions
each sit in their own segment, so concentrating onto them leaves the segment count — and
hence the variance — untouched; that is why the earlier table read as a flat line and was
mistaken for a dead end. A contiguous window spans few segments, so `u` approaches "does the
drawn donor match the target through this window", whose variance approaches `p(1-p)`.

**AND THERE IS AN EXACT CEILING, which settles what is possible rather than what we tried.**
The tilt's gain is `eta_tau * Var_Q(u)` and a single draw's noise is `sd(u) = sqrt(Var)`, so

    gain / noise per draw  =  eta_tau * sqrt(Var)  <=  eta_tau / 2

because `u` in [0,1] bounds `Var` at 1/4. That is a bound on EVERY utility of this form, at
every `beta`, forever:

| `tau` | 0.10 | 0.25 | 0.50 | 0.90 | 0.99 |
|---|---|---|---|---|---|
| `eta_tau` | 0.1003 | 0.2554 | 0.5493 | 1.4722 | 2.6467 |
| ceiling on gain/noise per draw | 0.0502 | 0.1277 | 0.2747 | 0.7361 | 1.3233 |
| draws for 2 sigma at the ceiling | 1,589 | 245 | **53** | 7 | 2 |

**So the mechanism's effect is ALWAYS sub-noise in a single release** — that is what a TV
bound of `tau` means, and it is not a defect. What varies is how many draws it takes to SEE
it. At our current configuration `sd(u)` = 2.0115e-03 gives 1.1049e-03 sigma per draw and
**3,276,480 draws** for 2 sigma. A contiguous window of 100 gives 5.0197e-02 sigma per draw
and about **1,600 draws** — 4 hours at 9.4 s per diploid release. The theoretical best at
`tau` = 0.5 is **53 draws**. We are sitting **250x below the achievable ceiling**, not at it.

⚠ **AND R5's OPTION 2 WAS ALSO WRONG.** It dismissed "shrink the scope `tau` covers" as
"k independent releases of one genome, and they COMPOSE, which standing constraint 4
forbids". That conflates two different things. Releasing k windows each with its own `tau`
does compose and is forbidden. Defining ONE release whose utility sums to 1 over ONE window
does NOT compose — the chain still runs genome-wide, one path is released, and the positions
outside the window are simply a pure prior draw. It is a narrower UTILITY, not a second
release. That distinction is the door the earlier analysis closed by mistake.

**What this does NOT settle, and it is the real question for the PI.** A utility over a
contiguous window means the release preserves target information THERE and is a pure prior
draw everywhere else. That is a different product from whole-chromosome fidelity, and
choosing between them is a question about what the release is FOR, not about `beta_v`. The
measured trade is now on the table instead of being hidden behind a flat line:
whole-chromosome fidelity is unmeasurable at any feasible n; regional fidelity is measurable
at roughly 1,600 draws. ⚠ Both are PRE-REGISTRATION items under standing constraint 5.


## 2026-09-30 (later) — A-FIX-mechanism-audit: an independent audit of the same day's build

**Hypothesis.** A test suite written by the agent that wrote the code tests that agent's
MODEL of the code, not the spec. 88/88 passing should not be taken as evidence of
completeness.

**Method.** Six auditors, one per spec dimension (`baseline_model`, `utility`, `ffbs`, the
toy gate, the eight standing constraints, and an INVERSE audit of the spec against the
code). Each was given `docs/` as the authority and told explicitly that the docstrings are
the claim under audit, not evidence for it — because the same agent wrote the confident
prose and the tests that check it. Every gap rated blocker or material went to three
refuters with distinct lenses (did the auditor misread the SPEC; did it misread the CODE; is
the deviation ALREADY DISCLOSED), and survived only on >= 2 of 3 failing to refute. Then a
completeness critic, asked what the six dimensions never opened.
**Provenance.** 127 agents, 0 errors, 41 minutes. 30 gaps survived, 10 were refuted, 29
minor/known were recorded without verification.

**RESULT: six real defects, all of which the 88-test suite passed over.** Fixed, with the
first four verified by MUTATION before the fix was accepted — the lesson of this project's
own `ratio_spread` finding, applied to the fixes rather than the gates. Full entries in
`docs/bugs.md`; in brief:

| defect | how it would have shown up |
|---|---|
| `epsilon_tau` returned `2*arctanh(tau)/delta_u` | at `delta_u` = 2, exactly `arctanh(tau)` — the expression `terms.md` forbids in capitals |
| diploid `log_z_p` never bracketed | a doubled diploid exponent ACCEPTED at every `tau`; 0.763 published against a stated [0, 0.549] |
| diploid seeds were `seed` and `seed + 1` | one release's haplotype 1 bit-identical to the next release's haplotype 0, 0 of 305,886 differing |
| `release()` took the alleles as a free argument | a target-pruned matrix of the same shape accepted — the constraint-1 shape the directory split cannot catch |
| `phi_v` never called by the mechanism | the documented, tested function on no live path; the indicator re-inlined twice with different null guards |
| `sum_v beta_v^2` reported as `mean(beta_v)` | 3.269192e-06 against the true 4.556723e-06; the two coincide only for uniform weights |

**AND IT OVERTURNED THE DAY'S HEADLINE.** Corrections are in place above. The claim "the
released path is bit-identical at `tau` = 0, 0.1, 0.5 and 0.9 — 0 of 611,772 positions
differ" was measured on ONE seed under COMMON RANDOM NUMBERS, and at ten disjoint-stream
draws **8 of 10 move**, by a mean of 417 state positions per haplotype. The critic's sharpest
observation is one no auditor made: **the docs PREDICTED that result two paragraphs before
reporting it as a discovery** ("the expected number of positions the tilt moves is of order
`eta_tau`"), so the measurement could not have surprised anyone and should never have been
the headline row. The substance survives untouched, because it rests on the exact derivative
measurement instead.

**AND IT CLOSED R5's LEADING OPTION.** "The lever is `sum beta^2`" was quantitatively
unsupported; measured, concentrating `beta_v` 3,000-fold buys 16.5x the variance and stays
55x below the noise, because **94.8% of `Var_R(u)` is between-position covariance** from the
mosaic's segment structure. What is left is structural and is the thing for the PI: exact
`ffbs` sampling REQUIRES `phi_v` to decompose per position, and that requirement is what
forces `u` to be a low-variance average over T positions.

**WHAT THE CRITIC NAMED AS THE SINGLE MOST IMPORTANT THING STILL MISSING, and it is right.**
**The diploid release path has no correctness gate, and it is the only path that has ever
been run.** All 18 releases in the first chr21 run carry `ploidy: 2`. Every falsification
gate in the project — exact TV, `ratio_spread`, the sharp pair, the state-path joint,
clause (f), the KILL CRITERION — is HAPLOID. And it cannot simply be written: `exact.py` has
no diploid primitive, so closing it needs a diploid `enumerate_baseline` and
`release_distribution` over the product support first. The bracket hole is patched and
regression-tested, but a patched assertion is not a gate.

**Lessons worth keeping.**
1. **A row whose columns are all identical deserves suspicion before it deserves a synthesis
   paragraph.** `utility_retained` printing 0.6940 at every `tau` from 0 to 0.9 was the
   tell, and it was written up as the finding instead.
2. **A rounded MEAN can conceal the movement its own per-draw records show.**
   `donor_switches` printed 14,169 at every `tau` while the records held
   14,147/14,211/14,147 at `tau` = 0.9 against 14,148/14,212/14,148 below it.
3. **Three transcriptions of a table are not a check on it.** `calibration.py` claimed the
   docs contract verified it against `docs/terms.md`; that test recomputes from `math.atanh`
   and never imports the module, and `test_mechanism.py` compares against a third hardcoded
   copy. No edge to the document existed, which is how the `epsilon_tau` defect survived. The
   replacement PARSES the table out of the file.
4. **Two implementations of one definition, with no test between them, is the same defect
   this project criticises upstream** — and we shipped three of them in a day (`phi_v` vs
   two inlined copies; `exact.utility_vector` vs `utility.u_of_output`; the calibration
   table thrice).


## 2026-09-30 — A-THY-toy-enumeration, A-IMP-cohort-hmm, A-IMP-utility, A-IMP-ffbs-sampler

**THE MECHANISM EXISTS AND IS GATED.** `src/mechanism/` holds `baseline_model`, the
bounded utility and an exact `ffbs` sampler; `A-THY-toy-enumeration` passes all its
clauses plus a sixth we had to add. 46 new tests, 87 in the suite, all passing. The
first RanPanMixer numbers exist and they are reported below with a result nobody
expected.

**Order of work, deliberately.** The toy and the enumeration reference were built
FIRST, before any sampler, so the sampler was written against a checkable target.
That is possible because clauses (b), (c) and (d) of the gate need only the
DEFINITIONS of `R_D` and `u_path` — not an implementation of them — so the toy can
falsify Theorem 1 with no sampler in existence. Only clause (a) needs `ffbs`.

### The toy: K = 4, T = 7, reads nothing

16,384 state paths -> **2,304 with positive mass** -> **96 output paths in `Y_D`**.
Hand-written in `src/mechanism/toy.py`: no graph, no VCF, no cohort artifact, no
PanMixer path, no `data/`. Every feature earns its place and `toy_cohort()` asserts
each one at construction, because a toy that has quietly stopped exercising a
feature still reports green: two donors differing at EVERY position, a donor that is
another donor minus one cell (so two state paths collapse to one output path), a
missing donor allele, a three-allele position, a zero-cM step, a 0.02 cM step, and
non-uniform support.

### ⚠ THE PRE-REGISTERED TV GATE HAS NO POWER. We added a sixth clause.

Exit-gate clause (b) is `tv_empirical <= tau` by exact enumeration. It holds — and it
is nearly useless. Measured worst exact TV over all 3,571,128 pairs of admissible
inputs per `tau`:

| `tau` | 0.01 | 0.05 | 0.10 | 0.25 | 0.50 | 0.90 |
|---|---|---|---|---|---|---|
| worst exact TV | 0.002738 | 0.013706 | 0.027488 | 0.069974 | 0.150024 | 0.387987 |
| as a share of `tau` | 27.4% | 27.4% | 27.5% | 28.0% | 30.0% | 43.1% |

The bound is never approached, so it cannot detect a mis-scaled tilt. Demonstrated by
MUTATION: apply the tilt with `epsilon_tau` in place of `eta_tau` — twice the intended
exponent, a one-token error, and by `CLAUDE.md` constraint 2 the highest-severity error
class in the project — and then ask which gates fire.

| gate | `tau` 0.01 | 0.05 | 0.10 | 0.25 | 0.50 | 0.90 |
|---|---|---|---|---|---|---|
| `tv_empirical <= tau` * | PASSES | PASSES | PASSES | PASSES | PASSES | PASSES |
| `log_z_p` in [0, `eta_tau`] * | PASSES | PASSES | PASSES | PASSES | **fails** | **fails** |
| `ratio_spread <= exp(epsilon_tau)` * | **fails** | **fails** | **fails** | **fails** | **fails** | **fails** |

**Provenance.** `tests/test_toy_enumeration.py`, K = 4, T = 7, 96 output paths, exact
enumeration; 2,916 admissible inputs deduped to 2,673 distinct release distributions
(exact, since `Q_p` depends on `p` only through its utility vector).

- *`tv_empirical`* — a real bound and a bad test. Broken mechanism passes at every `tau`.
- *`log_z_p`* — blind below `tau` = 0.5, because `log Z_p` is about `eta_used * E_R[u]`
  and `E_R[u]` = 0.455 on the toy, so a doubled exponent still lands at 0.91 of the
  bracket. The defensible `tau` range is exactly where it cannot see.
- *`ratio_spread`* — fires everywhere, because the normaliser CANCELS out of it and a
  mis-scaled exponent has nothing to hide behind.

**SYNTHESIS.** `ratio_spread` — `max_y Q_p(y)/Q_q(y)` over `min_y` of the same ratio — is
now a PRIMARY gate and is defined in `docs/terms.md`. Theory: the ratio is
`exp(eta_tau (u_p - u_q)) * Z_q/Z_p`, the `Z` factor is common to every `y` and cancels,
leaving `exp(eta_tau * range_y(u_p - u_q)) <= exp(2 eta_tau) = exp(epsilon_tau)`. It is
SHARP: `toy.sharp_pair()` constructs the input pair that saturates it — D1's alleles
against D2's, which differ at every position, both reachable — and the measured spread
equals `exp(epsilon_tau)` to 1e-9 at every `tau` on the grid (1.020202, 1.105263,
1.222222, 1.666667, 3.000000, 19.000000). A gate sitting exactly on its boundary detects
any inflation of the exponent immediately. The lesson generalises past this one error: we
had a gate list where the headline clause was the weakest member, and only a mutation
test revealed it. `test_a_doubled_tilt_exponent_is_caught_by_exactly_one_gate` asserts
the whole table above, so if any cell ever changes, this entry is what has gone stale.

### Finding — RENORMALISE cannot be a renormalised transition row

`missing_policy` = RENORMALISE says "drop that donor from the state distribution at that
position", which invites masking row `j` of the transition to the allowed set and
renormalising. **That is not well defined.** Across a step with `Delta_x = 0` the
transition is the identity, so a donor available at `t-1` and unavailable at `t` has a
masked row summing to **exactly zero** — no chain to renormalise, and a division by zero.
Measured on our chr21 build: **1,058 zero-distance steps inside `chain_span`** against
4.26% mean missingness, so this is reachable on the real data, not a contrivance.

So the mask enters as a **0/1 potential** and `R_D` is the Li-Stephens chain CONDITIONED
on never occupying an unavailable donor — which is also exactly the form `docs/plan.md`
had already pre-registered, with no row renormalisation in it. ⚠ **The consequence is a
trap.** The conditioning constant is part of `sum_t log c_t`, so `log_z_p` is that sum
MINUS the baseline's own. Measured: `baseline_log_norm` = **-584.8135** on chr21, against
a bracket of [0, 0.5493] at `tau` = 0.5 — out by a factor of about a thousand. The
tempting "fix" of widening the bracket would have discarded the cheapest correctness
check in the project. It is now `baseline_log_norm` in `docs/terms.md`, is
target-independent by construction, and is computed once per cohort and cached.

### Cross-checks that could have failed and did not

- **The O(TK) rank-one forward pass against an O(TK^2) dense reference.** Max absolute
  alpha difference **1.11e-16** on the toy; agrees to 1e-12 at K = 60, T = 40 on random
  data. `exact.py` writes the dense transition from the published closed form and imports
  nothing from `hmm.py` — asserted by a test, because a check that shares code with what
  it checks cannot fail, which is the defect we documented in PanMixer's own vacuous
  shape assertion.
- **`log_z_p` two ways.** The forward pass minus `baseline_log_norm` against
  `log sum_y R_D(y) exp(eta u)` over the enumerated support: agree to **< 1e-12** across
  the whole `tau` grid and 417 inputs.
- **The published closed form against the rank-one decomposition.**
  `p_stay = e^(-d/n) + (1-e^(-d/n))/n` against `(1-s) + s/K`, checked at construction at
  K = 4, 44, 86, 88 to 1e-15. Constants re-derived from the Methods (Ne = 10,000,
  r = 1.26, `d = Delta_x * Ne * r`), never from PanMixer's code.
- **`ffbs` against exact enumeration.** Output distribution at `tau` = 0, 0.25, 0.9 over
  40,000 draws: max per-path z 1.96 / 2.16 / 2.67, TV(empirical, exact) 0.008999 /
  0.008504 / 0.011945 against an analytic expectation of 0.009751 / 0.009732 / 0.009284.
  The state-path JOINT too, not only the output marginals — 2,304 paths, 30,000 draws,
  TV 0.028866 against 0.030701 expected. The joint matters: a sampler drawing each
  position from its posterior marginal independently would pass a marginal check and
  produce mosaics with the wrong recombination structure.
- **`E[donor switches]` under the prior = 7,021.0** on chr21 at K = 88, against the
  **7,155** that the 2026-09-28 model-C decision predicted from an independent
  calculation. The architecture choice is confirmed by the built model.
- **`u_ceiling` = 0.815455** over 253,846 of 305,886 chain positions, independently
  reproducing the **0.815** that deferred decision R1 predicted from weight arithmetic.

### ⚠⚠ THE RESULT: AT chr21 SCALE THE TILT CHANGES NOTHING

First mechanism run. HG00096, diploid, chr21, `beta` scope `chain`, 3 draws per `tau`.
Metrics are rows, `tau` values are columns; the full table is
`logs/A_HG00096_chr21_v1.md`. The rows that matter here:

| metric | `tau`=0 | 0.05 | 0.10 | 0.25 | 0.50 | 0.90 |
|---|---|---|---|---|---|---|
| `utility_retained` * | 0.6940 | 0.6940 | 0.6940 | 0.6940 | 0.6940 | 0.6940 |
| `u_ceiling` * | 0.8155 | 0.8155 | 0.8155 | 0.8155 | 0.8155 | 0.8155 |
| `log_z_p` * | 0.0000 | 0.0348 | 0.0697 | 0.1774 | 0.3815 | 1.0226 |
| `log_z_p` as a share of `eta_tau` | 0 | 69.5% | 69.5% | 69.5% | 69.5% | 69.5% |
| released positions differing from `tau`=0, seed 1000 | 0 | 0 | 0 | 0 | 0 | **0** |
| `donor_switches` (both haplotypes) | 14169 | 14169 | 14169 | 14169 | 14169 | 14169 |
| `seconds_per_release` | 9.4 | 9.5 | 9.4 | 9.8 | 9.4 | 9.4 |

**Provenance.** `chr21/hprc-v1.0-pggb/K88/T305886` · 44 individuals / 88 haplotypes ·
`indicator/w=inv_support/beta=chain`, `beta` over 305,886/305,886 chain positions ·
3 draws per `tau`, seeds 1000-1002 · `baseline_log_norm` -584.8135 ·
`chain_span` PROVISIONAL · deferred decisions R1-R4 UNSETTLED · v1 numbers.

**COLUMNS** · *`tau`=0* is the floor, an untilted draw from `baseline_model`, a real
cohort haplotype mosaic carrying zero target information. · *0.05-0.25* is the
defensible range. · *0.90* is far past any privacy claim we would make and is included
because a bound must be checked where it is loose.

**ROWS** · *`utility_retained`* does not move — not "moves little", does not move. ·
*`u_ceiling`* is 0.8155 and not 1, which is deferred decision R1 measured rather than
predicted. · *`log_z_p`* DOES move, exactly proportionally to `eta_tau`, at a constant
69.5% of it — so the tilt is genuinely applied at the right magnitude. ·
[INCORRECT] - *released positions differing* is **0 of 611,772** (305,886 positions x 2 haplotypes) at
every `tau` up to 0.9, for a fixed seed. · *`donor_switches`* likewise identical. ·
[CORRECTION - 2026-09-30]: **FALSE, and it was true only of the single seed checked.**
Re-measured at ten seeds with disjoint streams, haploid, `tau` = 0 against `tau` = 0.9:
**8 of 10 draws move**, by a mean of **417.2 state positions** and **54.4 allele positions**
per haplotype, range 0 to 1,641. Seed 1000 -- the one seed the original check used -- is one
of the two that happen not to move. Two separate errors produced the claim: the paths were
compared under COMMON RANDOM NUMBERS (`_draw_one` opens a fresh generator from the same seed
at every `tau`), so agreement was the expected outcome rather than a finding; and n = 3 was
too small to see an effect the docs themselves predicted at "order `eta_tau`" per draw. Note
also that one flipped position re-seeds the whole downstream path through the backward
recursion, so the movement is 151x the "order `eta_tau`" prediction -- that prediction bounds
the number of INDEPENDENT flips, not the number of differing positions. ⚠ *`donor_switches`*
was likewise NOT identical: the released records show 14,148 / 14,212 / 14,148 at
`tau` <= 0.5 against 14,147 / 14,211 / 14,147 at `tau` = 0.9, and the table's rounded MEAN
(14,169 either way) concealed it. **THE SUBSTANCE OF THE FINDING SURVIVES AND IS UNCHANGED**
-- it rests on the exact derivative measurement below, not on this row: expected utility
moves 2.94e-06 across the whole `tau` grid. The correct statement is that the tilt VISIBLY
RESHUFFLES the released path and the reshuffling is UTILITY-NEUTRAL, which is a different and
more interesting claim than "nothing changes".
*`seconds_per_release`* is 9.4 s for a diploid release, so compute is not the constraint.

**SYNTHESIS.** The mechanism is correct and the guarantee holds, and it holds because
almost nothing is happening. Three exact measurements, all sampling-free, pin why. Using
`d/d_eta log Z_p = E_Q[u]` and `d2/d_eta2 log Z_p = Var_Q(u)` at `tau` = 0.5 on
HG00096 haplotype 1: `E_Q[u]` = **0.695128481**, `Var_Q(u)` = **4.046e-06**, hence
`sd(u)` = **2.011e-03**.
[CORRECTION - 2026-09-30]: the original text continued "— which matches the observed
seed-to-seed spread across three draws (0.695589, 0.695120, 0.691236)". **That
corroboration was spurious and matched the wrong number.** The diploid seeds were derived
as `seed + h`, so `release(seed=s)` haplotype 1 was BIT-IDENTICAL to `release(seed=s+1)`
haplotype 0 — measured, 0 of 305,886 positions differ — and the three "independent" draws
were four correlated haplotype draws. The overlap inflated the observed spread about 2.2x,
which landed it coincidentally on the HAPLOID `sd(u)` when the diploid prediction is
`sd(u)/sqrt(2)` = 1.42e-03. Two errors cancelling. Fixed: the streams now come from
`SeedSequence(seed).spawn(2)`. Ten disjoint-stream haploid draws give 2.199e-03 against the
exact 2.011e-03, which is a real corroboration. The utility the tilt can buy is
`E_Q[u] - E_R[u]` over [0, `eta_tau`], which is about `eta_tau * Var` =
**2.222e-06**. So the signal is 2.2e-06 against a draw-to-draw noise of 2.0e-03, a
ratio of about **1:900**. Resolving it would need n ~ (2.0e-03/2.2e-06)^2 ~
**8.3e+05 draws**, at 9.4 s each — about 90 days of compute for one target at one `tau`.

**AND THE REASON IS EXACT, NOT EMPIRICAL.** The total tilt pressure summed over
positions is `sum_v eta_tau * beta_v = eta_tau * 1 = eta_tau`, because `beta` sums to 1.
So the expected number of positions at which the tilt changes the draw is of order
`eta_tau` — **0.549 at `tau` = 0.5, 1.47 at `tau` = 0.9** — and it CANNOT grow with T.
Measured per-position tilt at mean `beta`: 3.28e-07 / 1.80e-06 / 4.81e-06 at
`tau` = 0.1 / 0.5 / 0.9. A mechanism that reweights 305,886 positions by one part in a
million, and whose entire budget summed over those positions is 0.549, moves the
released path at a fraction of one position.

**WHAT THIS IS NOT.** It is not a bug, and four things say so independently: `log_z_p`
tracks `eta_tau` at a constant 69.5%, which is the first-order identity
`log Z_p ~ eta_tau * E_R[u]` with `E_R[u]` = 0.694; `ffbs` provably draws from
`Q^tau_p` on the toy, where the per-position tilt is 1.28 at `tau` = 0.9 rather than
1.000005; the two independent `log_z_p` routes agree to 1e-12; and the seed-to-seed
spread is exactly the `sd(u)` computed from the second derivative.

[INCORRECT] - **WHAT THE LEVER IS, and it is not `tau`.** (The full original paragraph is preserved verbatim below the correction, because it runs to fourteen lines and the verifier requires the correction directly beneath this marker.)
[CORRECTION - 2026-09-30]: **TWO ERRORS, AND THE CONCLUSION INVERTS.** The quoted
`sum_v beta_v^2` was the wrong statistic, and concentrating `beta` is measurably NOT the
lever. Detail below, then the original text verbatim.

*The number was the wrong statistic.* `sum_v beta_v^2` is **4.556723e-06**, not 3.27e-06.
The quoted figure is `mean(beta_v)` = `1/T` = **3.269192e-06** — the value `sum beta^2`
takes only for UNIFORM weights. Two quantities that coincide in the uniform case were
conflated, and the coincidence is why it looked plausible.

*And concentrating `beta` is NOT the lever — measured directly.* `Var_Q(u)` computed from
`d2 log Z_p / d eta2` under a `beta_v` held uniform over the top-n positions by inverse support:

| `beta_v` scheme | sum of squared weights | `Var_Q(u)` | gain at `tau`=0.5 | `E_Q[u]` |
|---|---|---|---|---|
| inv_support (current) | 4.557e-06 | 4.046e-06 | 2.22e-06 | 0.6951 |
| uniform over all T | 3.269e-06 | 4.203e-06 | 2.31e-06 | 0.7075 |
| top 100,000 | 1.000e-05 | 1.140e-05 | 6.26e-06 | 0.6610 |
| top 10,000 | 1.000e-04 | 1.114e-05 | 6.12e-06 | 0.2926 |
| top 1,000 | 1.000e-03 | 1.113e-05 | 6.12e-06 | 0.0873 |
| top 100 | 1.000e-02 | 6.676e-05 | 3.67e-05 | 0.0494 |

Concentrating the budget **3,000-fold** raises `sum beta^2` by **2,200x** and `Var_Q(u)` by
only **16.5x** — leaving the gain **55x below** the per-draw `sd(u)` of 2.011e-03 — while
collapsing `E_Q[u]` from 0.695 to 0.049, because a utility defined on 100 positions is a
different and far weaker quantity. **The reweighting does not work, and cannot.** The reason
is that **94.8% of `Var_R(u)` is between-position COVARIANCE** from the mosaic's segment
structure, not the diagonal: `sum_v beta_v^2 p_v(1-p_v)` = 2.115e-07 against
`Var_Q(u)` = 4.046e-06. The segment count is set by the chromosome's genetic length and
`Ne*r`, which no choice of `beta` touches. The note's own closing caveat had this right and
contradicted its headline three paragraphs above.

**So R3 is NOT the answer to R5**, and R5's option list is wrong where it says concentrating
`beta` "raises `sum beta^2` directly" as a fix. R3 remains worth settling on its own merits
(call rate is not informativeness), but it will not make the mechanism measurable. What
remains is the structural tension: `phi_v` must decompose per position for `ffbs` to sample
`Q^tau_p` EXACTLY, and that requirement is what forces `u` to be a near-uniform average over
T positions, which is what makes `Var_R(u)` tiny, which is what makes the tilt powerless. A
utility with real variance would not factor, and an approximate draw carries no bound at all
— which is precisely why the k-mer Jaccard was ruled out as the in-sampler potential. **That
tension, not `beta_v`, is the thing to take to the PI.**

**THE ORIGINAL PARAGRAPH, PRESERVED VERBATIM** (the `[INCORRECT]` marker above points here;
`docs/memory.md` is append-only and the wrong version is what makes the reasoning legible):

> **WHAT THE LEVER IS, and it is not `tau`.** `docs/terms.md` `beta_v` says "if the
> achieved utility is too low, the lever is `tau` — a weaker privacy claim — not a
> different segmentation". Measured, `tau` is a very weak lever: `tau` = 0.9 flips about
> 1.5 positions of 305,886. The quantity that governs the achievable gain is
> `Var_R(u)` — how much the reachable utility varies across `R_D` — and that is set by
> how CONCENTRATED `beta` is, not by `tau`. `sum_v beta_v` is pinned at 1, but
> `sum_v beta_v^2` is not: it is **3.27e-06** under the near-uniform
> `w_v = 1/support(v)`, whose median support is 88 and which the 2026-09-28 note already
> described as "nearly uniform on this cohort". **This makes deferred decision R3 the
> highest-value open question in the project, not a cosmetic one.** R3 asks whether
> `beta_v` should emphasise rarity (`-log f`) instead of call rate; the answer now also
> determines whether the mechanism does anything measurable at all. ⚠ And it is a
> PRE-REGISTRATION item under constraint 5, so it must be settled before
> `A-LCK-preregister`, not after seeing an attack result.

⚠ **This does NOT contradict the 2026-09-28 note "the budget is invariant in T".** That
note is about the total tilt between the best- and worst-matching path, which is
`exp(eta_tau)` however the path is partitioned, and it is correct. What the note does not
say, and what a reader would wrongly infer, is that the ACHIEVABLE gain is that large:
the tilt can only move mass within the range of `u` that `R_D` actually reaches, and that
range is `E_R[u]` +/- 2.0e-03, not [0,1]. Whether `Var_R(u)` itself shrinks with T is NOT
established here — the mosaic's segment count is set by the chromosome's genetic length
rather than by T — and must not be claimed either way.

**Provenance for the three exact quantities.** Central differences of `log_z_p` at
`eta_tau(0.5) +/- 0.05`, three forward passes, no sampling; `log_z_p` itself is a
difference of two sums of magnitude 585 over 305,886 terms, so it carries about 1e-09 of
floating-point error against the quantities reported.

### Two things built that were not asked for, and why

- **`u_ceiling` is now a PRIMARY term.** `utility_retained` was defined against a
  ceiling of 1 ("the target's own path, u = 1 by construction"). On a real external
  target under `beta` scope `chain` the ceiling is **0.815455**, because the target has
  no call at 16.4% of the weight. Quoting 0.6940 against 1 understates the mechanism;
  quoting it against 0.8155 without saying so overstates it. Both now have to be
  reported, and `u_ceiling` is stamped on every release record.
- **A diploid release is ONE release.** `u_diploid = (u1 + u2)/2` factorises exactly —
  `exp(eta (u1+u2)/2) = exp((eta/2) u1) exp((eta/2) u2)` — so the two haplotype chains
  are independent at HALF the exponent and `log Z_p` is the sum of the two, landing in
  [0, eta/2] + [0, eta/2] = [0, `eta_tau`]. One private input, one `tau`, one bracket,
  which is what "one release per genome" means. Running each haplotype at the full
  `eta_tau` is the SUM aggregator in disguise; measured on the toy at `tau` = 0.5, the
  two `log_z_p` values then sum past `eta_tau`, so the test can tell them apart.


## 2026-09-29 (audit) — A-FIX-independence-audit: four auditors, and two criticals of my own

**Method.** Adversarial workflow `wf_814ac6b2-5f2`, four auditors on independence, cohort
correctness, target correctness and provenance, each required to verify by execution. One
auditor's transcript was briefly read as "died" and its lead chased by hand; all four in fact
returned.

**What held.** The cohort arrays are correct: every quantity re-derived independently with
bcftools/awk matched, and all 12 artifacts rebuild **byte-identically**. Sidecar coverage
84/84 with 0 stale digests, axis-digest sensitivity correct on all three axes, six of six
guards firing, and the dev fixture confirmed a clean subset of OUR axis (its digest
reproduces from our VCF; PanMixer's gives a different one).

### Two criticals in code I had written hours earlier

**1. The split-ALT join fabricated homozygous-reference calls.** It inferred "this haplotype
is reference" from the absence of a `1` at every panel record for our alts — valid only if
the panel covers ALL of them. Where it covered 2 of 4, a haplotype carrying alt 3 was written
as reference: **6,309 fabricated calls per target**, booked as reason 0 `called` under a
policy documented to mean "honest, unknown". Restricted to **fully covered** rows: **1,856**,
not 8,377 — matching the auditor's independent count exactly — where a `(0,0)` is conclusive.
Measured after the fix: 957 of the 1,856 are genuinely homozygous reference, 899 carry a
non-reference allele.

**2. The recommended `beta_v` normalisation rested on a false premise, and it was not the
recovery's fault.** `docs/plan.md` justified normalising over the panel-matched set with
"identical for all ten targets" — but that was attached to `in_chain & (target_reason == 0)`,
which is **NOT** invariant: **4 distinct sets across ten targets**, because a contradiction at
a poly-A repeat drops a row for one target and not another. Normalising over it is precisely
the constraint-1 violation the same section warns voids Theorem 1.
Fixed by emitting **`readable.npy`**: a row is readable if the panel has a usable record for
it, and both tests read only the panel's **site list**, the same file for every target.
**Verified 1 distinct set across all ten, locked by a test.** Normalising over
`in_chain & readable` gives max `u` = **1.000** over **253,846 of 305,886** chain positions
(83.0%), against **0.815** over all of them.

### Two more criticals, found independently by two auditors

- **`git_rev` stamped HEAD with no dirty check**, so **64 of 84** artifacts named a commit
  whose code provably could not emit the parameters they carried. Now `<sha>` when clean and
  `<sha>-dirty+<diff digest>` otherwise. Provenance that confidently names the wrong cause is
  worse than none.
- **`store.check_axis` was dead code** — defined, unit-tested, never called from `src/` — so a
  cohort directory mixing two site axes was consumed without error. Now invoked in both
  stages, and it immediately caught a shortcut in one of my own tests that wrote
  `in_chain.npy` with a bare `np.save`.

### Major

- **The `vg` digest pin was theatre.** No code read `tools.vg.path`, and the test hashed the
  *config's* path rather than the resolved binary — an impostor on `$RANPANMIXER_VG` ran while
  the test passed. `deconstruct` now resolves vg through the manifest, verifies the digest of
  the binary that will run, and records its path, version and digest. vg's version decides ALT
  ordering and therefore every allele index, so this is the one pin that matters.
- **`CLAUDE.md` still said PanMixer supplies the cohort data and to use their env.** That is
  the mechanism that re-seeds the coupling — this file records it happening before with
  constraint 1. Rewritten, and held at exactly 150 lines.
- **All ten target VCFs and the attack DB were cut from PanMixer's `starting_data`**, with
  that path in every header. Re-cut from `data/raw/`; verified no header names their tree.
- **The sidecar `sha256` was write-only.** Added `store.verify`; all artifacts pass.
- **The support/`INFO/AN` canary never enforced** (`strict=False` at its only call site), so
  real row drift would write `support.npy` and exit 0. Now strict, result recorded.
- **`deconstruct` and `drop_chm13` wrote no sidecars**, breaking the lineage from GFA to
  arrays. Both attach provenance now, including the vg version.
- **Stale PanMixer-derived counts corrected** in `terms.md`, `plan.md` and the manifest:
  n_sites **340,849** (not 340,824), T **305,886** (not 305,887), dropped_below **34,371**
  (not 34,345), **767,376** declared alleles, **41** sites at support<=1 (not 89).

### Process failures of mine worth recording

- I wrote a `vg` sha256 into the manifest with its **tail fabricated** from a truncated
  terminal display. Corrected, and a test now compares declared tool digests against the
  binaries on disk.
- I **committed while the docs contract was failing**, because a `&&` chain meant to gate the
  commit did not: piping pytest through `tail` masked its exit code. Fixed in the following
  commit; check `PIPESTATUS` or capture the exit code explicitly.

---


## 2026-09-29 (later) — A-IMP-preprocess: independence from PanMixer, and a real join bug

**Goal.** Dylan's direction: the tool must be COMPLETELY INDEPENDENT of PanMixer, which is
here only as a published example of a problem in the same space. Plus re-run chr21 through
the committed code path rather than by hand, and recheck everything.

### Independence — done and verified

- **Our own conda env** `ranpanmixer` (`environment.yml`, plus `requirements.txt` for a
  venv). Nothing from their `environment.yaml`. `src/` imports only stdlib, numpy and
  pyyaml; scipy is present for the sampler. **pytest is in our env**, so both suites run
  normally (it is absent from `panmixer`, which is why `test_preprocess.py` is self-running).
- **The dev fixture is rebuilt from OUR VCF.** It had been sliced from PanMixer's.
- **`vg` relocated** to `data/tools/vg`, pinned by sha256 in the manifest, overridable via
  `$RANPANMIXER_VG`. Verified by re-running chrY deconstruct through it: 34,014 records.
- **The 1000G panel staged** to `data/raw/`, recorded as evaluation-only.
- `grep` finds **no** reference to `external/PanMixer` in `src/`, `configs/` or `tests/`, and
  tracing every provenance sidecar shows PanMixer's chr21 VCF sha256 `5b880a3a…` appears in
  **no** artifact's inputs.

### chr21 re-run through the committed path

8:33 / 9.475 GB, **340,849 records, 45 samples** — reproducing the by-hand run exactly
(8:22 / 9.477 GB / same counts). The CLI performed its own decompression, confirming the gz
fix. So the artifacts are now provably the output of committed code.

### ⚠ A REAL BUG: the strict join could never match a multi-allelic record

**The 1000G panel is entirely decomposed into biallelic records** — every one of its
1,002,752 chr21 records declares exactly one ALT — while **7.9% of ours (27,021) are
multi-allelic**. A strict (POS, REF, full-ALT-string) join therefore **cannot succeed on any
of them**, for a purely notational reason, and the first implementation silently discarded
target information we actually held.

Recoverable by matching each of OUR alts against its own panel record and mapping back to
OUR allele numbering:

| | records |
|---|---|
| strict match | 258,610 |
| **recoverable via split alts** | **8,377** (all multi-allelic; 8,283 inside `chain_span`) |
| genuinely absent from the panel | 73,862 |

Implemented as `split_alt_join`, on by default. Effects, measured across all ten targets:
sites matched **258,610 -> 266,987**; chain positions with a target call **82.4% -> 85.1%**;
and once the accounting was made consistent, representable non-ref calls **80.0% -> 83.5%**
(a mean of 2,269 calls per target were being miscounted as unrepresentable).

⚠ **The decomposition assumption is not universally safe, and the counter proves it.**
**3 contradictions** across the ten targets (0.04% of recoveries) — a haplotype reading as
carrying two different alts. Examined: chr21 **28,186,675**, REF `C`, our alts
`CA,CAAA,CAA,CAAAA` — a **poly-A repeat** where the panel's biallelic records are NOT
mutually exclusive, and haplotype 0 is `1` at both `C->CA` and `C->CAA`. The code counts it
and **leaves the row unresolved rather than guessing**, and those calls fall back to
`pos_only`. That is the correct behaviour and the contradiction counter is what makes it
visible.

**Second-order bug this exposed:** the recovery initially improved `path.npy` while
`representability.json` still classified those calls as unrepresentable — the two artifacts
disagreed about the same fact. Fixed by deferring the classification until the reconstruction
resolves.

### Three other items closed

- **`support(v) == 0` is excluded from the chain**, not patched. `w_v = 1/(support+1)` was
  the tempting fix and is wrong: at support 0 no donor has an allele, so the sampler has
  nothing to emit and RENORMALISE has an empty pool — the position is unusable, and a finite
  weight it can never earn would just hide that. Removes 0 positions today (the one such
  site is already outside the span), so the guard no longer rests on that luck.
- **`target_reason.npy`** distinguishes `no_record_in_target_callset` from
  `record_present_but_no_call`. Measured: the inert 19% is **53,861 code 1 and ZERO code 2**
  — a phased panel has no missing genotypes inside its own site list. So the unmatched set
  is identical for all ten targets and provably carries no per-person variation, which is
  the property needed to normalise `beta_v` over the panel-matched positions legitimately.
- **A fabricated checksum, caught by the check written to catch it.** A `vg` sha256 was
  written into the manifest with its tail invented from a truncated terminal display. It
  looked plausible. Corrected, and a test now compares every declared tool digest against
  the binary on disk — verified to fail against a corrupted one. A checksum that was never
  computed is the worst violation of the traceability rule, since being checked is its only
  purpose.

**Provenance.** 40 tests pass in the `ranpanmixer` env. Audit workflow `wf_814ac6b2-5f2`
launched four auditors; all four were cut off before returning, but one flagged the
multi-allelic lead that this entry's main finding came from — chased and confirmed by hand.

---


## 2026-09-29 — A-IMP-preprocess: the preprocessing pipeline, built and run end to end

**Goal.** Implement preprocessing for the cohort and one target, verify every step, and
audit the result. First code this project has that is its own rather than PanMixer's.

**Result: it runs, and every artifact cross-checks.** 17 unit tests, 28 consistency checks
on the real output, all passing. Still no mechanism and no RanPanMixer number.

### Our own chr21 build, and how it differs from PanMixer's

`vg deconstruct` on the chr21 PGGB GFA, **8:22 wall / 9.5 GB peak at -t 32** — far cheaper
than the 1-4 h estimated. ⚠ **Sample grouping was the risk and it is settled**: paths are
PanSN `SAMPLE#HAP#CONTIG` with each assembly split across many contigs, and `-H` is
deprecated in vg 1.68, but a chrY smoke test (2:01, 34,014 records) produced 17 clean sample
columns. chr21 produced **45** (44 donors + chm13), as it should.

**Our VCF is NOT the published one, and every derived number must say which it is.**

| quantity | OURS (vg 1.68) | PanMixer's (vg 1.36) |
|---|---|---|
| chr21 records | **340,849** | 340,824 (+25) |
| sample columns before chm13 | 45 | 45 |
| missing matrix cells | **1,189,258** (3.965%) | 1,278,653 (4.26%) |
| bare `.` GT fields | **0** | 237,594 |
| declared alleles | **767,376** | 767,324 |
| T inside `chain_span` | **305,886** | 305,887 |
| dropped below span / above | **34,371** / 592 | 34,345 / 592 |
| sites at support 0 | 1 | 1 |
| max alleles at one record | 90 | 90 |

**⚠ AND THE ALLELE INDICES DIFFER — the most important finding of the run.** "340,849 vs
340,824" understates it badly. On the strict (POS, REF, ALT) key only **94.9%** of records
match: 17,454 are ours-only and 17,429 theirs-only. But comparing on (POS, REF) and then
asking whether the ALT *set* matches shows the cause:

| | records |
|---|---|
| same ALT set, order may differ | **339,358** |
| genuinely different ALT set | 386 |
| (POS, REF) only in ours | 1,105 |
| (POS, REF) only in theirs | 1,080 |

So the 5.1% mismatch is **almost entirely ALT REORDERING**: `T -> TAA,TA,TAAA` in ours
against `T -> TAA,TAAA,TA` in theirs, the same variant with the alleles listed in a
different order. **Allele index 1 in our VCF is not allele index 1 in theirs** at those
sites. Any PanMixer-derived artifact joined positionally to ours — the allele matrix, a
released `new_haplotypes`, the attack database masks — would be silently wrong there, with
no error and no shape mismatch. This is precisely the hazard `site_axis_digest` exists to
catch, and it is now a measured fact rather than a worry: the two VCFs' axes genuinely are
different objects. **Never mix artifacts across the two builds.**

⚠ **vg 1.68 calls 89,395 more genotypes than 1.36 on the same graph** — and emits **zero**
bare-`.` GT fields where 1.36 emitted 237,594, writing `.|.` instead. Both counts were
verified independently with bcftools against both VCFs. Two consequences: the "haploid
genotype" cause of a `-1` does not merely happen to be always `.` on our data, it **does not
occur at all**; and the chrY divergence measured on 2026-09-28 (92.7% of tuples matching)
was atypical — chr21 differs by 25 records in 340,824, i.e. 0.007%.

### Verified, not assumed

- **`INFO/AN` == the non-missing count, exactly, at all 340,849 sites.** The prior
  measurement was on PanMixer's VCF; it holds on ours too. That makes the cross-check a
  free row-alignment canary, and it is not vacuous — the tiny fixture carries one
  deliberately wrong `AN` and the check fires on exactly that record.
- **Reason codes hand-verified.** Every count on the fixture matches hand calculation,
  including the precedence collision at pos 720 where S3 is both conflict-named and
  would-be-inherited: conflict wins, as pinned.
- **Code 4 is not "unexplained".** 571,182 uncategorised cells, of which **99.3% sit in a
  run of >= 10 consecutive missing** — they are assembly gaps, which deliberately get no
  per-cell code and live in `missing_runs.tsv` as intervals instead. The bucket is behaving
  exactly as designed.
- **The coordinate cut is doing its job.** 8,327 records sit in [10,326,676, 12,968,320] —
  inside the genetic map but outside `chain_span` — and every one carries cM **0.584144
  exactly**. That is the centromeric plateau the cut exists to exclude, observed directly.
- **Representability on all 10 targets** (in-region, so a partial fixture cannot report its
  own window as a failure): mean **80.0% exact / 10.2% invisible / 21.5% of information
  invisible**, against the 2026-09-27 figures of 80.4% / 9.9% / 20.5% measured on six
  samples and the published VCF. Independent reproduction.

**⚠ NEW FINDING — representability is ancestry-dependent.** The two AFR targets are the
worst on both axes: NA19350 and NA20294 carry the most non-reference calls (78,595 and
79,060 against ~60,000 for EUR/EAS) and have the lowest exact-match rate (**78.9%** against
81.1% for HG01468, AMR). The cohort represents a EUR target better than an AFR one. This is
exactly what the deliberately multi-continental target set was chosen to expose, and it is
a utility-ceiling effect that a single-population test set would have hidden entirely.

### Two real bugs, both found by running on real data

1. **`_csv.Error: field larger than field limit`** — `load_sites` could not read its own
   output. The hypervariable record's ALT string is **23.2 MB** and Python's csv module
   refuses fields over 128 KB. Fixed with `csv.field_size_limit`. Found in seconds because
   the dev slice was chosen to contain that record.
2. **`sites.tsv` had no provenance sidecar** — the one artifact that defines what every row
   means was the only one written outside the `store.save()` chokepoint, because
   `write_sites_tsv` streams it to compute the digest in the same pass. Added
   `store.attach_prov()` as the narrow, deliberate exception.

Both now have regression tests. A third apparent bug was **my own audit check** being
wrong, not the code: "all cM outside the span are NaN" is false by design, because the
centromeric plateau is inside the map and outside the span.

**One fixture error of my own**, caught by the canary it was meant to test: record 522 had
`AN=4` where the truth is 5. Corrected, so only the deliberate error at 1300 remains.

### Provenance

`data/cohort/chr21/` 240 MB, 12 artifacts each with a sidecar, one axis digest
`251b554f...` across all of them and all ten target paths. `data/targets/<id>/chr21/` holds
`path.npy` + `representability.json` per target. Dev fixture at `data/dev/` is the identical
code path via `--root`. `pyyaml` was added to the `panmixer` conda env (absent upstream).

---


## 2026-09-28 — A-FIX-panmixer-framing + A-DAT-genetic-map-parm + A-DAT-input-set + A-DAT-evaluation-fit

**Goal.** Start the build. Create a clean branch, decide the chain architecture, settle the
preprocessing input and output sets, and correct the documentation where the PanMixer
comparison had been forced onto a problem it does not describe.

**Method.** Two adversarially-verified workflows (10 and 12 agents) plus hand verification.
Every number below is measured in this project on chr21 unless labelled otherwise.

### Decisions (Dylan)

1. **Model C is the chain architecture** — one chain over the chromosome, units = individual
   VARIANTS. Model B (units = LD blocks) is the retained fallback; **model A is ruled out**.
   This answers open question 8, raised 2026-09-22 and deferred to the PI.
2. **The scope decision is narrowed.** *"Maximise reuse of PanMixer, change exactly ONE
   thing"* imported PanMixer's threat model — target-as-cohort-member — into a design whose
   target is EXTERNAL to everything. Dylan: *"Boundaries aren't a weak function of the
   target's genotypes because the target is completely outside of any other dataset."*
3. **Write fresh, do not port.** Copy functions from PanMixer where they are right, but the
   file structure and pipeline are ours, so nothing is inherited without a decision.
4. **The pipeline starts from our own inputs and produces all outputs.** `vg deconstruct` is
   pipeline step 1. Accepted knowingly: our VCF will not match PanMixer's.
5. **`missing_policy` = RENORMALISE.** The cost — penalising poorly-assembled donors for a
   defect that is not theirs — was accepted deliberately, to be noted in the writeup rather
   than modelled around at this stage.
6. **Target mapping is a separate stage and a separate command**, so a new target does not
   re-run the cohort stage.
7. **`chain_span` for chr21 starts at 12,968,320**, flagged for revisit.
8. **The k-mer Jaccard is deferred** to a post-build evaluation/reporting tool.

### The framing audit

91 findings filed across `CLAUDE.md`, the five docs and the three primer sections;
adversarial verification **upheld 76, rejected 15, and found 51 further instances** the
auditors missed. The rejections were almost all the predicted failure mode — a correct
description of PanMixer read as a claim about us.

Highest-leverage: **`CLAUDE.md` constraint 1** ended with an unscoped *"Recompute
leave-one-out"* inside the hard rules, in the file auto-loaded every session — so the error
regenerated every session. Also `terms.md`'s `loo_cohort`, which defined OUR cohort as
`HPRC \ {target}` at Status PRIMARY; `plan.md`'s exit gate 5, which made PanMixer benchmark
compatibility a blocking gate of the sampler; and `03-panmixer.md:2363`, which asserted
block boundaries are a function of the target's genotypes — flatly false of our setting and
internally contradicted by line 1220 of the same file.

### The genetic map, and the p arm

**The defect found is worse than the one we were looking for.** `np.interp` clamps by
default, so all 26,018 variants below the map's start (10,326,676) take its first value —
and the map's own first interval (10,326,676 -> 12,968,320, spanning the centromere) carries
**dcM = 0.000000**. Together: **34,345 consecutive positions with identical cM**, so
`Delta_x = 0`, `P(switch) = 0`, and the sampler emits **one donor's real haplotype verbatim
across 7.25 Mb** — departure (b) at ~200x the scale of the worst single block. A privacy
failure, not a utility one, and silent.

**No map fixes it.** deCODE 2019 (native hg38, highest resolution) starts at **14,143,522**,
worse — it would leave 47,805 variants (14.03%) unmapped against 26,018 (7.63%). The 1000G
hg38 map has a single value-0 interval over chr21:0–10,326,676. Bherer 2017 is GRCh37 and
starts later still. pyrho hg38 (Spence & Song 2019) DOES cover from **5,088,754 bp** and
would give every variant a distinct cM — but it spends **9.98 cM (ACB) / 9.15 (CEU) across a
single 3.4 Mb interval defined by only two markers**, attributing 23–26% of chr21's genetic
length to a region with zero observed crossovers. Interpolation across a void.

**Why the region is like that — the p arm.** Every chromosome has a centromere dividing it
into a short arm (**p**, French *petit*) and a long arm (**q**). Five human chromosomes —
13, 14, 15, 21, 22 — are **acrocentric**: the centromere sits near one end, so the p arm is
small and is almost entirely repetitive DNA plus tandem arrays of **rDNA** genes. Those
arrays are near-identical BETWEEN the five, so assemblers cannot place them uniquely and
GRCh38 fills them with N. chr21's p arm is 0–~12 Mb.

Five independent legs support excluding it: GRCh38 chr21 5,010,000–10,814,560 is placeholder
model sequence (25 contigs, fabricated N-gaps, 22 of exactly 50 kb); 1000G's strict
accessibility mask rules **26,014 of 26,018** (99.98%) inaccessible; our cohort matrix is
**36.76%** missing there against 6.57% elsewhere; deCODE's pedigree map assigns **0.0 cM/Mb
to every 1 Mb bin from 0–13 Mb**; and a four-generation CEPH pedigree recorded *"Not a
single allelic recombination was observed on the p-arm"* across 107 transmissions.

⚠ **Two smaller zero-distance defects remain and the cut does not touch them**: 592 variants
above the map end give a 593-position zero-distance run at the telomere (UNDECIDED), and 618
interior map intervals with dcM = 0 hold ~1,650 kept variants (correct to leave — real map
plateaus). ⚠ **What the cut costs is unmeasured** — the share of a target's `-log f`
discarded has not been computed.

### The input set

- **PanGenie drops out entirely.** Exactly two live consumers, both deleted by our design;
  no evaluator reads it. Saves a 5.06 GB download.
- **The external panel drops out of the mechanism** under model C — no blocks, no PLINK.
  It survives only as the attack database.
- **`INFO/AN` is bit-identical to `support(v)`.** `np.array_equal` is True across all
  340,824 chr21 records. Support needs no computation, and the comparison doubles as a free
  row-alignment check.
- **Do NOT compute LD blocks from the 44-sample cohort** (for the B fallback): 7,858 blocks
  covering 50.4% of variants against 14,137 / 74.2% from the panel; largest block 2,562
  against 618. Only 44% of SNPs clear MAF 0.05 in 88 haplotypes.
- **The haploid-genotype cause of `-1` does not exist on chr21.** All 237,594 non-piped GT
  fields have exactly one distinct value, `.`. Four causes, not five.
- **The graph is published only as GFA** for PGGB freeze1 — no GBZ, no snarls. chr21 =
  434,772,428 B; whole genome 15.64 GB compressed / 86.5 GB uncompressed.
- **`vg` 1.68 does not reproduce `vg` 1.36.** On chrY it recovers 99.55% of the published
  VCF's sites by POS but only **92.7% by (POS, REF, ALT)**, concentrated at multi-allelic
  sites. Also: the published VCF carries 45 samples, our local chr21 carries 44 — `chm13`
  removed by PanMixer's own step — so a rebuild does not drop in as a replacement.
- **The graph buys little the VCF lacks**: 98.37% of consecutive top-level snarls share a
  boundary node, so donor switching is provably a valid graph walk from `AT` alone; all
  340,824 alleles are explicit sequences. The one irreplaceable use is mapping a target's
  reads — and HPRC marks short-read mapping **"untested"** for PGGB.
- **The graph does carry the p arm**: the CHM13-backbone VCF has **185,414 records** in
  `chm13#chr21:1-10 Mb` against 22,746 for GRCh38. But that recovers sequence, not genetic
  distance, and missingness there is 43.9%.

### The evaluation

**Nothing in PanMixer's suite tests the claim we make.** Read every evaluator:
- `af_loss` and `ld_loss` reconstruct an edited cohort (`new_pangenome[subject_index] =
  released_array`). We never edit the cohort, so they are not defined for our release
  object. ⚠ A `--dont_replace` flag exists that skips the substitution, so this is
  conditional rather than structural — but the decision stands regardless: Dylan,
  *"we don't use af_loss for our project at all."*
- **The gap score is the one that transfers**, after three changes: "self" becomes the
  external target's true path; the attack DB excludes the target; the AF source is disjoint
  from the cohort. ⚠ Not reusable verbatim — `score_genotypes_remove_missing` writes
  `test_hist.png` to the working directory on every call, 11x per row.
- **Their shape assertion is vacuous.** Both sides equal `site_mask.sum()` by construction,
  so a differently-ordered release is silently mis-scored rather than erroring. This is the
  row-order hazard, located.
- **Two of their metrics are broken**: `ld_loss`'s window is **5 bp, not 5 kb** (7,629 pairs
  on chr21 instead of ~10.6M), and the haplotype gap score sets `ref_af = alt_af`, so
  reference matches are weighted by the alt frequency.
- The test we need is **two-target indistinguishability** — new entry in `docs/terms.md`.

### Two implementation facts worth not rediscovering

- **The budget is invariant in T.** `u = sum_v beta_v phi_v` with `sum beta = 1`, so the
  total tilt between a perfectly-matching and a non-matching path is `exp(eta_tau)` however
  the path is partitioned. At tau = 0.5: per-position tilt 1.0000055 at T = 100,757 against
  1.0000016 at T = 340,824, total 1.7321 either way. Finer T does NOT spread the budget
  thinner.
  [INCORRECT] - ⚠ What it DOES cost is float32 headroom: the per-position tilt under C is only **13.5x** the float32 resolution against 45.7x under blocks, so float64 stored alphas are load-bearing. And the thin-budget problem is a property of **tau**, not of T.
  [CORRECTION - 2026-10-10]: the 13.5x is computed over the 340,824-record SITE AXIS, not
  the chain. The chain is the denominator — `eta_tau * mean(beta_v) = eta_tau / T` — so at
  T = 305,886 it is 0.549306/305,886 = 1.80e-06 and **15.1x**, not 13.5x. The 45.7x for a
  block partition (T = 100,757) is correct. ⚠ **The conclusion is now stronger than this
  bullet claims, not weaker:** under decision B (2026-10-08) `beta_v` sums to 1 genome-wide,
  which drops the per-position tilt BELOW float32 eps, so float64 is not a 15x margin but
  the only thing preventing a silently untilted draw. Recorded in `docs/bugs.md`.
- [INCORRECT] - **Never materialise `A`.** Li-Stephens `A = (1-s_t)*I + (s_t/K)*J` is rank-one plus diagonal, so the forward step is O(K) from one scalar sum and the backward draw likewise. A materialised `A` per position would be 340,824 x 88 x 88 x 8 B = **~21 TB**.
  [CORRECTION - 2026-10-10]: the arithmetic in that last sentence is wrong by 1000x.
  340,824 x 88 x 88 x 8 B = 21,114,728,448 B = **~21.1 GB**, not 21 TB. The conclusion
  (never materialise `A`) is unaffected, which is exactly why the slip survived six weeks —
  a number supporting a correct conclusion is never checked. The same error stood in
  `src/mechanism/ffbs.py` and `src/mechanism/hmm.py` and was fixed there 2026-10-09.

**Provenance.** Workflows `wf_1dd6bfa5-5a8` (inputs, 10 agents) and `wf_d0195c3c-52d`
(framing audit, 12 agents), both adversarially verified; `wf_785b6f66-755` applied the primer
corrections. Doc corrections in commit `1f8ea3f`. No mechanism has been run and no
RanPanMixer number exists.

---


## 2026-09-27 — A-THY-chain-architecture: three models, measured, and a toy that shows the tilt vanishing

**Why.** Dylan asked whether "the HMM is per block and the chain never crosses a block
boundary" is the right framing for OUR mechanism, explicitly from first principles rather
than from the docs. It is not: that sentence is a FINDING ABOUT PANMIXER'S CODE which
leaked into our spec as though it were a design principle. Nothing in the proposal or in
Theorem 1 requires it.

**First-principles position.** `baseline_model` is a model of what a plausible haplotype
looks like: a mosaic whose segment boundaries fall where recombination happened, at a rate
set by genetic position in centiMorgans. **LD blocks appear nowhere in that.** They are a
statistical summary of where recombination has historically been rare — downstream of the
same process the transition model already encodes. Using blocks INSTEAD of the map
substitutes a coarse proxy for the real quantity; using both double-counts.

**The three options, measured on chr21** (acquired GRCh38 map, Ne = 10,000, r = 1.26,
n = 86; analytic segment n/(Ne*r) = 0.0068 cM ~ 4.5 kb):

| model | chain | E[switches] | mean segment | alphas |
|---|---|---|---|---|
| A: per-block chains, no inter-block (pinned) | 100,757 separate | 104,379 | 0.39 kb | 69 MB |
| B: one chain, units = blocks | 100,757 | 3,532 | 11.60 kb | 69 MB |
| C: one chain, units = variants | 340,824 | 7,155 | 5.73 kb | 234 MB |

- Columns: one map, one formula, our real positions.
- Rows: A is the pinned spec; B and C are candidates.
- Synthesis: **A is 11x too scrambled and is the worst of the three.** Its fresh-uniform
  reset per block is algebraically the `d` -> infinity limit of the switch formula. Measured
  inter-block switch mass is median **0.0022** against the **0.9884** independence implies,
  and 98.5% of block pairs are below half of it. B discards **51%** of the recombination
  (the within-block half); C is nearest the analytic segment length.

**Supporting measurements.** 89,087 of 100,757 block-dictionary entries (88.4%) are orphan
variants in NO called plink block, so block boundaries are mostly bookkeeping, not
recombination sites. Map resolution is 407 bp median against 60 bp variant spacing, so the
map already resolves finer than blocks do. The largest block is 618 variants over 39.5 kb
which C recombines through 73% of the time and B emits as one donor verbatim.

**Numerical findings.** An unnormalised forward pass underflows float64 after about **183
positions**; per-column normalisation (`log_z_p` = sum of log `c_t`) removes it entirely and
is preferred over log-space (cheaper, and `alpha_t` is already the distribution backward
sampling needs). **Storage must be float64**: float32 resolution 1.2e-07 against a
per-position tilt of 5.5e-06 is only 46x of headroom.

**The toy (4 donors, 7 variants, 3 blocks), model B.** Forward pass reproduces `Z_p` exactly
against enumeration over all 64 paths (1.330113 both ways); `log_z_p` = 0.285264 inside
[0, 0.5493]; E[u] moves 0.5000 -> 0.5386 at tau = 0.5.

**The tilt vanishing, shown directly.** Tiling the same 3-block pattern to longer chains,
holding tau = 0.5:

| T blocks | max psi | E[u] gain |
|---|---|---|
| 3 | 1.2654 | +3.86e-02 |
| 300 | 1.0024 | +3.72e-04 |
| 30,000 | 1.0000235 | +3.73e-06 |
| 100,755 | 1.0000070 | +1.11e-06 |

Exactly 1/T, and **identical under all three architectures** — the chain choice does not
touch the budget. The same toy at T = 3 gives E[u] gains of +0.0402 (A), +0.0386 (B),
+0.0371 (C): A buys marginally the most utility precisely BY discarding the linkage
constraint, which is also why its output is least haplotype-like.

**Status: OPEN, deferred to Dylan's PI**, written up as a decision in `docs/plan.md`. If C
is chosen the factorisation argument becomes unnecessary and the singleton special case
disappears, making C simpler than what is pinned rather than harder.

## 2026-09-27 — A-DAT-genetic-map: map acquired, and exit gate 3 is mis-specified

**Why.** The GRCh38 genetic map was the last unacquired input for `A-IMP-cohort-hmm`.
Earlier this session I reported the SHAPEIT4 b38 maps as returning 302 — that was wrong,
a redirect I failed to follow. Both candidates return 200 and were downloaded.

**Result 1 — the two candidate maps are the same map.** SHAPEIT4 and Beagle both yield
chr21 with **44,618 positions, identical positions, identical cM, maximum absolute
difference 0**. The choice is format and provenance only, not content. Two independent
distributions agreeing exactly is evidence the map is not a single-source artifact. Both
registered in `docs/data.md`; the SHAPEIT4 copy is the working file, being pinnable by
commit SHA rather than served from a lab webserver.

[INCORRECT] - **Result 2 — 7.63% of our chain positions sit outside the map.** Map span is 10,326,676 - 46,680,243 bp; our chr21 variants span 5,713,651 - 46,699,788. So **26,018 variants (7.63%) fall below the map start** and 592 (0.17%) above.
[CORRECTION - 2026-10-10]: re-measured directly from `data/cohort/chr21/positions.npy`
against the extracted map, the count below the map start on our 340,849-record build is
**26,044**, not 26,018 — a difference of 26 records (0.1%), which does not move the 7.63%
or any conclusion in this entry. The 592 above the map end and the variant span are both
confirmed exactly. **26,018 is PanMixer's vg-1.36 340,824-record figure**, and its presence
here is why two different totals for the cut circulated: 26,018 + 8,327 = 34,345 against
our 26,044 + 8,327 = **34,371**, the latter being what the cohort build actually drops.
The full measured split now lives in the `chain_span` and `map_coverage` entries of
`docs/terms.md`. The in-span portion — 8,327 positions on the map's own zero-cM plateau,
all at cM 0.584144 — is the part decision C does NOT reach. That region is the
acrocentric p-arm and heterochromatin, where no recombination map exists for a real
reason. Interpolation clamps, giving flat cM, so **296 multi-variant blocks there have
P(mosaic) = 0 exactly** and can only ever be one donor copied verbatim. This is the hazard
flagged from Pickrell's README, confirmed on our data; it is a property of the region, not
of the map we picked.

**Result 3 — the transitions work, and exit gate 3 as written does not test them.**
Map interpolated onto all 340,824 chr21 positions, published formula (Ne = 10,000,
r = 1.26, d = `Delta_x` * Ne * r, total switch mass = (n-1)(1-exp(-d/n))/n, n = 86):

| quantity | value |
|---|---|
| median `Delta_x` between adjacent chain positions | 1.811e-05 cM (about 12 bp) |
| per-step switch mass, median | 0.0026 |
| per-step switch mass, p75 / p90 | 0.0098 / 0.0304 |
| steps with switch mass in [0.1, 1] | 3.0% |
| expected donor switches across chr21 | **3,622.6** |
| P(block is a mosaic), median over 11,670 multi-variant blocks | 0.0866 |
| P(block is a mosaic), median over the 3,523 blocks with >= 20 variants | 0.2747 |
| PanMixer shipped, same measurement | 5.8e-12; 300/300 blocks single-donor |

- Columns: one map, one formula, our real chain positions.
- Rows: per-step figures are between ADJACENT variants; block figures span a whole block.
- Synthesis: median per-step switch mass is **4.5e8 times** PanMixer's, and the model gives
  mosaic segments of a few kb, consistent with the analytic n/(Ne*r) = 86/12,600 =
  0.0068 cM. **But exit gate 3 asks for "switch mass of order 0.1 to 1 at realistic
  separations", and that target is biologically wrong at this density.** Adjacent chr21
  variants are about 12 bp apart; switch mass near 0.1 there would mean roughly one
  ancestral recombination per 12 bp. The measured 0.0026 is the CORRECT value.

**Exit gate 3's second clause is mis-specified too, and more subtly.** It requires that
"a sampled block is demonstrably NOT a single donor copied verbatim". With the correct map
the median multi-variant block still has only **8.7%** chance of containing a switch, so
about 91% of blocks WILL be one donor copied verbatim. That is NOT PanMixer's failure
repeated: theirs was 100% from an arithmetic bug at 5.8e-12; ours is ~91% because **an LD
block is by definition a stretch with little recombination**, and v1's chain scope is one
block. The mosaic emerges ACROSS blocks — precisely the scope v1 defers. Inside v1 the tilt
acts on donor SELECTION per block, not on within-block recombination.

**OPEN — needs restating before implementation (Dylan).** Exit gate 3 should be rewritten
against quantities meaningful at block scope. Candidates: expected donor switches per
chromosome (3,622.6); median P(mosaic) among blocks with >= 20 variants (0.2747); or
agreement between measured segment length and the analytic n/(Ne*r). Also open: whether to
extrapolate cM below 10.33 Mb or accept 296 permanently single-donor blocks.

## 2026-09-27 — A-DAT-representability-split: 9.9% invisible, carrying 20.5% of the information

**Why.** The earlier entry today put the invisible share at 19.6%. Its median allele
frequency came out at 0.28 — implausible, since a variant that common should appear in 88
haplotypes. Re-measured by splitting on whether the graph has ANY record at the position.

**Method.** Same three 1000G samples outside the HPRC 44 (HG00096, HG00097, HG00099), chr21.
Three categories per target non-reference call: exact (POS, REF, ALT) match into the graph;
POS present but no exact match; POS absent entirely. Allele frequency is leave-one-out over
the 3,202-sample panel with the target's own two haplotypes removed, singleton-floored at
0.5 counts. Information is `sum of -log f`, which is the published linkage score of Eq. 5.

| category | share of sites | median AF | MAF < 1% | share of information | ratio |
|---|---|---|---|---|---|
| exact match — measured | 80.2% | 0.47 | 0.3% | 68.3% | 0.85x |
| POS only, alt differs — measured | 9.9% | 0.35 | 1.1% | 11.1% | 1.13x |
| **no record — INVISIBLE** | **9.9%** | **0.20** | **16.9%** | **20.5%** | **2.07x** |

- Columns: pooled over three targets; per-sample spread was under 2 points in every cell.
- Rows: only the third is invisible. The second HAS a chain position, because a record
  exists at that POS — the mismatch is a REF/ALT spelling difference, typically a
  multi-allelic graph context, and `phi_t` scores it as an ordinary mismatch.
- Synthesis: **one tenth of a real target's variants have no record at all, and they carry
  one fifth of its identifying information** — a 2.07x over-representation, with rare
  variants enriched about 56-fold (16.9% under MAF 1%, against 0.3% among exact matches).

**The framing that matters (Dylan).** This loss is NOT ON THE `tau` DIAL. It is fixed by the
graph's coverage and happens before `psi_t` is evaluated, so no `tau` recovers it: at
`tau` -> 1, with no privacy guarantee at all, the release still cannot express those alleles.
**The utility ceiling is therefore about 79.5% of the target's information, not 100%**, and
every target_fidelity number must be read against it. The same fact read from the privacy
side: 20.5% of the linkage score is annihilated unconditionally, for every target, for free.
Bad for a clinician, good for a privacy claim, and controlled by neither of us.

**Caveat for reporting.** The middle category is measured but unfair to the mechanism: the
target's exact allele is unavailable for a join reason, so `phi_t` records a mismatch the
sampler could not have avoided. Separate it when reporting utility, or the mechanism is
charged for a representation artifact.

## 2026-09-27 — A-DAT-representability: a fifth of a real target never reaches the mechanism

**Raised by Dylan.** A target variant absent from both the graph and every external callset
never enters any VCF. The release is a path through `G`, so that variant cannot be emitted —
but it also has no record, hence no chain position, hence no `phi_t` term and no `beta_t`
weight. It is removed AND unmeasured. Is anything accounting for it?

**Answer: no, and the magnitude is large.** Measured on six 1000 Genomes samples that are
not among the HPRC 44, chr21, using the strict (POS, REF, ALT) mapping:

| metric | per target |
|---|---|
| non-reference calls | 60,561 |
| with an exact record in the graph | 48,678 (80.4%) |
[INCORRECT] - | **invisible: no record, no chain position** | **11,883 (19.6%)** |
[CORRECTION - 2026-09-27]: 19.6% conflates two categories. Only about half of it — **9.9%** —
has no record at all. The other **9.9%** has a record at that POS with a different REF/ALT
spelling, so a chain position DOES exist and `phi_t` does measure it. See the
A-DAT-representability-split entry above. The 80.4% exact-match figure is unchanged.

Per-sample spread was tight: HG00096 80.3%, HG00097 80.1%, HG00099 80.3%, HG00100 80.6%,
HG00101 80.5%, HG00102 80.5%.

- Columns: one typical external target.
- Rows: "invisible" means no exact-match record exists, so no chain position is created.
- Synthesis (⚠ the invisible share is corrected to 9.9% above): **both bounds move the wrong way.** 80.4% is an UPPER bound on retainable
  fidelity, because having a record is necessary but not sufficient — the allele must also
  be carried by a donor, and 8,823 declared alleles have zero cohort support (2026-09-27
  entry above). And the invisible share — 9.9% as corrected, not the 19.6% written here — is
  itself a LOWER bound, because these are 1000G samples
  whose variants are by construction already known; a genuinely novel variant is invisible
  to this measurement too. The measurement is also SNV-biased, since 1000G short-read calls
  under-ascertain exactly the structural variation the graph is richest in.

**The consequence for our utility claim.** `u(p, y) = sum beta_t phi_t` is defined over
chain positions, and chain positions come from VCF rows. So **u can equal 1.0 while the
release retains only about 80% of the target's variation.** Every downstream instrument
inherits the blind spot: AF loss, LD loss and read mapping all operate on the graph's site
set and none can see a variant with no row.

**Not inherited from PanMixer.** Their target is a cohort member, so their haplotypes ARE
the VCF rows by construction and nothing they carry is inexpressible. This gap belongs to
the external-target setting and has no upstream precedent.

**Privacy is untouched.** The invisible variants cannot be emitted, so they are maximally
protected; and the mechanism's behaviour does not branch on whether the target had one, so
nothing leaks. This is a measurement and utility defect, not a privacy defect.

**Fix, and the shape it must take.** Dylan proposed adding a record at those locations.
That cannot be done inside the mechanism: a target-specific record makes `output_support`
target-dependent, which is the §4.4 Alice/Bob counterexample — disjoint supports, TV = 1,
no guarantee. The safe inversion is to record at MAPPING time as an audit artifact that
never touches `output_support`, `baseline_model`, `beta_t`, `phi_t` or the tilt — the same
pattern as the REASON array. Reporting rule now in `target_fidelity`: the quantity has a
ceiling below 1 independent of `tau`, and must never be quoted without it.

## 2026-09-27 — A-DAT-cohort-support: alleles v1 cannot emit, and what v1 still needs externally

**Question (Dylan).** If v1 routes every block through the tilted sampler, what external
data is still required? And since the target is external and never enters the cohort,
how are target-unique variants handled?

**Result 1 — external dependencies shrink to two for the mechanism.**

| dependency | role in v1 | status |
|---|---|---|
| 1000 Genomes panel | `plink --blocks` boundaries only | present (30x GRCh38) |
| GRCh38 genetic map | `Delta_x` in cM for transition constants | **NOT ACQUIRED — blocks exit gate 3** |
| 1000 Genomes panel | evaluation: attack DB, AF/LD site sets, Beagle reference | present |
| PanGenie callset | **nothing in v1**; needed only to run the RESTRICT `missing_policy` baseline | optional |

Dropping anchors removes PanGenie from the mechanism's critical path entirely. Note also
that block boundaries do not affect correctness — the per-block factorisation gives the
same `tau` for any segmentation with `sum beta_b` = 1 — so plink blocks are inherited for
HEAD-TO-HEAD COMPARABILITY with PanMixer, not out of necessity.

**Result 2 — alleles declared in the VCF that no donor carries.** Measured from
`external/PanMixer/starting_data/chr21/pangenome.npy` (44 x 340,824 x 2) against
`num_alleles.npy`, counting support over all 88 haplotypes:

| quantity | count | denominator | share |
|---|---|---|---|
| declared allele slots | 767,324 | — | — |
| declared alleles with ZERO cohort support | 8,823 | 767,324 slots | 1.15% |
| sites with >= 1 unsupported declared allele | 6,809 | 340,824 sites | 2.00% |
| sites where the unsupported allele is REF | 4,885 | 340,824 sites | 1.43% |
| sites monomorphic across all 88 haplotypes | 4,713 | 340,824 sites | 1.38% |
| sites with no called haplotype at all | 1 | 340,824 sites | 0.00% |

- Columns: counts, then the denominator each share is taken over. Slots and sites are
  different denominators and must not be mixed.
- Rows: "declared" means the allele index is below `num_alleles[site]`, i.e. the VCF lists
  it. "Zero cohort support" means no haplotype of the 88 carries it.
- Synthesis: v1's output support is cohort-supported paths, so **the sampler cannot emit
  these 8,823 alleles.** On 1.43% of sites the unreachable allele is REF, so a target
  carrying reference there cannot be matched by any donor state. The consequence is a
  `phi_t` penalty — **pure utility loss, and bounded at about 2% of sites. No privacy
  consequence**: `Y_D` is fixed and target-independent, so common support (condition 1)
  holds regardless. Monomorphic sites are released deterministically, which is also safe
  because the determinism is identical for every target and cancels in the likelihood
  ratio; determinism only breaks the bound when it depends on the input.

**Framing to keep.** In PanMixer the target is DATA — a column inside the matrix every
cohort statistic is computed over, which is how departure (c) arises. In v1 the target is
a PARAMETER: it never enters the matrix, `support_D`, or the state set, and touches the
mechanism only through `phi_t` comparing the target's block against each donor's, inside `psi_t`. Adding it to the cohort would
break the `R_D` cancellation in Step 3 of the Theorem 1 proof — not a weaker bound, no bound.

**Three senses of "target-unique", only one of which is a real gap.** (i) An allele absent
from `G`: unreachable by the mechanism because `p = Map(g, G)` is a traversal of `G`; lost
at mapping time, a limitation of the setting. (ii) A combination unique to the target: the
ordinary case, and exactly what the mosaic produces — conditional on the transition
constants actually recombining (exit gate 3). (iii) An allele in `G` with no cohort
support: the measured 1.15% above.

## 2026-09-27 — A-DAT-lv-nesting: the biallelic SNP mask is blind to snarl level

**Question.** Does `build_biallelic_snp_mask.py` (Step 11) use the `LV` tag? No — and
no PanMixer Python reads the INFO column at all (`grep -rn "LV=\|INFO" --include=*.py`
over the pinned checkout returns nothing). Step 11 re-reads `pangenome.vcf.gz` directly,
so INFO is present in the line it parses; it splits with `maxsplit=5` and tests only REF
and ALT.

**Method.** One pass over `external/PanMixer/starting_data/chr21/pangenome.vcf.gz`,
parsing `LV=` out of INFO per record, cross-tabulated against
`biallelic_snp_mask.npy`. All 340,824 records carry an LV tag; none missing.

| snarl level | all records | share | inside the mask | share of mask |
|---|---|---|---|---|
| LV=0 (top level) | 309,335 | 90.76% | 252,039 | 91.83% |
| LV=1 | 28,029 | 8.22% | 20,153 | 7.34% |
| LV=2 | 3,228 | 0.95% | 2,141 | 0.78% |
| LV=3 | 222 | 0.07% | 119 | 0.04% |
| LV>=4 | 10 | 0.00% | 6 | 0.00% |
| **nested (LV>0)** | **31,489** | **9.24%** | **22,419** | **8.17%** |

- Columns: all chr21 records, versus the 274,458 records the mask marks True.
- Rows: `LV` is the snarl-tree level written by `vg deconstruct -a`, 0 = top level.
- Synthesis: **9.24% of records are nested.** This independently reproduces the 9.2%
  recorded on 2026-09-22 for records whose POS lies inside an earlier record's span —
  two different methods (tag-based here, coordinate-based there) agreeing to rounding,
  which is a useful cross-check on both.

**Why it matters.** The mask is the site set for the gap-score attack (Step 21), which
assigns each site a Hardy-Weinberg genotype probability and sums `-log` across sites,
treating loci as independent evidence. 8.17% of that set is nested, and a nested allele
is conditional on its parent traversal, so those sites are not independent evidence.
The published Methods DO state a conditional-independence approximation for nested
variants — `p_hat(a_nested) ~ p(a_parent) * p(a_nested)` — but on the SAMPLING side
only; nothing equivalent appears in the scoring or attack path. **Magnitude unmeasured:**
whether this materially moves the gap score is not established here, and the
correlation may be weak in practice. Structure only.

**For v1.** `LV` and `PS` are two of the fields the REASON array should carry (see
`missing_policy`): they are what distinguishes "this donor's parent allele does not
contain the child bubble" from the other four causes of a `-1`.

## 2026-09-26 — A-DAT-missingness-scope: what 11.31% and 0.35% are means OVER

**Why.** `docs/terms.md` (`missing_policy`) and `docs/plan.md` quoted "11.31% versus
0.35%" with no measurement entry behind either, and the primer then described 11.31%
as the mean "across every variant" / "once every variant is a chain position". A
second-round primer audit re-measured it; re-measured again here, independently.

**Method.** `external/PanMixer/starting_data/chr21/pangenome.npy` (44 x 340,824 x 2, int16), missingness
per record = share of the 88 haplotype entries equal to `-1`, averaged over the class.
Anchors = the 273,475 keys of `pangenome_to_thousand_g_alignments.pickle` (the PanGenie
overlap, Correction 5). Block membership from `blocks_dict.json` (100,757 entries, every
record covered exactly once). Target included (the matrix as shipped).

| position class | n | share of variants | mean missingness | records > 25% missing |
|---|---|---|---|---|
| all variants | 340,824 | 100.00% | 4.263% | 30,800 |
| anchors, all | 273,475 | 80.24% | 0.351% | 925 |
| non-anchors, all | 67,349 | 19.76% | 20.150% | 29,875 |
| variants in multi-variant blocks | 251,737 | 73.86% | 1.582% | 8,108 |
| non-anchors in multi-variant blocks | 29,966 | 8.79% | **11.306%** | **7,555** |
| anchors in multi-variant blocks | 221,771 | 65.07% | 0.268% | 553 |
| variants in blocks with >= 2 anchors | 243,327 | **71.39%** | 0.501% | 1,274 |
| anchors in blocks with >= 2 anchors | 221,319 | 64.94% | **0.264%** | 538 |
| variants in multi-variant blocks with <= 1 anchor | 8,410 | 2.47% | 32.853% | 6,834 |
| singleton blocks | 89,087 | 26.14% | 11.840% | 22,692 |

- Columns: n and share are record counts; missingness is the class mean of the
  per-record missing share; the last column counts records over a quarter missing.
- Rows: `anchors in blocks with >= 2 anchors` is what PanMixer's HMM actually steps
  over. `all variants` is what v1 steps over, because v1 sends singletons through the
  same sampler as T = 1 chains. `variants in blocks with >= 2 anchors` is what PanMixer
  SAMPLES by donor copying today; everything else is drawn from allele frequencies.
- Synthesis: **11.31% is the mean over the 29,966 non-anchor records inside
  multi-variant blocks** — the positions v1 adds INSIDE blocks PanMixer already chains.
  It is not the mean over every variant (4.26%), and 0.35% is the mean over all
  anchors, not over the anchors PanMixer's HMM steps over (0.26%). The honest
  like-for-like is **v1's positions 4.26% (30,800 over 25%) against PanMixer's HMM
  positions 0.26% (538 over 25%)**, about 16x (derived). The "32x" (11.31 / 0.35) is a
  ratio of two differently scoped means and is quotable only with both scopes named.

**Second finding — the donor-copied share today is 71.4%, not 73.9%.** 73.86% is the
share of variants in multi-variant blocks. But the sampler sends a multi-variant block
with <= 1 anchor to the allele-frequency route too (8,410 variants, 2.47%), so the share
PanMixer samples by copying a donor is **71.39%** — the figure the 2026-09-18 entry
already gives ("10.4% of blocks ... 71.4% of variants"). `docs/plan.md` v1 REPLACE
item 2 said 73.9%; corrected there the same day. The 8,410 are also the worst-behaved
class in the table (32.85% mean missingness), which is presumably why they have <= 1
anchor.

## 2026-09-26 — A-DAT-primer-provenance: four primer numbers that `docs/` never logged

**Why this entry exists.** A three-critic audit of the expanded primer flagged
`274,458`, `177,865` and `3,077.60` as fabricated, and `8,661` as invented precision,
because `grep` over `docs/` returned nothing for any of them. All four were already in
the primer as committed on 2026-09-23. Each was re-measured from its artifact before
anything was deleted. **All four are exact.** The defect was that they were measured
during a session and quoted in teaching material without being logged here. They are
logged now, so they are quotable.

**1. Biallelic SNP records on chr21: 274,458 of 340,824 (80.53%).** Measured by
loading `external/PanMixer/starting_data/chr21/biallelic_snp_mask.npy` (dtype bool,
shape (340824,)) and summing it. The mask is built by
`external/PanMixer/starting_data/scripts/src/build_biallelic_snp_mask.py`: a record
qualifies iff it has exactly one ALT and both REF and ALT are a single base in ACGT.
**Relation to the other counts:** 65,176 records are indel/structural (see the
2026-09-23 conflicts entry), so 275,648 records are single-base; 274,458 of those are
biallelic, so **1,190 are multi-allelic SNP records**. The two numbers are consistent;
one is a subset of the other. `238,052` (69.8%) is a DIFFERENT quantity — the sites
the gap-score attack keeps after the mask AND the 1000G match — and must not be quoted
as the mask count.

**2. The capacity-0.5 obfuscation run, HG00438, chr21, seed 123.** From
`logs/pm_obf.log` (sentinel `logs/pm_obf.status` = 0), verbatim optimizer and stacker
lines:

| metric | capacity 0.1, seed 123 | capacity 0.5, seed 123 |
|---|---|---|
| selected moves | 47,070 | 177,865 |
| utility_loss / max | 928.92 / 9288.85 (10.0%) | 3,077.60 / 9288.85 (33.1%) |
| pmi_gain (nats, as reported) | 115,315.84 | 120,443.42 |
| alleles changed / total | 68,335 / 650,541 (10.50%) | 106,549 / 649,242 (16.41%) |

- Columns: two capacities from the same run script, same target, same seed.
- Rows: `moves` and `utility_loss` are the optimizer's own log line; `changed` is
  the stacker's. The allele denominators differ slightly between runs (650,541 vs
  649,242); the log does not say why, and it is recorded rather than explained.
- Synthesis: raising the budget five-fold (0.1 -> 0.5 of max loss) spends only 33.1%
  and buys 4.4% more reported privacy (120,443.42 / 115,315.84 = 1.0445). The knapsack
  saturates. Seed 456 at capacity 0.1 gave 47,069 moves and 928.88 — seed-stable.

**3. PanMixer's Python: exactly 8,661 lines.** `git ls-files '*.py' | xargs cat | wc -l`
inside the pinned checkout. Every `.py` on disk is tracked, so the untracked count is
the same. `docs/data.md` previously rounded this to ~8,660.

**4. The published Equation (5) carries its minus sign.** Rendered
`paper/s41467-026-77591-0_reference.pdf` **PDF page 10** (not 9) at 300 dpi:
`L(g_i, g*_target) = − Σ log f_j`, the sum running over the loci j in `S(g_i, g*_target)`. The minus is typeset and
unambiguous. The "sign discrepancy" the primer once carried as unresolved was an
artefact of text extraction dropping the glyph (this PDF has no text layer at all), not
a defect in the paper. With the sign, the equation agrees with the paper's own prose
("emphasizes rare allele matches"). Same page states the read-mapping donors as
HG00138, HG00635, HG01112, HG02698 and **NA18853**, and the HMM constants
Ne = 10,000, r = 1.26, Δx in cM, d = Δx·Ne·r, p = e^(−d/n) + (1 − e^(−d/n))/n.

**Lesson.** "Not in `docs/`" means *unlogged*, not *false*. A checker that cannot
find a number must re-measure it from the artifact before recommending deletion; the
audit here would otherwise have deleted four true measurements and published a
retraction of a correct one. See the matching entry in `docs/bugs.md`.

## 2026-09-25 — A-DAT-missingness: `-1` conflates at least four different things

**Goal.** Dylan asked how a donor can fail to traverse a parent bubble, given that a
variant ought to have a value for every haplotype. Establish what `-1` actually means
before choosing a `missing_policy`.

**The answer for TOP-LEVEL variants is that it cannot.** Every haplotype path runs
along the GRCh38 backbone, so a top-level variant has an allele for everyone; the
2.25% missingness there is assembly gaps and conflicts, not structure.

**For NESTED variants it can, and routinely does.** A child bubble sits inside a
PARTICULAR parent allele. A haplotype taking a different parent allele never visits
that part of the graph, so there is no allele to report. It is NOT APPLICABLE, not
missing — and VCF has no way to say so, therefore `.`, therefore `-1`.

Clean measured example: parent row 37624 at POS 13,205,133 (REF 13 bp; five ALTs of
15, 1, 17, 11 and 13 bp) with child row 37625 at POS 13,205,136, three bases inside
it. The 87 haplotypes carrying parent alleles 0, 1, 3, 4 or 5 all have a child
allele. The one haplotype carrying **allele 2 — the 1 bp allele, i.e. the region is
deleted** — has `-1`. The child variant does not exist on that haplotype because the
sequence it sits in was deleted.

**How much of nested missingness this explains**, over 5,735 nested records with any
missing child call:

| cause | share |
|---|---|
| entirely inherited from a missing PARENT call | 12.1% |
| among haplotypes WITH a parent call, perfectly determined by which parent allele they carry | **66.1%** |

**A second, unrelated cause: assembly gaps.** Run lengths of consecutive missing
positions (strand 0, all 44 subjects): 28,052 runs, median 1, mean 20.0, max 23,478.
**54.4% are isolated single positions** (structural / not-applicable), while **349
runs of length >= 100 cover 448,283 entries** — those are contig breaks, not structure.

**A FIFTH cause, from the converter rather than the data.**
`external/PanMixer/tools/common/VCFtoNP.py` maps a HAPLOID genotype (one allele, no
`|`) to `-1` on strand 1. Measured on 150,000 chr21 records: **3.135% of sample GT
fields are haploid**, and 5.581% of the two allele slots are an explicit `.`. So some
of the matrix's missingness is a conversion decision, not a property of the data.

**The causes are distinguishable — until conversion throws them away.** `LV` and `PS`
separate "not applicable" from "inherited"; run-length separates assembly gaps; the
`CONFLICT` tag names conflicts; GT arity identifies the haploid case. But `VCFtoNP`
reads only `fields[1]` (POS) and `fields[9+i].split(':')[0]` (the GT subfield) — the
entire INFO column is discarded. By the time anything reaches the matrix, all five
look identical, and `1/support(v)` and `pmi[hap == -1] = 0` treat them alike.
**Recommendation for v1: emit an auxiliary REASON array during conversion.** Every
cause is recoverable at that moment and none is recoverable afterwards.

**So `-1` means at least five different things**: (1) not applicable, the parent
allele does not contain this bubble; (2) inherited from a missing parent; (3) an
assembly gap; (4) a CONFLICT, where the sample had multiple paths that disagreed.
A single `missing_policy` is being asked to cover all five, and they do not want the
same treatment: "not applicable" argues for renormalise, an assembly gap arguably
argues for something closer to wildcard. **Consider splitting the policy by cause —
the LV/PS tags and run-length both distinguish them, and both are currently discarded
at VCF-to-numpy conversion.**

**Provenance.** chr21, `external/PanMixer` @ `c182c38`, measured 2026-09-25.


## 2026-09-25 — Design note: if the chain comes from the graph, what happens to LD blocks

**Verified fact.** PanMixer never loads a genetic map. Searching the whole codebase,
every `recomb` hit is an output FILENAME for the RECOMB conference; the only genetic-map
references are the unwired `external/PanMixer/starting_data/scripts/get_genetic_maps.sh`
(which has a broken shebang, `!/bin/bash`, and is referenced only by `.gitignore`).
Its transition "distance" is a VCF row index.

**That explains the architecture.** With no recombination map, the LD blocks have to
carry the recombination information structurally: the per-block HMM treats everything
inside a block as one linked unit and everything across blocks as independent. Blocks
are a discretised stand-in for a distance-dependent transition model.

**Consequence for our design.** In a Li-Stephens copying model the linkage is carried
by the DONOR PANEL and the TRANSITION PROBABILITIES, not by where block boundaries
fall — copying a stretch from donor j reproduces whatever allele combinations j
carries, including rare co-occurrences that exist only because they sit on one
ancestral chromosome. So with a genome-wide chain plus distance-dependent transitions,
LD blocks become redundant for the MODEL. The two jobs they do for PanMixer both
disappear for us: they scope its per-block HMM (we have one chain), and they are the
unit of its knapsack (we have no knapsack).

**Practical requirement this creates.** Distance-dependent transitions need a genetic
map in GRCh38 coordinates. The map repo PanMixer's unwired script points at
(the joepickrell 1000-genomes-genetic-maps repo) offers `interpolated_from_hapmap` and
`interpolated_OMNI`, both GRCh37-era. SHAPEIT4 publishes b38 maps (URL responds 302),
which is the likely source. **Not yet acquired or verified — an open task, not a
decision.** Falling back to bp distance with a constant rate is possible but crude.

**What LD blocks would still be for.** Evaluation comparability with PanMixer
(reporting per-block results on their partition), and possibly grouping for the
`beta_t` weights. Neither is a modelling need.


## 2026-09-25 — A-DAT-anchors: why anchors exist, and whether we need them

**Goal.** Dylan asked why the HMM runs only over anchor variants rather than every
variant in a block, and whether we can drop the concept. Establish the rationale
before deciding.

### The paper's rationale is about TOP-LEVEL SNPs, not external matching
Published Methods, on choosing the branch: *"dependent on the number of top-level
SNPs that exist inside the LD block. If 2 or more exist, enough SNPs exist to build
an HMM, and we can synthesize more accurate new haplotype blocks. If only one
top-level SNP exists, we do not have enough useful information to build an HMM."*
Top-level means on the GRCh38 backbone path (LV=0). The >=2 threshold is simply
that a Markov chain needs two positions to have a transition between them.

### Why restrict to top-level SNPs — the data supports it
Missingness per record class on chr21 (fraction of the 88 haplotypes with no called
allele):

| class | n | mean missing | records >50% missing |
|---|---|---|---|
| top-level SNP (the paper's V_SNPs) | 253,000 | **2.25%** | 421 |
| top-level non-SNP | 56,335 | 3.11% | 722 |
| **nested (LV>=1)** | 31,489 | **22.53%** | 3,515 |

Nested variants are ~10x more often uncalled, because a donor whose path does not
traverse the parent bubble has no meaningful allele at the child — and the matrix
cannot distinguish "not applicable" from "missing". Using them as chain positions
would feed the emission model noise. **That is a real justification.**

Two candidate rationales that do NOT hold:
- **Compute.** Anchors 273,475 vs all variants 340,824 — only **1.25x** more
  positions (24.1M vs 30.0M operations at K=88). Negligible.
- **Multi-allelic emissions.** Median ALT count is 1 even for non-SNP records; the
  90-allele tail is real but rare.

### But the CODE does not implement the paper's set
| | records |
|---|---|
| paper's V_SNPs (LV=0 and SNP) | 253,000 |
| code's anchors (PanGenie overlap) | 273,475 |
| overlap | 234,532 |
| **anchors that are NOT top-level SNPs** | **38,943** (14.2% of anchors) — of which **13,436 are NESTED** and 25,507 are non-SNP |
| top-level SNPs that are NOT anchors | 18,468 (7.3%) |

So the implementation admits 13,436 nested variants as chain positions — precisely
the class the paper's own rationale excludes — while dropping 18,468 top-level SNPs
that qualify. The PanGenie overlap is a proxy for "well-behaved position" that is
neither necessary nor sufficient.

### Implication for our design
Drop the PanGenie-overlap criterion: it has no principled basis and it is what makes
rare non-anchor variants invisible to `eps_j` (measured earlier: in 32% of blocks the
ignored non-anchor rarity exceeds `eps_j` itself). Do NOT simply promote every
variant to a chain position — the 22.5% missingness at nested records is a real
obstacle. The principled replacement is to take positions from the **graph
structure** (top-level / snarl-tree) rather than from an external callset, which is
the same answer open question 8 is circling for `T`.

**Provenance.** chr21, `external/PanMixer` @ `c182c38`, published paper Methods,
measured 2026-09-25.


## 2026-09-24 — A-DAT-noop-moves: 31% of selected obfuscation moves change nothing

**Method.** Walked block 1077 end to end for HG00438 on chr21, then swept the whole
capacity-0.1 run comparing each selected move's released alleles to the original.

**Block 1077** (rows 34816-34820, POS 13,000,333-13,000,621): 5 variants, all 5
anchors, so the HMM branch. Target carries `[0,0,0,0,0]` on both strands.
- `eps_j` = **0.0711 nats** per strand, computed from the target's ORIGINAL block.
- support(v) = 87 of 88 for every variant, so 1/support = 0.01149 each.
- The sampled replacement came back `[0,0,0,0,0]` — **identical to the target**,
  because the degenerate transitions copy one donor verbatim and that donor happens
  to carry the same common alleles.
- Changed variants = 0, so **`eta_j` = 0.00000** (it would be 0.05747 if all 5 changed).

So the move offers positive privacy at **zero** utility cost. The LP takes it
unconditionally. Confirmed in the run: `xsol[625]` is 1 on both strands, and the
released alleles at those rows are unchanged.

**Swept across the whole run: of 47,070 selected moves, 14,526 changed NOTHING
(30.9%).** Every one contributed its `eps_j` to the reported `pmi_gain` of
115,315.84 while altering no allele and costing no utility.

### Follow-up — the inflation is 7% of eps, not 31%
The 31% figure counts MOVES; it is not the share of privacy credit. Measured on a
3,000-move sample of the 47,070 selected:

| selected moves | n | median eps_j | mean eps_j |
|---|---|---|---|
| changed something | 2,112 | 2.1904 | 3.3103 |
| **changed nothing** | 888 | **0.3898** | 0.5964 |

**No-op moves carry only 7.0% of the total eps.** They concentrate on LOW-eps blocks
— common patterns that many donors share, so a randomly drawn donor often carries
the target's own alleles. Block 1077 is typical: 82 of 88 donors carry exactly the
target's pattern, giving eps_j = 0.0711, and -log(82/88) = 0.0706 confirms the
reading of eps_j as the surprisal of the target's block under the cohort.

Note also what is NOT inflated: for a no-op, `eta_j` is correctly **0**, because it
sums `1/support(v)` only over variants that differ. So the utility axis is honest;
it is the privacy numerator that gains ~7% for nothing, and the privacy-per-utility
ratio that is consequently overstated.

**Consequence.** `pmi_gain` — the numerator of the reported "Privacy Risk" axis —
is inflated by moves that did not touch the genome. The effect is systematic rather
than incidental: it follows directly from `eta_j` being computed as
`sum(utility_loss * (new != orig))` while `eps_j` is always positive and is computed
from the original block regardless of what was sampled. It compounds the sampler
degeneracy: the more often the sampler returns the target's own alleles, the more
free privacy the LP books.

**Provenance.** chr21, HG00438, capacity 0.1, seed 123, `external/PanMixer` @
`c182c38`, measured 2026-09-24.


## 2026-09-24 — A-DAT-af-table: where f_v comes from, and two corrections

**Goal.** Dylan asked what population the allele frequencies in `eps_j` are computed
from, and whether "no anchors" implies a degenerate frequency.

**The table is built by a three-way priority** per variant
(`external/PanMixer/starting_data/scripts/src/get_af.py`:27-78). Measured on chr21:

| branch | source of counts | target included? | variants |
|---|---|---|---|
| 1. variant matches the 1000G map | **1000G only** | not under Phase 3; **YES under our 30x repin** | 258,610 (75.9%) |
| 2. matches PanGenie but not 1000G | PanGenie **plus all pangenome haplotypes** | yes | 25,325 (7.4%) |
| 3. novel (neither) | pangenome only, then REF padded to `NOVEL_DENOMINATOR` | yes | 56,889 (16.7%) |

`NOVEL_DENOMINATOR = 2 * (N_1000G + N_pangenome)` = 2*(3202+44) = **6,492** under our
repin, so a novel allele carried by one haplotype gets f_v = 1.54e-4, NOT 1. The
padding exists precisely to stop a novel allele looking common because only its own
carriers were counted; the REF allele absorbs it.

**"Anchor" and "has population data" are different things.** The anchor set is the
PanGenie map; the AF priority also consults the 1000G map. **10,460 non-anchor
variants still receive 1000G frequencies.** A block with zero anchors can be fully
populated with real frequencies.

Median −log f_v where HG00438 carries a NON-REF allele: branch 1 = 0.63 nats,
branch 2 = 0.43, branch 3 = **5.41**. At novel sites the target usually carries REF
(median f_v 0.998), so most novel sites contribute ~0 — the signal is in the 11,760
novel sites where it carries a non-REF allele.

### Correction 7 — the `#FIX ME` zeroing does NOT fire on our data
[INCORRECT] - Alleles like these have population frequency 0 in any external callset, so `-log(0)` -> infinity -> zeroed by that `#FIX ME`. **PanMixer assigns the most identifying sites on the chromosome a privacy score of exactly zero.**
[CORRECTION - 2026-09-24]: Measured for HG00438 across all 319,092 chr21 sites where it carries a called allele: **ZERO have f_v = 0**, so the `#FIX ME` branch never executes. The reason is that the target's own allele is counted in every branch — via the pangenome counts in branches 2 and 3, and (under our 30x repin) via the panel itself in branch 1, since 39 of 44 HPRC donors including HG00438 are in the 30x panel. The zeroing is a LATENT defect whose trigger is panel-dependent: under the shipped Phase 3 panel, where 0 of 44 HPRC donors appear, a branch-1 variant whose allele no 1000G sample carries WOULD give f_v = 0 and fire it. The conclusion that the hypervariable site contributes zero privacy is still correct, but the MECHANISM is different — see Correction 8.

### Correction 8 — the hypervariable site scores zero by non-anchor exclusion, not by zeroing
The 257 kb record (row 51791) gets a perfectly good frequency: HG00438 carries
alleles 89 and 88, each f_v = 1.540e-4, i.e. **8.78 nats** of self-information. That
value is never used. The record is **not an anchor** (`r in align` is False), and its
block 1639 has 92 records with 82 anchors, so the block takes the **HMM branch** —
where the forward algorithm scores only the anchors. The most identifying record on
the chromosome contributes nothing because it is not in the PanGenie callset, which
is a completely different failure from the `-log(0)` zeroing.

### New risk — our own repin weakened target-independence
Branch 1 draws counts from the 1000G panel alone, and the code comment asserts "None
of the pangenome subjects appear in 1000g_phased" — true for Phase 3 (0 of 44
overlap, measured) but **false for our 30x substitution (39 of 44 overlap)**. So our
rebuilt `allele_frequencies.npy` includes the target's own alleles in the
frequencies used to score it. For PanMixer that is merely a stronger version of an
existing self-inclusion. [INCORRECT] - **For US it violates constraint 1**, so any leave-one-out
cohort must also drop the target from the 1000G panel, not only from the pangenome.
[CORRECTION - 2026-09-28]: The measurement is right and the implication is wrong. 39 of 44
HPRC donors including HG00438 really are in the 30x panel — that is true of the pipeline as
built, and it matters for the BENCHMARK, where a cohort member stands in as target. It is
not a constraint-1 violation of our MECHANISM, whose target is external to every panel by
construction. Two further reasons this is now moot: under model C the mechanism uses no
external panel at all, and it computes no allele-frequency table (decided 2026-09-28). The
overlap still needs removing from the ATTACK DATABASE — and note the fix is incomplete, because
all 39 overlaps lie in the 698 samples the 30x release ADDED beyond Phase 3 to complete trios,
so dropping the 39 named samples leaves their first-degree relatives in.

**Provenance.** chr21, HG00438, `external/PanMixer` @ `c182c38`, measured 2026-09-24.


## 2026-09-23 — A-DAT-conflicts: what CONFLICT records do downstream, and how blocks handle spanning variants

**Goal.** Dylan asked whether all CONFLICT records behave alike, what they imply for
blocking, and what actually happens in the released output when a long variant and
records inside its span fall in different blocks.

### Finding 1 — a CONFLICT becomes MISSING DATA, and missingness is not neutral
At a CONFLICT record vg writes that sample's genotype as `.|.`, so `VCFtoNP` stores
**-1**. Measured on chr21:

| | value |
|---|---|
| CONFLICT records | 1,855 (0.54%) |
| mean haplotypes missing at a CONFLICT site | **29.4 of 88** |
| mean haplotypes missing elsewhere | 3.6 of 88 |
| matrix entries missing overall | 1,278,653 (4.26%) |
| sites with >=1 missing entry | 52,926 (15.5%) |

Missingness then drives BOTH of PanMixer's per-move scores, in the same direction:
- **Privacy -> 0.** `external/PanMixer/starting_data/scripts/src/get_support_and_pmi.py` zeroes the score wherever the target's
  haplotype is -1 (`pmi_1[haplotypes[:,0] == -1] = 0`).
- **Utility cost -> up.** support(v) counts non-missing entries, so 1/support(v)
  rises. Median utility weight at CONFLICT sites 0.01818 vs 0.01136 elsewhere —
  **1.60x more expensive to change**.

**Consequence: these sites are worth zero privacy and cost more, so the knapsack
will essentially never select them. They are released VERBATIM.** That is a leak
concentrated exactly where the target's structure is unusual.

### Finding 2 — conflicts are NOT uniform; they concentrate in complex regions
| property | CONFLICT records | all records |
|---|---|---|
| nesting LV=0 | 59.1% | 90.8% |
| nesting LV=1 | **37.5%** | 8.2% (a 4.6x enrichment) |
| nesting LV=2 | 3.4% | 0.9% |
| median size (max REF/ALT) | 4 bp (mean 602, max 311,406) | mostly 1 bp |
| median alleles per record | 4 (max 81) | 2 |

They also land in the SIMPLE part of the block structure: the 1,855 records occupy
1,634 blocks, of which **1,469 are singletons** and only **1.8% are HMM-path blocks**
(baseline 10.4%). So they are mostly alone in their own block and mostly handled by
the independent per-variant allele-frequency sampler, not the haplotype model.

### Finding 3 — the read-mapping evaluation never sees any of this
`external/PanMixer/tools/downstream/utility_out/vg_prep.py` builds the evaluation graph with `bcftools view -v snps`, which keeps
288,020 of 340,824 chr21 records (84.5%) and drops the 52,804 indel/structural
records where conflicts and nesting concentrate. So PanMixer's read-mapping utility
metric is computed on a graph that excludes the hardest regions by construction.
`af_loss`, `ld_loss` and the gap score read the matrix/VCF directly and DO see them.

### Finding 4 — what the released output actually is when a variant spans blocks
Read from `external/PanMixer/tools/panmixer/stacker.py`:37-70. The output is `new_haplotypes.npy`,
shape (n_sites, 2) — **one allele index per RECORD**. It is built as
`new_variants = original_haplotypes.copy()`, then for each selected (block, strand)
the record indices in that block are overwritten from the presampled replacement.

So there is no sequence assembly step, and no coherence check. A long parent record
and the records positionally inside its span are decided **independently**: if the
parent's block is not selected the parent keeps the target's TRUE allele, while
child records in a selected block are replaced. The released VCF then asserts both
at once. Whether those assignments jointly correspond to a walk through the graph is
never computed, checked, or asked.

### Finding 5 — the MECHANISM sees everything; only the EVALUATION filters
An earlier note implied indel/structural records are excluded. That is wrong as a
statement about the mechanism. Measured on our HG00438 chr21 run at capacity 0.1:

| stage | records seen | non-SNP | CONFLICT |
|---|---|---|---|
| mechanism (PMI + knapsack + stacker) | 340,824 (100%) | 65,176 | 1,855 |
| `af_loss` / `ld_loss` (strict 1000G mask) | 258,610 (75.9%) | 21,550 | 37 |
| gap score (posref mask) | 268,851 (78.9%) | 29,965 | 110 |
| gap score after the biallelic-SNP mask | 238,052 (69.8%) | **0** | 49 |
| read mapping (`bcftools view -v snps`) | 288,020 (84.5%) | mixed records kept | **1,167** |

So indels and structural variants ARE scored, selected and rewritten. In that run
**15,206 non-SNP records changed — a 23.3% rate, HIGHER than the 14.9% for SNPs** —
and the 257 kb hypervariable record was rewritten on both strands
(alleles [89, 88] -> [35, 64]).

**Sequence actually rewritten: 3,129,810 bp = 3.13 Mb, 6.7% of chr21.** 222 changed
alleles exceed 1 kb; the largest single rewrite is 257,485 bp. PanMixer reports the
same run as "68,335 of 650,541 alleles changed (10.5%)" — because `1/support(v)` is
per RECORD, the 3.13 Mb never enters the accounting.

Correction to an intermediate figure used while investigating: dropping non-SNPs
does NOT remove most conflicts. `bcftools view -v snps` keeps mixed multi-allelic
records, and conflicts are enriched in those (median 4 alleles), so **1,167 of
1,855 conflicts survive** the read-mapping filter — 37% are dropped, not 84%.

**Net effect on the scoring.** The biggest edits by base pairs land on exactly the
records the read-mapping metric cannot evaluate (its graph is SNP-only), while
`af_loss` and `ld_loss` see 75.9% of records and the gap score's biallelic path
sees no non-SNP records at all. No metric in the suite is sensitive to the 3.13 Mb.

**Provenance.** chr21, `external/PanMixer` @ `c182c38`, measured 2026-09-23.


## 2026-09-22 — A-DAT-hypervariable-sites: measured the worst-case site class, flagged as a risk and an opportunity

**Goal.** Dylan asked whether hypervariable structural sites could affect how well
our mechanism performs, and whether they are somewhere we could improve on
PanMixer. Measure the class before reasoning about it.

**Method.** Found the most multi-allelic site on chr21 and characterised it and its
class from the VCF, the graph traversals (`AT` field) and the numpy matrix.

**Result — the worst site on chr21 is `grch38#chr21:14569980`, graph bubble
`>102270111>102277685`, LV=0 (top level).**

| property | value |
|---|---|
| alleles declared | 90 (89 ALT) |
| distinct alleles carried by the cohort | **88, across 88 haplotypes — every haplotype unique** |
| REF length | 257,439 bp |
| ALT lengths | 257,124 – 514,768 bp |
| size of the single VCF line | **27,983,949 bytes (~28 MB)** |
| `AC` for allele 1 | 0 — an orphan left by `view -s ^chm13` without `--trim-alt-alleles` |
| nested child records (LV>=1, PS = this bubble) | **2,464** |
| total matrix rows for this one region | **2,485** (0.7% of chr21's 340,824 rows) |
| graph nodes in the bubble | 7,575 distinct; traversals 4,971–10,008 nodes |
| allele 58 | 10,008 nodes over 5,014 distinct = ~4,994 revisits -> a **tandem duplication** of the whole region |
| inversions | none (0 traversals use reverse orientation) |

Class frequency on chr21: **7** sites with >=88 alleles, **155** with >=44, **618**
with >=20.

**Why `vg deconstruct -a` emits both levels.** The nested decomposition DOES happen
— 2,464 children exist. The top-level record is emitted IN ADDITION because some
alleles cannot be expressed as a combination of the children: allele 58 traverses
the region twice, and no set of independent per-site substitutions can encode "this
whole region is duplicated". VCF alleles are also flat strings with no syntax for
"allele = this combination of child variants", so the top-level traversals must be
spelled out in full. Hence 28 MB.

**Provenance.** `external/PanMixer/starting_data/chr21/pangenome.vcf.gz` ·
`num_alleles.npy` · `pangenome.npy` · measured 2026-09-22.

### Why this matters to US — four consequences, recorded as risks to watch
1. **Double counting.** The same sequence is represented twice: once in the
   top-level record's giant allele strings and again across its 2,464 children. All
   2,485 are rows in the matrix, so a target's variation here contributes to
   `eps_pmi` at BOTH levels. The paper's nested approximation
   p(nested) ~ p(parent) x p(nested) addresses child-vs-parent dependence but not
   the fact that the parent is itself a scored row.
[INCORRECT] - 2. **The most identifying sites score ZERO.** Every allele here has population frequency 0 in any external callset, so `-log(0)` -> infinity -> zeroed by the `#FIX ME` in `get_support_and_pmi.py`. The sites where each haplotype is unique — a perfect fingerprint — contribute nothing to PanMixer's privacy score.
[CORRECTION - 2026-09-26]: The conclusion stands and the mechanism is wrong. The zeroing branch never fires on this data (Correction 7: 0 of 319,092 sites have f_v = 0, because the target's own allele is counted in every frequency branch). These sites score zero because the privacy score reads anchor alleles only and they are not anchors (Correction 8). Marked in place on 2026-09-26; the copy under Correction 7 was made without marking this original.
3. **Utility weighting ignores scale.** `1/support(v)` is per VARIANT. A move here
   rewrites up to half a megabase yet is costed like a one-base SNP (and, since
   support is 1, is among the most expensive — but for the wrong reason).
4. **Array width is set by the worst site.** `allele_frequencies.npy` is
   (340824, 90): every site is padded to the widest one in the chromosome.

### Opportunity for OUR mechanism — to evaluate, not yet a claim
These sites are exactly where a bounded, path-level guarantee should beat a
per-block edit rule, because our `tau` bound holds regardless of how rare an allele
is, whereas `eps_pmi` collapses to 0 precisely here. Three things must be checked
before claiming any advantage:
- whether `output_support` can contain a *valid* traversal here at all (the
  AF-sampling branch draws variants independently and would produce combinations no
  haplotype carries — see the graph-validity open question);
- whether `phi_t` should be length-weighted rather than per-variant, given a single
  block can span 257 kb;
- how to avoid double-counting the top-level record and its 2,464 children in both
  `u_path` and any privacy measurement.
**Nothing here is measured about our mechanism — it does not exist yet.**

### Follow-up 2026-09-22 — the overlap is STRUCTURAL, not an edge case
Dylan asked for the fraction of chr21 affected. Measured over all 340,824 records,
taking each record's span as [POS, POS + len(REF) - 1]:

| measure | value |
|---|---|
| records whose POS lies inside an earlier record's span | 31,487 (**9.2%**) — matches the LV>=1 nesting fraction exactly |
| ...inside a record of REF length >= 10 kb | 19,211 (5.6%) |
| ...inside a record of REF length >= 100 kb | 4,608 (1.4%) |
| blocks containing >=1 overlapped record | 15,540 of 100,757 (**15.4%**) |
| variants living in an affected block | 78,937 (**23.2%**) |
| HMM-path blocks (>=2 anchors) affected | 1,229 of 10,502 (**11.7%**) |

**Verdict: not documentable as an edge case.** Roughly a tenth of records and a
quarter of variants-by-block are involved, so block disjointness — which both the
knapsack's additive accounting and our FFBS chain assume — fails broadly, not locally.

**The duplication's structure, measured.** Allele 58 of the parent bubble walks
10,008 node steps over 5,014 distinct nodes; 4,994 nodes appear exactly twice and 20
appear once. At step ~5004 the walk reaches `>102277684` (beside the closing
boundary) and jumps back to `>102270112` (beside the opening boundary). So it is a
**whole-region tandem repeat**, and the two copies are NOT identical — the 20
singleton nodes are ~10 interior positions where the copies diverge. That is exactly
why the interior variants cannot be given one allele per haplotype for this sample,
and why `CONFLICT=HG01891` appears at POS 14577263.

**Why VCF cannot express this compactly.** The graph already stores it compactly —
the traversal simply visits the nodes twice, duplicating no sequence. The ~28 MB is
purely a VCF export artifact: a VCF allele is a flat string, so a cyclic traversal
must be linearized. And a duplication CANNOT be re-expressed as a plain insertion
without loss, because the inserted copy carries its own variation and **VCF positions
are reference coordinates** — there is nowhere to place a variant that lives inside
inserted sequence. This is the limitation pangenome graphs exist to remove.

### Clarification 2026-09-22 — CONFLICT is mostly NOT caused by revisiting
Dylan asked whether the two-alleles-at-a-child-variant situation arises whenever a
repeat revisits a bubble. Measured on chr21:

| | count |
|---|---|
| records whose own AT contains a revisiting traversal | 37 (0.011%) |
| records carrying a CONFLICT tag | 1,855 (0.54%) |
| both | 11 |
| revisit but NO conflict | 26 |
| **conflict but NO revisit** | **1,844** |

CONFLICT names all 45 samples, roughly 200-470 records each — far too widespread to
be explained by tandem duplications, of which only 37 records chr21-wide show any
sign. The dominant cause is therefore something else; the most plausible candidate
is assembly fragmentation (each HPRC haplotype is many contigs, and two contigs of
one sample traversing the same region give that sample multiple paths), but **this
cannot be determined from the VCF alone and has not been verified.**

⚠ **Limitation of that test.** A record's `AT` lists the traversals of ITS OWN
bubble, so it cannot reveal that a sample passed through that bubble twice — that
is only visible at the enclosing parent. So the cross-tabulation UNDERSTATES the
link for child records inside a duplicated region. The mechanism is real; its
genome-wide share is not what the 11/37 overlap suggests.

Note the asymmetry it exposes: the 257 kb parent HAS a revisiting allele but is NOT
flagged CONFLICT, because at the parent level HG01891 has exactly one allele (58)
whose walk happens to loop. The conflict only materialises at CHILD records, where
that single parent allele corresponds to two passes that may disagree.

**Design implication for `T`.** A snarl-tree segmentation would make this one region
ONE block instead of 370 overlapping ones, and blocks disjoint by construction. The
usual objection — that snarls are not LD blocks — is weak for OUR design, because in
a Li-Stephens copying model the linkage is carried by the donor panel, not by the
block boundaries; and switching donors at snarl boundaries is graph-valid by
construction, which is what `output_support` needs anyway.


## 2026-09-18 — A-EVL-readmapping: reproduced the published read-mapping BASELINE, close but not exact

**Goal.** The sharpest available test of our installation. The paper's read-mapping
experiment is chr21-only and uses five EXTERNAL donors, so scope cannot excuse a
miss, and the BASELINE (unobfuscated graph) needs no obfuscation sweep.

**Method.** Extracted chr21 reads for the paper's five donors (HG00138 EUR,
HG00635 EAS, HG01112 AMR, HG02698 SAS, NA18853 AFR) from the 1000G 30x CRAMs;
built the graph per `vg_prep.py` and aligned single-end per `quick_align.py` with
vg 1.68.0; scored with `vg stats -a` using the metric definitions in
`external/PanMixer/plot_util/plot_helpers.py` (each of perfect / gapless / MAPQ60 over total aligned).

**Result — the right regime, but a systematic offset.**

| metric | published baseline | ours | difference |
|---|---|---|---|
| perfect | 77.83% | 79.05% | +1.22 |
| gapless | 95.62% | 95.96% | +0.34 |
| MAPQ60 | 77.01% | 79.23% | +2.22 |

Supporting signals that the pipeline itself is sound: ~10.9-11.4M reads aligned per
donor, matching the scale of the authors' own sample output preserved in
`external/PanMixer/tools/downstream/utility_out/giraffe_parser.py` comments (total
aligned 11,521,845; that sample works out to 77.68% perfect / 95.82% gapless, close
to the published figures, so those comments come from this pipeline). That sample
also shows `Total paired: 0`, confirming the authors mapped SINGLE-END as we did.
The per-donor ordering is biologically sensible: NA18853 (AFR) is lowest on every
metric, as expected for an African genome against a reference-biased graph.

**Provenance.** chr21 graph VCF 287,147 SNP records · vg 1.68.0 · 5 donors ·
`logs/pm_vg.log`, `logs/pm_vg2.log` · read slices in
`external/PanMixer/read_fastqs/`.

### Correction 6 — my explanation for the offset was wrong, and the test refuted it
[INCORRECT] - the offset is likely caused by the reference build: I used chr21 from the 1000G GRCh38 analysis set, which hard-masks 2.05 Mb more of chr21 (8,671,409 N vs 6,621,364 N) than the UCSC hg38 chr21 that clean_fa.py implies; masking the false-duplication regions removes mapping ambiguity and should inflate perfect alignment and especially MAPQ60
[CORRECTION - 2026-09-18]: Tested by rebuilding the entire graph and re-aligning all five donors against UCSC hg38 chr21. The result is essentially IDENTICAL: 79.05 / 95.96 / 79.23 (UCSC) versus 79.03 / 95.96 / 79.26 (analysis set) — a shift of at most 0.03 points, against an offset of 1.2-2.2. The hypothesis is refuted. The reason is clear in hindsight and should have been predicted: the reads were extracted from CRAMs that were themselves aligned to the analysis set, so no reads originate in the extra-masked regions under EITHER graph. The masking difference was real but could not matter. **The measured N-count difference stands; only the causal claim was wrong.**

### Hypotheses tested, and what remains unexplained
1. **Reference build masking** — REFUTED by direct rebuild (above).
2. **`convert_2_vcf` genotype rewriting** — MEASURED, at most a minor contributor.
   That converter rewrites any sample field lacking `|` to `.|.`
   (`external/PanMixer/tools/common/convert_2_vcf.py`:71-73), which would shrink allele support and
   drop variants at the `-c 1` step, lowering mapping rates. On 200,000 chr21
   records, 2.38% of sample fields lack `|`, but **zero records lose all support**,
   so few if any variants would be dropped. Direction is right; magnitude is not.
3. **`vg stats` version** — UNTESTABLE, and the attempt produced a finding of its
   own: the published Methods say alignment used Giraffe v1.68.0 and the GAMs were
   analysed with `vg stats` **v1.36.0**. Those two are incompatible. Running the
   real v1.36.0 binary on a v1.68.0 GAM aborts with "obsolete, invalid, or corrupt
   protobuf input". The stated combination cannot have been used as written.

[INCORRECT] - **VERDICT.** The installation runs the pipeline correctly and lands in the right regime with the right internal structure, but this is **NOT an exact reproduction** of the published baseline, and the residual +1.2 to +2.2 point offset is **unexplained**. Do not describe the read-mapping baseline as reproduced.
[CORRECTION - 2026-09-26]: Verdict unchanged; the range was misstated. The table directly above gives +1.22 (perfect), +0.34 (gapless) and +2.22 (MAPQ60), so the offset spans **+0.34 to +2.22** points. "+1.2 to +2.2" dropped the gapless column. Caught by the 2026-09-26 primer audit.

**What this does NOT establish.** Nothing about whether the paper's numbers are
wrong. An offset of this size is consistent with an undocumented difference in read
extraction or filtering, in the exact graph VCF the baseline was built from, or in
the tool versions actually used — all invisible from the paper and the released code.


## 2026-09-18 — A-EVL-published-validation: checked our run against the PUBLISHED paper

**Goal.** Dylan supplied the PUBLISHED PanMixer paper (Nature Communications,
DOI 10.1038/s41467-026-77591-0, "Article in Press", 12 pages). All prior analysis
in this project was made against the bioRxiv preprint. Validate our installation
against what the authors actually report, and re-check our findings against the
version of record.

**Method.** Extracted the published text, diffed it against the preprint, and ran
a 9-agent workflow (4 readers + adversarial verifiers + synthesis) over both
versions and the code. Then hand-verified every load-bearing claim.

**VERDICT: the installation is running correctly**, in the precise sense that it
faithfully executes the released code at `c182c38` and produces internally
coherent results. Four independent signals, none guaranteed in advance:

| check | published | ours | verdict |
|---|---|---|---|
| cohort size | 44 individuals | 44 / 88 haplotypes | match |
| capacity -> normalized utility loss | budget is a fraction of max loss | 928.92/9288.85 = 10.000% at capacity 0.1 | exact |
| per-subject cost, longest chromosome | 5 min / 14 GB (Table 1) | 68 s / 1.7 GB on chr21; x5.47 to chr2 = 6.2 min / 9.3 GB | consistent |
| determinism / invariants | not reported | byte-identical per seed, 0 allele-bound violations | pass |

**Provenance.** `paper/s41467-026-77591-0_reference.pdf` (sha256 `2bbba736...`) ·
`external/PanMixer` @ `c182c38` · chr21 on the 30x GRCh38 panel · HG00438 ·
workflow `wf_c1f52889-672`, 9 agents, 0 errors.

**What the verdict does NOT say.** We have reproduced NO reported result. Every
headline number — `eps_private` = 0.002, average utility loss 0.28, AF divergence
0.006/0.006/0.005, LD 0.003, read mapping 77.82/95.61/77.01 — is unrun. "Correct
installation" means "we run their code as shipped", NOT "the code implements the
published method"; we have four confirmed places where it does not.

### Correction 4 — the "47 individuals" discrepancy was preprint-only
[INCORRECT] - The cohort is **44 individuals / 88 haplotypes**, not the paper's "47" (the HPRC VCF carries 45 samples, one being `chm13`, and the pipeline drops it and chrX).
[CORRECTION - 2026-09-18]: The PREPRINT says 47; the PUBLISHED paper says 44 (page 1 line 20, and "n = 44" in the Fig. 2 caption and Table 1). There is no discrepancy against the version of record — the published cohort size matches our measurement exactly. Recorded because the original line reads as though the authors' paper disagreed with their own data, which is true only of the superseded preprint.

### Correction 5 — the HMM's anchor set is the PanGenie overlap, not the 1000G overlap
[INCORRECT] - Anchors are pangenome records that exactly match a PanGenie "bi_all" record ... the HMM-vs-AF routing is therefore exposed to the GRCh37/GRCh38 build mismatch along with everything else derived from the 1000G match set
[CORRECTION - 2026-09-18]: The first half is right and the implication is wrong. `hmm.py:32-33` loads `pangenome_to_thousand_g_alignments.pickle`, and in `get_mappings.py:28` the variable `thousand_g_alignments` is bound to **`PG.vcf.gz`** — the PanGenie callset — NOT to 1000G. The name is a trap. Because PanGenie is GRCh38, the per-block HMM-vs-AF ROUTING is NOT degraded by the build mismatch. What the mismatch does corrupt is the LD block BOUNDARIES (plink runs on the 1000G panel), the attack database, and the `af_loss`/`ld_loss` site sets. Two different anchor notions; do not conflate them. The paper meanwhile defines anchors as V_1000G (page 8 lines 944-945), which the code does not implement.

### The build inconsistency is in the PUBLISHED paper, not only the released code
The published Methods pin V_1000G to "the Phase 3 phased 1000 Genomes dataset" as
a SUBSET of V_0, the GRCh38 backbone (page 8 lines 941-945), and the NEW Data
availability section pins the exact GRCh37 EBI URL, release 20130502 (page 11).
The backbone is stated as GRCh38 three times. The full text contains zero
mentions of GRCh37, liftover, 30x, high-coverage or NYGC. Our 1,025 vs 258,610
chr21 match measurement shows the stated containment cannot hold as written. The
preprint had no Data availability section, so the published version sharpens the
inconsistency rather than resolving it. Note also that the official GRCh38
liftover of Phase 3 was WITHDRAWN in 2021, with 1000 Genomes redirecting users to
the 30x GRCh38 collection — the panel we substituted, and the one the repo's own
unwired `get_blocks_grch38.sh` downloads.

### Numbers that CHANGED between preprint and published — use the published ones
- `eps_private`: 0.001 -> **0.002**. Its DEFINITION also changed, from "the
  smallest privacy risk value at which all targets could no longer be linked" to
  "the **largest**" — the preprint's wording was backwards.
- Wasserstein AF divergence at `eps_private`: 0.004 / 0.004 / 0.002 ->
  **0.006 / 0.006 / 0.005** (all alleles / SNPs / rare SNPs).
- The published version adds a membership-inference analysis absent from the
  preprint (0 mentions in the preprint, 4 in the published version).
Any reproduction target must come from the published version.


## 2026-09-18 — A-DAT-grch38-repin: chr21 rebuilt on GRCh38; PanMixer executed and measured

**Goal.** Make PanMixer runnable here, on build-correct data, and replace the
inferred claims from the ingest with measurements on real chr21.

**Method.** Wrote a blocking local `sbatch` shim (`scripts/sbatch`) so the pinned
checkout runs UNMODIFIED, and re-pointed the hardcoded `1000g_phased.vcf.gz` slot
at the 30x GRCh38 panel. Ran the shipped preprocessing for chr21, then
`obfuscate.py` for HG00438 at two capacities and two seeds, then diagnostics.

**Result — the repin works.** `get_mappings`, PanMixer's own matching code, reports
**258,610 strict / 268,851 relaxed** matches, against 1,025 / 2,176 on the wired
GRCh37 panel. That is exactly the number an independent `bcftools` comparison
predicted, so two independent routes agree. The attack database is now
(3202, 258610, 2) instead of one built on 1,025 spurious matches.

**Result — PanMixer runs and is well-behaved.** 68 s per subject-chromosome, peak
RSS 1.7 GB. Same seed gives byte-identical `new_haplotypes.npy` and `xsol.npy`;
different seeds differ. Output shape (340824, 2), **zero** allele-bound violations,
VCF correctly bgzipped and indexed. Capacity is respected: 0.1 -> 10.0% utility
loss. At capacity 0.5 it spends only 33.1%, i.e. the knapsack saturates — all
positive-value moves fit — so capacity is not a tight knob at the top end.

**Result — block structure.** 14,137 plink LD blocks but 100,757 entries in
`blocks_dict.json`, **88.4% singletons**. Only **10.4% of blocks** take the HMM
path, though those hold **71.4% of variants**. [INCORRECT] - For us: `T_blocks` ~ 100,757 on
chr21 with `K_states` = 88, so FFBS is ~9M operations.
[CORRECTION - 2026-09-28]: The block structure is correct; reading it as OUR chain length is
not. Model C was adopted 2026-09-28 — one chain position per VARIANT — so T is 340,824 raw
and 306,480 after the `chain_span` cut, roughly 3x the figure above, and FFBS is ~30M
[CORRECTION - 2026-09-29]: 306,480 was the ONE-SIDED cut (p arm only). The telomeric end
was closed on 2026-09-29 for the same reason — 592 variants above the map end formed a
593-position zero-distance run — so the figure is **305,887** (340,824 - 34,345 - 592),
89.75%. The FFBS cost is unchanged in magnitude.
operations. `K_states` = 88 is right and stands. LD blocks play no role in our mechanism at
all under C.

**Provenance.** `external/PanMixer` @ `c182c38` · chr21 rebuilt on the 30x GRCh38
panel (sha256 `a925c112...`) · HG00438 · capacity 0.1 and 0.5 · seeds 123, 456 ·
`logs/pm_chr21.log`, `logs/pm_obf.log`.

### Confirmed on real data — the sampler is inert
300 of 300 randomly chosen HMM-path blocks returned a block **identical to a single
donor haplotype**. The synthetic-panel result from the ingest holds on real chr21.

### Correction 3 — I overstated the eps_pmi defect, and the measurement refutes it
[INCORRECT] - the forward panel's self-inclusion (obfuscate.py:130,132) has flattened the privacy score, so eps_j is pinned near log(2N) = 4.4773 nats for essentially every HMM block regardless of true rarity
[CORRECTION - 2026-09-18]: Measured on 150 real HMM-path blocks for HG00438: eps_j runs 0.071 to 4.508 with median 1.242, and only 9/150 sit within 0.05 nats of log(88) = 4.4773. It is NOT pinned. (The value 4.4773 itself remains correct — it is log(2N) and is still the ceiling; what was retracted is the claim that eps_j SITS there.) What self-inclusion actually does is impose a CEILING at ~log(2N): because the target always matches itself, p(h_bj) >= 1/2N, so eps_j <= ~4.48 however rare the block truly is. Re-scoring the same blocks with the target REMOVED from the emission panel gives median 1.299 but max 78.1. For the 9 ceiling-bound blocks the leave-one-out value has median 12.06 and max 78.12 — an understatement of up to 17.4x. So the defect is real but its shape is the opposite of "flattened": typical blocks are almost unaffected (median ratio 1.03x) while the RAREST blocks — the most identifying, exactly the ones a privacy mechanism must prioritise — are truncated hardest. The knapsack therefore systematically under-values the highest-risk blocks.

**What this does NOT establish.** Still nothing about whether the PAPER's numbers
are affected: eps_j enters the LP only through a ranking, and we have not measured
how much the selection changes when the ceiling is removed. Do not claim their
frontier is wrong; claim only that the shipped scorer truncates rare blocks.


## 2026-09-18 — A-DAT-panmixer-ingest: PanMixer pinned, its data acquired, and its mechanism mapped

**Goal.** Bring the upstream dependency in, pin it, run it, and establish how much
of it RanPanMixer can build on — the proposal names PanMixer as the source of its
target-independent prior, its blocks, and its utility information.

**Method.** Cloned G2Lab/PanMixer, pinned `c182c38d5bc8bb6f00f4b0b101207c4a009ca045`
(the only release, 2026-06-06, MIT). Relocated the checkout to `/data` because
PanMixer writes inside its own tree. Built its `environment.yaml` env. Downloaded
its three input datasets (16 GB) and its two external tools. Fetched the paper
(bioRxiv DOI 10.64898/2026.02.16.706152, CC-BY-NC-ND 4.0). Then ran a 22-agent
workflow: 8 readers over the paper and six code subsystems, each followed by an
independent adversarial verifier, plus 4 checks of what OUR proposal asserts about
PanMixer, a synthesis and a completeness critic. Finally re-verified every
load-bearing claim by hand — reading the cited lines, executing the shipped HMM
class on a synthetic panel, and measuring cross-build variant overlap directly.

[INCORRECT] - **Result.** The cohort is **44 individuals / 88 haplotypes**, not the paper's "47" (the HPRC VCF carries 45 samples, one being `chm13`, and the pipeline drops it and chrX).
[CORRECTION - 2026-09-26]: 44 / 88 is right and matches the PUBLISHED paper exactly; "47" is the superseded preprint only (Correction 4). Marked in place on 2026-09-26; the copy under Correction 4 was made without marking this original.
Two findings dominate, both independently reproduced:

1. **The released preprocessing pipeline matches GRCh38 against GRCh37.** Measured
   on chr21: 1,025 of 340,824 pangenome records match the wired 1000G Phase 3 panel
   exactly on (POS,REF,ALT) — 0.30%. The repo ships an unwired GRCh38 path.
2. **The Li-Stephens sampler never recombines.** Four compounding defects, present
   identically in the standalone tool AND in the pipeline copy the paper's
   precomputation imports. Executing the shipped class: total switch mass 5.8e-12
   at a 1-index gap; the sampled block equals exactly one donor haplotype. With the
   paper's own Fig. 6 constants it would be 0.13-0.99.

Both are in `docs/bugs.md` with evidence. The consequence for us is the headline:
**"reuse PanMixer's cohort HMM as R_D" is not available.** What we inherit is its
DATA, its BLOCK structure, its forward recursion's algebra, and its EVALUATION
suite — roughly 55-65% of ~8,660 lines by volume, but under 10% of the mechanism.

**Provenance.** `external/PanMixer` @ `c182c38` · HPRC v1.0 PGGB GRCh38 VCF
(sha256 `ead25541...`) · 1000G Phase 3 chr21 (sha256 `1942e070...`) · PanGenie
Zenodo 7669083 (sha256 `e8c0c48d...`) · paper sha256 `99ac1bd2...` · workflow run
`wf_a30688a6-3d7`, 22 agents, 0 errors.

**What this does NOT establish.** Nothing about RanPanMixer's own mechanism — no
line of it exists yet. It also does not establish that the PanMixer PAPER's
published numbers are wrong: the defects are in the RELEASED code, and the unwired
GRCh38 script plus the code branches preferring 3,202-sample artifacts suggest the
authors' own runs used a path the release does not wire in. We have not reproduced
their figures, so we cannot say which code produced them. **Do not state or imply
that their results are invalid.** What we CAN say is narrower and sufficient: the
released pipeline as documented does not reproduce the paper, and its sampler as
shipped is not the model the paper describes.

### Control — the build mismatch is measured, not inferred

Downloaded the GRCh38 30x panel (the one the repo's unwired `get_blocks_grch38.sh`
fetches) and repeated the identical comparison against the same pangenome chr21
records. Same pangenome, same method, only the panel's build differs.

| metric | Phase 3 (wired, GRCh37) | 30x (unwired, GRCh38) |
|---|---|---|
| exact (POS,REF,ALT) matches | 1,025 | 258,610 |
| relaxed (POS,REF) matches | 2,176 | 268,851 |
| match rate of 340,824 pangenome records | 0.30% | 75.88% |

252x. The deficit is the reference build, not the data. The residual ~24% is
expected — the pangenome carries graph-only and structural variants absent from a
SNV/INDEL/SV panel — but has not been characterised and should not be assumed benign.

**Still not established:** which code produced the paper's numbers. The unwired
GRCh38 script and the branches preferring 3,202-sample artifacts are consistent
with the authors having run a corrected path, but that is an inference, not a
measurement. Do not write that PanMixer's results are invalid.

### Correction 1 — the proposal's characterization of PanMixer is wrong in four places
These were checked against paper and code by dedicated agents and by hand.

[INCORRECT] - The HMM implementation by PanMixer uses a transition model with one probability for remaining in the same donor state and another probability for switching donors, and the per-target running time is O(TK).
[CORRECTION - 2026-09-18]: The stay/switch FORM is real and its forward recursion is genuinely O(K) per step (obfuscate.py:140-152). But PanMixer's SAMPLER is O(A*K^2) in time and memory — `get_transition_probabilities_without_subject` materialises every dense K'xK' matrix for a block up front (:117-124). More importantly the chain runs WITHIN a block over anchor SNPs, with a fresh uniform start state per block and NO transition model between blocks; there is no T-blocks x K-states chain, no backward pass and no FFBS anywhere in the repo. Our across-block chain and backward sampler are new work, not inherited.

[INCORRECT] - PanMixer protects individuals contributing to the pangenome by selecting and applying block-level haplotype obfuscations. The proposed mechanism protects an external target individual's mapped path without modifying the graph.
[CORRECTION - 2026-09-18]: "Without modifying the graph" is not the contrast. PanMixer does not modify graph topology either — its release is the same cohort VCF with one sample column rewritten. The real contrast is WHO is protected (a cohort member vs an external individual) and WHETHER the release is written back into the reference.

[INCORRECT] - PanMixer-compatible LD blocks or graph subpaths and the cohort HMM can serve as the target-independent prior R_D.
[CORRECTION - 2026-09-18]: Only the block BOUNDARIES (plink on an external panel) and the transition FORM are target-independent as shipped. The allele-frequency table includes the target at PanGenie-matched and novel sites (get_af.py:41-53,68-73); the support counts include the target (obfuscate.py:208); the forward/PMI panel includes the target (:130,132); and the site-to-block assignment ranges over sites private to the target. Every one must be recomputed leave-one-out before R_D satisfies the theorem's hypotheses.
[CORRECTION - 2026-09-28]: The catalogue above is accurate about PANMIXER and the final sentence is wrong about US. It reads as a precondition on our R_D; it is not. Our target is EXTERNAL to G, D, every panel and the attack database, so none of those artifacts can contain it and there is nothing to recompute: K = 88, and constraint 1 holds by construction. Leave-one-out is a requirement of the HEAD-TO-HEAD ARM ALONE, where a cohort member stands in as target and externality must be simulated. This entry is directly contradicted by the 2026-09-25 entry in this same file ("In v1 the target is a PARAMETER: it never enters the matrix, `support_D`, or the state set"). See CLAUDE.md constraint 1 and `loo_cohort` in docs/terms.md, both corrected in commit 1f8ea3f.

[INCORRECT] - PanMixer samples a candidate constrained to differ from the target's own block.
[CORRECTION - 2026-09-18]: The PAPER says this (p.10 lines 653-654); the released CODE implements no such constraint — `sample_block_prior` is an unconditional prior draw with no comparison to the original and no rejection loop. So its candidate generator is already target-independent given the panel. PanMixer's target-dependence comes instead from the knapsack SELECTION and from releasing unselected blocks verbatim. Note the consequence: because a no-op resample is possible while eps_j stays positive, the LP can buy privacy credit at zero utility cost from moves that change nothing.

### Correction 2 — an error this session made and caught
[INCORRECT] - the shipped transition code gives a total switch mass of 1.4e-12 at a 1-index gap, quoted as the sampler's behaviour
[CORRECTION - 2026-09-18]: That figure belongs to `get_transitions` (the FORWARD pass, no factor 4). The SAMPLER uses `get_transition_matrix_without_subject`, which includes `* 4 *` and gives 5.8e-12. Measured directly. The degeneracy conclusion is unchanged, but the number was attributed to the wrong function. Relatedly, the first draft claimed the forward pass "omits" a factor 4 that the sampler "includes", implying 4 is correct; paper Fig. 6 gives `d = Delta_x * Ne * r` with NO factor 4, so the SAMPLER is the deviation.

### Decision — the integration seam is `new_haplotypes.npy`
Every PanMixer attack and utility evaluator consumes one artifact: a
`(n_sites, 2)` int16 array, `-1` = missing, row-aligned to the chromosome's VCF.
Our sampler will emit exactly that file, making it effectively a third stacker
strategy, so their entire benchmark runs on our output unchanged. This is the
single highest-leverage design decision from this session and it removes most
implementation-difference confounds from the head-to-head.

### Decision — the head-to-head is leave-one-out, aligned on measured attack success
`eps_pmi` and `tau` are formally incomparable, and so are `capacity` and `tau`.
The comparison aligns the arms on MEASURED empirical attack success and reports two
frontiers. A hard limit to state in any write-up: the PGGB graph topology was built
from all HPRC assemblies, so even a full leave-one-out only APPROXIMATES our
external-target threat model on their data.

### Finding — their cohort-level utility metrics cannot see target fidelity
A pure-prior draw (our tau = 0) is a real cohort haplotype: it would score
near-perfectly on U_AF, U_LD and external read mapping while carrying zero
information about the target. Any comparison must carry a `target_fidelity` axis or
our maximal-privacy point looks like a free lunch. The critic established that such
a measurement already exists and is mechanism-agnostic
(`MIA_privacy.py:80-83,140`), correcting the synthesis's claim that it had to be built.

