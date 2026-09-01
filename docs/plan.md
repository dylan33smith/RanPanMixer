# plan.md — the board

**Last updated:** 2026-08-29

Read at session start. This file exists so a new session never has to grep
`docs/memory.md` to know where things stand.

---

## Current state

- **Branch:** `main`. Remote `origin` -> `https://github.com/dylan33smith/RanPanMixer`.
- **Phase:** `A` — foundation. Goal of the phase: an end-to-end sampler that is
  provably correct on a graph small enough to enumerate exhaustively.
- **Status: specification only.** The repository contains the proposal, this
  documentation set, and the docs verifier. **There is no source code, no graph,
  no cohort and no target genome.** No number has been produced by this project,
  so no number is quotable from it.
- **What is established** comes entirely from `paper/Private_Genome_Path_Release.pdf`
  and is theory, not measurement: tilting a target-independent cohort model
  `baseline_model` by `exp(eta_tau * u_path)` bounds the total variation distance
  between the release distributions of any two input paths by `tau`, and therefore
  bounds any equal-prior binary attacker at `p_succ_bound`. Theorem 1 additionally
  gives the stronger pointwise likelihood-ratio bound `epsilon_tau`. With PanMixer's
  `stay_switch` transitions, `ffbs` sampling costs O(TK) per haplotype.
- **The binding constraint right now** is that no graph or cohort has been chosen.
  That decision fixes `T_blocks`, `K_states`, the block segmentation, and hence the
  entire runtime and utility story. Everything else in the backlog is downstream
  of it, except the toy-scale correctness work, which is deliberately not.
- **The one thing most likely to go wrong** is not the proof — it is an
  implementation that quietly makes `output_support` or `baseline_model` depend on
  the target, and thereby reports excellent utility under a guarantee that no
  longer holds. Nothing observable in the output distinguishes that case. It has
  to be prevented structurally, not detected afterwards.

## Headline result

**None.** No experiment has been run. This section stays empty until
`A-THY-toy-enumeration` produces the first measured table; it is not a placeholder
for expectations.

## Reporting contract

Two results are comparable only if their config stamps match on:
**graph+cohort id · `tau` · utility id (`phi_t` k and w, `beta_t` scheme)**.

Metrics are rows, configurations are columns. Every table carries a provenance
line and is followed by one bullet per column, one bullet per row, and a prose
synthesis. Every number in prose traces to a table cell. Emit tables from the
scoring script; never hand-assemble one.

The floor is `tau` = 0 (an untilted draw from `baseline_model`, target-independent
by construction). The ceiling is the target's own path, `u_path` = 1. Both are
columns, not prose asides.

## Ledger — Phase A

| ID | Intervention | Endpoint | n | Result | Verdict | memory |
|---|---|---|---|---|---|---|
| A-FIX-docs-system | Adopt the six-file documentation system and write the verifier | verifier passes on a greenfield tree | n/a | 6 files + 20 checks, 19 active | done | 2026-08-29 |

**Provenance:** no runs. The single ledger row is documentation work and produced
no measurement.

## In Progress

Nothing is in progress. The next item is `A-DAT-graph-cohort`.

## Backlog

Ordered. Each item names what must pass before it starts and what result closes it.

### [A-DAT-graph-cohort] Select the graph, cohort and block segmentation — NEXT
- **Why:** fixes `T_blocks`, `K_states`, `output_support` and every runtime claim.
  Nothing downstream is meaningful before it, and choosing it late means
  re-deriving all of it.
- **Prerequisite:** none.
- **Exit gate:** a `docs/data.md` row per artifact with state `OK`, stating n
  cohort individuals, `K_states` before and after duplicate collapse, `T_blocks`,
  and the segmentation policy. An external target genome held out of `D` and
  demonstrably not used in graph construction.
- **Decide with it:** whether to reuse PanMixer's LD blocks or segment ourselves.
  Reusing them is the paper's stated plan and the cheaper path.

### [A-THY-toy-enumeration] Brute-force correctness gate on a toy graph
- **Why:** this is the ONLY test that can falsify the implementation of Theorem 1.
  On a support small enough to enumerate, `release_distribution` can be computed
  exactly and compared to the sampler's empirical distribution, and the true TV
  between two input paths can be computed rather than estimated. At scale neither
  is possible — `tv_empirical` is biased upward and cannot falsify anything.
- **Prerequisite:** `A-IMP-cohort-hmm`, `A-IMP-utility`, `A-IMP-ffbs-sampler`.
  Deliberately NOT dependent on `A-DAT-graph-cohort`: it runs on a synthetic
  toy graph and should be built before any real data arrives.
- **Exit gate:** (a) sampler empirical distribution matches exact enumeration
  within Monte-Carlo error; (b) exact TV between induced distributions is <= `tau`
  across a grid of `tau` and many input-path pairs; (c) `log_z_p` lies in
  [0, `eta_tau`] on every draw; (d) `tau` = 0 reproduces `baseline_model` exactly.
- **Kill criterion:** any exact-TV value above `tau` on the toy support. That is
  an implementation bug — the theorem is proven — and everything stops until found.

### [A-IMP-cohort-hmm] The target-independent prior
- **Why:** `baseline_model` is the object the whole guarantee rests on.
- **Prerequisite:** none for the toy version.
- **Exit gate:** rho and A built from `G` and `D` alone, with a test that fails if
  any target-derived quantity reaches the constructor. Zero-probability transitions
  encode graph-incompatible concatenations. `stay_switch` forward values match a
  dense-transition reference to numerical tolerance.

### [A-IMP-utility] Bounded utility
- **Why:** an unbounded `u_path` removes the guarantee rather than weakening it.
- **Prerequisite:** none for the toy version.
- **Exit gate:** `phi_t` in [0,1] and `u_path` in [0,1] asserted at runtime, not
  documented; `beta_t` summing to 1 asserted; a recorded `utility_id` capturing
  k and w; `diploid_utility` using MEAN, with a test that catches a SUM.

### [A-IMP-ffbs-sampler] Algorithm 1
- **Why:** the mechanism.
- **Prerequisite:** `A-IMP-cohort-hmm`, `A-IMP-utility`.
- **Exit gate:** per-position normalization (no underflow at realistic `T_blocks`);
  `log_z_p` accumulated from the forward normalizers and range-asserted; O(TK)
  `stay_switch` path agreeing with the O(TK^2) dense path; no argmax anywhere.

### [A-LCK-preregister] Freeze the evaluation before looking at it
- **Why:** the proposal is explicit that the utility function must be chosen
  before inspecting target-specific attack outcomes. Choosing it afterwards
  makes every subsequent privacy claim unfalsifiable.
- **Prerequisite:** `A-THY-toy-enumeration` passing.
- **Exit gate:** a committed file pinning the `tau` grid, the utility id, the
  primary endpoint and the attacker set, dated, before any attack is run.

### [A-EVL-utility-curve] Utility versus tau
- **Why:** the practical question the mechanism exists to answer.
- **Prerequisite:** `A-LCK-preregister`.
- **Exit gate:** `utility_retained` across the pinned `tau` grid, with floor and
  ceiling columns, reported on BOTH measurement paths (all blocks and switched
  blocks only) separately labelled — see `docs/terms.md`.

### [A-ATK-baseline] A concrete attacker
- **Why:** a one-sided implementation check. `attacker_accuracy` above
  `p_succ_bound` proves a bug; below it proves nothing.
- **Prerequisite:** `A-LCK-preregister`.
- **Exit gate:** at least one attacker that is strong on the `tau` -> 1 control
  (so a null result is powered rather than uninformative), and no attacker
  exceeding `p_succ_bound` at any pinned `tau`.

## Blocked

| Item | Blocked on |
|---|---|
| everything involving real data | `A-DAT-graph-cohort`: no graph, cohort or target genome selected |
| any runtime or scaling claim | `T_blocks` and `K_states` are unknown until the graph is chosen |
| PanMixer component reuse | no checkout, no pinned commit — see `docs/data.md` |

## Dropped — with reasons

| Item | Why |
|---|---|
| Viterbi / MAP decoding of the tilted HMM | Not a draw from `release_distribution`; carries no guarantee while producing plausible-looking output. Retired in `docs/terms.md` before use. |
| Per-target pruning of the HMM state space | Makes `output_support` target-dependent and voids Theorem 1. The proposal flags it as needing a separate privacy analysis; it is out of scope for Phase A. |

## Open questions

1. **Which graph and cohort.** See `A-DAT-graph-cohort`. Everything waits on it.
2. **What `tau` is defensible.** `tau` = 0.10 bounds equal-prior identification at
   0.55. Whether that is an acceptable operating point is a policy question this
   project can inform but not settle; it must be pre-registered either way.
3. **How to certify `output_support` independence.** Constraint 1 in `CLAUDE.md`
   is currently enforced by discipline. It should be enforced by construction —
   ideally the target is simply not in scope where `baseline_model` is built.
4. **What `phi_t` should be.** The proposal offers weighted k-mer Jaccard,
   PanMixer utility information, and mapping-seed preservation as candidates. The
   choice must be made and frozen before attack outcomes are inspected.
5. **Composition across releases.** Constraint 4 says one release per genome. What
   a system should do when the same target legitimately needs a second release is
   unanswered, and the current answer — "reuse the first" — is an operational
   constraint, not a theoretical result.
