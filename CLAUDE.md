# RanPanMixer — Agent Contract

A randomized **path release mechanism**: given a public pangenome graph `G` built
from a public cohort `D`, release an external target's mapped path in sanitized
form, with the guarantee that any two input paths induce output distributions at
total variation distance at most a pre-specified `tau`.

This file is the **contract**: mission, environment, hard rules, conventions.
It carries **no findings, no results, no numbers** — those live in `docs/`.

---

## Project Core

- **Mission:** implement and validate the cohort-tilted `tau`-private path release
  of `paper/Private_Genome_Path_Release.pdf`, on top of PanMixer's cohort HMM.
  Success = a sampler whose released paths retain measurable utility while the
  proven `tv_bound` holds and no implemented attacker beats `p_succ_bound`.
- **What this is NOT:** PanMixer protects cohort members by editing the graph.
  This protects an *external target's mapped path* and never modifies the graph.
- **Environment:** single Linux workstation. No data acquired yet — every path in
  `docs/data.md` is `PLANNED` until it is built and its row is updated.
- **Stack:** Python. NumPy/SciPy for the sampler. PanMixer is pinned at `c182c38`
  under `external/PanMixer` (a symlink to `/data` — it writes inside its own tree
  and assumes SLURM, which this machine does not have). Its conda env is `panmixer`.
  Check before assuming any library is installed.

## Documentation Architecture

- `CLAUDE.md` — **auto-loaded every session.** This contract. **Zero findings.**
- `docs/plan.md` — **read at session start.** Current state, active work, ledger.
- `docs/terms.md` — **search before naming any metric.** Definitions + provenance.
- `docs/data.md` — **read before touching data, runs or paths.**
- `docs/memory.md` — **never read on startup; `grep` it.** Chronological ledger.
- `docs/bugs.md` — **`grep` by symptom.** [Symptom] -> [Proven fix].

Also present: `paper/` (the proposal and manuscript drafts), `archive_docs/`
(read-only source material), `tests/` (code tests + the docs contract).

## Standing Constraints (hard rules; rationale lives in `docs/memory.md`)

1. **The output support `Y_D` and the baseline `R_D` are built from `G` and `D`
   only — never from the target.** Any target-dependent pruning of the HMM state
   space voids Theorem 1. This is the easiest way to buy utility and silently
   lose the guarantee; it is the single rule most likely to be broken by accident.
   ⚠ On PanMixer's data, target inclusion is the DEFAULT — its allele frequencies,
   support counts and scoring panel all include the target, because there the
   target is a cohort member. Recompute leave-one-out; see `docs/data.md`.
2. **The utility must be globally bounded in [0,1].** If a raw utility has range
   `delta_u != 1`, the calibration is `eta_tau = arctanh(tau)/delta_u`. An
   unnormalized utility does not weaken the bound — it removes it.
3. **Sample; never argmax.** Viterbi/MAP decoding is not a draw from the release
   distribution and carries no guarantee. The canonical decoder is `ffbs`.
4. **One release per genome, cached and reused.** Independent re-releases of the
   same path compose and degrade the guarantee. Never resample per query.
5. **`tau`, the utility function and the endpoints are pre-registered before any
   attack outcome is inspected.** They do not change mid-phase.
6. **An attacker that fails is not evidence the bound holds.** Measured attacker
   accuracy is a LOWER bound on distinguishability. Only the proof gives the upper
   bound; a weak attacker proves nothing about the mechanism.
7. **Never mix two `tau` values, two utility functions, or two cohorts in one
   comparison.** Utility is only comparable at matched `tau`; privacy is only
   comparable at matched utility.
8. **The guarantee covers the released path only.** It does not protect the raw
   target genome, the graph, or any other released statistic about the same person.

## Agent Behavior & Prohibitions

- **Verify before acting.** Never guess paths, shapes or counts. Use `ls`/`grep`/`docs/data.md`.
- **No sweeping changes.** No global `sed` or multi-file refactor without permission.
- **Strict documentation limits.** Do not create new doc files; these six are the set.
- **Prevent definition drift.** Search `docs/terms.md` before naming a quantity.
  Never invent a synonym for something already named.
- **Report the failure, not the workaround.** A missing tool or unmet assumption
  must never silently become a negative result.

## IMPORTANT: Results Reporting Format

**Metrics are ROWS. Configurations are COLUMNS.** Never the transpose — a wide
table of arms invites cherry-picking which metrics to show.

1. Report the full metric set every time, including rows that did not move.
   `n/a` with a reason; never an omitted row.
2. **Every number is stamped with its `tau`, its denominator and its `n`.** A
   utility or attacker figure without a `tau` is uninterpretable, not merely
   imprecise.
3. Row labels are the exact `docs/terms.md` identifier. Mark gate metrics with `*`.
4. Order rows by importance: primary -> gates -> structure -> context -> demoted.
5. Carry a provenance line: `<GRAPH> · <COHORT> · <TAU> · <UTILITY> · n draws`.
6. State the ceiling (`tau -> 1`, the untilted target) and the floor (`tau = 0`,
   a draw from `R_D` alone) as their own columns.
7. Every table is followed by THREE things as BULLET LISTS: (a) one bullet per
   COLUMN, (b) one bullet per ROW, (c) a prose SYNTHESIS.
8. EVERY NUMBER IN PROSE MUST BE TRACEABLE TO A TABLE CELL, or name its source.

## Execution & Long-Running Tasks

- Anything over ~60s runs detached with a status sentinel:
  `tmux new-session -d -s <n> '<cmd> > logs/<n>.log 2>&1; echo $? > logs/<n>.status'`
- Poll the sentinel. Never report a run finished without reading it.
- **Any fan-out MUST assert its own completeness and fail loudly.** A script that
  filters results cannot tell "found nothing" from "ran nothing".
- Never `pkill -f <pattern>` where the pattern matches your own command line.

## Naming

- **Work items:** `<phase>-<KIND>-<slug>`, `KIND` in `THY IMP DAT EVL ATK FIX LCK`,
  slug lowercase-hyphenated and meaningful. Current phase letter: `A`.
  Example: `A-IMP-ffbs-sampler`.
- **Directories/files:** `<phase>_<TARGET>[_<variant>]`. Never two names differing
  only in case. Put the varying parameter — above all `tau` — IN the filename.
  No brace/glob shorthand in docs (write `tau_0p10/`, `tau_0p25/`).
- Deprecated things are renamed `DEPRECATED_*` or deleted, and the choice is
  recorded in `docs/data.md`.

## THE IN-PLACE CORRECTION RULE

`docs/memory.md` is permanent. **Never delete or overwrite a historical entry.**
When something in it is proven wrong:

1. Prepend `[INCORRECT] - ` to the original line, preserving its text verbatim.
2. Insert directly below: `[CORRECTION - YYYY-MM-DD]: ` with the new finding.

Preserving the wrong version is what makes the reasoning legible later. Grep
`[INCORRECT]` to list everything this project has been wrong about.

**Corrections do not propagate themselves.** After writing one, hunt its stale
copies — grep the distinctive *number or phrase*, not the topic, across `docs/`,
`src/`, `tests/`, `paper/`, `CLAUDE.md`, and the assistant memory directory.

## YOU MUST EXECUTE: The Wrap-Up Protocol

Run when any unit of work completes — not at the end of the project.

1. **Archive to `docs/memory.md`.** Hypothesis, method, provenance, result.
2. **Compress to the ledger.** One row in `docs/plan.md`.
3. **Document fixes.** `[Symptom]` -> `[Proven fix]` into `docs/bugs.md`.
4. **Define new terms.** Full entry in `docs/terms.md`.
5. **Register artifacts.** A `docs/data.md` row for every new output.
6. **Reset the board.** Rewrite `docs/plan.md` Current State; bump `Last updated`.
7. **Verify.** `python -m pytest tests/test_docs_contract.py` must pass.
