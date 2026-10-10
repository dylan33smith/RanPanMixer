# RanPanMixer — Agent Contract

A randomized **path release mechanism**: given a public pangenome graph `G` built
from a public cohort `D`, release a target's mapped path in sanitized form such
that any two input paths induce output distributions at total variation distance
at most a pre-specified `tau`. This file is the **contract** — mission,
environment, hard rules, conventions — and carries **no findings, no results and
no numbers**; those live in `docs/`.

## Project Core

- **Mission:** implement and validate the cohort-tilted `tau`-private path release
  of `paper/Private_Genome_Path_Release.pdf`. Success = a sampler whose releases
  retain measurable utility while `tv_bound` holds and no implemented attacker
  beats `p_succ_bound`.
- **TWO ARMS, BOTH REQUIRED.** The paper has two sections, and the arm decides what
  is true of the target. **Never state a property of the target without naming its
  arm** — that confusion is what the 2026-09-28 scope correction was about.
  - **Arm 1 — the head-to-head.** The target is one of the 44 cohort members. We
    release the obfuscated paths of all of them, diploid, over the 22 autosomes,
    and run PanMixer's own evaluations so the comparison is direct. Requires
    **per-target leave-one-out**: BOTH of the target's haplotypes leave the donor
    panel. Every leave-one-out, reduced-`K` and cohort-distortion quantity belongs
    to this arm alone.
  - **Arm 2 — the external target.** The target is outside `G`, `D`, every panel
    and the attack database, so the full cohort is the donor panel and nothing is
    held out. Nothing in the mechanism is a function of the target but the tilt.
- **INDEPENDENT OF PANMIXER.** PanMixer is a published example of a problem in the
  same space — not a dependency, not a data source. We run our own
  `vg deconstruct` on the HPRC GFA we downloaded, and nothing in `src/`, `configs/`
  or `tests/` may read a path inside its tree. Arm 1 adopts its *question*, never
  its artifacts, and **both arms run on one consistent build** — the thing its own
  release got wrong.
- **Environment:** single Linux workstation. Every `docs/data.md` path carries its
  own state — read the row, do not assume.
- **Stack:** Python, NumPy/SciPy, env `ranpanmixer` (`environment.yml`; `pytest` IS
  in it). PanMixer is pinned under `external/PanMixer` (symlink to `/data`),
  read-only, for reference and Arm 1 only — its env is not ours. Check before
  assuming any library is installed.

## Documentation Architecture

- `CLAUDE.md` — **auto-loaded every session.** This contract. **Zero findings.**
- `docs/plan.md` — **read at session start.** Current state, active work, ledger.
- `docs/terms.md` — **search before naming any metric.** Definitions + provenance.
- `docs/data.md` — **read before touching data, runs or paths.** Every other
  directory in the repo is a row here.
- `docs/memory.md` — **never read on startup; `grep` it.** Chronological ledger.
- `docs/bugs.md` — **`grep` by symptom.** [Symptom] -> [Proven fix].

## Standing Constraints (hard rules; rationale lives in `docs/memory.md`)

1. **The output support `Y_D` and the baseline `R_D` are built from `G` and `D`
   only — never from the target.** Target-dependent pruning of the state space
   voids Theorem 1, and is the easiest way to buy utility and silently lose the
   guarantee. Arm 2 holds by construction; **Arm 1 holds only under
   leave-one-out**, and `G` itself still saw the target, so that arm's guarantee
   is conditional on the graph. State that; never claim it away.
2. **The utility must be globally bounded in [0,1].** If a raw utility has range
   `delta_u != 1`, the calibration is `eta_tau = arctanh(tau)/delta_u`. An
   unnormalized utility does not weaken the bound — it removes it.
3. **Sample; never argmax.** Viterbi/MAP decoding is not a draw from the release
   distribution and carries no guarantee. The canonical decoder is `ffbs`.
4. **`beta_v` sums to 1 over the whole scope `tau` covers — the 22 autosomes,
   genome-wide.** Each chromosome takes a SHARE of the budget. Normalising per
   chromosome silently makes 22 independent releases of one person.
5. **One release per genome, cached and reused.** Independent re-releases of the
   same path compose and degrade the guarantee. A `tau` grid is an evaluation
   device, never a release procedure. ⚠ Arm 1 releases 44 correlated paths whose
   donor panels are each other's true haplotypes: the per-release bound holds and
   the JOINT bound does not. Report that; do not quietly design around it.
6. **`tau`, the utility function, the endpoints, the attacker set and the attack
   database are pre-registered before any attack outcome is inspected.** They do
   not change mid-phase, and the file is dated and committed.
7. **An attacker that fails is not evidence the bound holds.** Measured attacker
   accuracy is a LOWER bound on distinguishability. Only the proof gives the upper
   bound; a weak attacker proves nothing about the mechanism.
8. **Never mix two `tau` values, two utility functions, or two cohorts in one
   comparison.** Utility is comparable only at matched `tau`, privacy only at
   matched utility. ⚠ Leave-one-out gives each Arm 1 target its own `beta_v`, so
   its 44 releases carry 44 utility functions — pool measured attack outcomes,
   never utilities.
9. **The guarantee covers the released path only.** Not the raw target genome, the
   graph, or any other released statistic about the same person.

## Agent Behavior & Prohibitions

- **Verify before acting.** Never guess paths, shapes or counts. Use `ls`/`grep`/`docs/data.md`.
- **No sweeping changes.** No global `sed` or multi-file refactor without permission.
- **Strict documentation limits.** Do not create new doc files; these six are the set.
- **Prevent definition drift.** Search `docs/terms.md` before naming a quantity.
  Never invent a synonym for something already named.
- **Report the failure, not the workaround.** A missing tool or unmet assumption
  must never silently become a negative result.
- **A gate that cannot fail is decoration.** Inject the error a gate exists to
  catch and assert which gates fire; a suite written beside its own code tests
  the author's model, not the spec.

## IMPORTANT: Results Reporting Format

**Metrics are ROWS. Configurations are COLUMNS** — never the transpose, which
invites cherry-picking. Report the full metric set every time (`n/a` with a
reason, never an omitted row). **Every number is stamped with its arm, its `tau`,
its denominator and its `n`** — one missing stamp makes a figure uninterpretable,
not imprecise. Row labels are exact `docs/terms.md` identifiers, gates marked `*`.
State the ceiling and the floor as columns, each against **its own arm's** ceiling.
Tables are emitted by the scoring script, carry a provenance line, and are followed
by one bullet per COLUMN, one per ROW and a prose SYNTHESIS; every number in prose
traces to a cell. Full checklist: `docs/plan.md`'s Reporting contract.

## Execution & Long-Running Tasks

- Anything over ~60s runs detached with a status sentinel:
  `tmux new-session -d -s <n> '<cmd> > logs/<n>.log 2>&1; echo $? > logs/<n>.status'`
- Poll the sentinel. Never report a run finished without reading it.
- **Any fan-out MUST assert its own completeness and fail loudly.** A script that
  filters results cannot tell "found nothing" from "ran nothing".
- Never `pkill -f <pattern>` where the pattern matches your own command line.

## Naming

- **Work items:** `<phase>-<KIND>-<slug>`, `KIND` in `THY IMP DAT EVL ATK FIX LCK`.
  Current phase letter: `B`. E.g. `B-IMP-genomewide`.
- **Directories/files:** `<phase>_<TARGET>[_<variant>]`, never two differing only
  in case. Put the varying parameter — above all `tau` and the arm — IN the name.
  No brace/glob shorthand in docs (write `tau_0p10/`, `tau_0p25/`).
- Deprecated things are renamed `DEPRECATED_*` or deleted, recorded in `docs/data.md`.

## THE IN-PLACE CORRECTION RULE

`docs/memory.md` is permanent. **Never delete or overwrite a historical entry.**
When something in it is proven wrong: prepend `[INCORRECT] - ` to the original
line, preserving its text verbatim, then insert directly below it
`[CORRECTION - YYYY-MM-DD]: ` with the new finding. Preserving the wrong version
is what makes the reasoning legible later; grep `[INCORRECT]` to list everything
this project has been wrong about. **Corrections do not propagate themselves** —
after writing one, hunt its stale copies by grepping the distinctive *number or
phrase*, not the topic, across `docs/`, `src/`, `tests/`, `paper/`, `CLAUDE.md`
and the assistant memory directory.

## YOU MUST EXECUTE: The Wrap-Up Protocol — at every completed unit of work

1. **Archive to `docs/memory.md`.** Hypothesis, method, provenance, result.
2. **Compress to the ledger.** One row in `docs/plan.md`.
3. **Document fixes.** `[Symptom]` -> `[Proven fix]` into `docs/bugs.md`.
4. **Define new terms.** Full entry in `docs/terms.md`.
5. **Register artifacts.** A `docs/data.md` row for every new output.
6. **Reset the board.** Rewrite `docs/plan.md` Current State; bump `Last updated`.
7. **Verify.** `python -m pytest tests/test_docs_contract.py` must pass.
