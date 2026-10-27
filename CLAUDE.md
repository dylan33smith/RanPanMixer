# RanPanMixer

A randomized **path release mechanism** for pangenome graphs. Given a public graph
`G` built from a public cohort `D`, it releases a target's mapped path in sanitized
form such that any two input paths induce output distributions at total variation
distance at most a pre-specified `tau`.

The release is a mosaic of cohort haplotypes drawn from a tilted Li-Stephens
posterior. It never contains the target's own data. The target need not belong to
the cohort.

## Layout

```
src/        the mechanism: HMM, FFBS sampler, utility, calibration, cohort/target IO
scripts/    runnable experiments and scoring, one concern per file
configs/    pinned run configurations
tests/      the gates; see "Testing" below
docs/       terms.md only -- the metric glossary the calibration gate parses
```

## Stack

Python with NumPy/SciPy. Environment `ranpanmixer`, defined in `environment.yml`.
`pytest` is included. Check before assuming any library is installed.

## The two arms

The arm decides what is true of the target. **Never state a property of the target
without naming its arm.**

- **Arm 1** — the target is one of the cohort members, released under **per-target
  leave-one-out**: both of the target's haplotypes leave the donor panel. Every
  leave-one-out, reduced-`K` and cohort-distortion quantity belongs to this arm.
  Used for head-to-head comparison against existing tools.
- **Arm 2** — the target is outside `G`, `D`, every panel and the attack database,
  so the full cohort is the donor panel and nothing is held out. Nothing in the
  mechanism is a function of the target but the tilt. This is the general case the
  method is for.

## Invariants

These are properties of the mechanism, not style preferences. Breaking one does not
degrade a result; it removes the guarantee.

1. **The output support `Y_D` and the baseline `R_D` are built from `G` and `D`
   only, never from the target.** Target-dependent pruning of the state space voids
   the theorem. Arm 2 holds by construction; **Arm 1 holds only under
   leave-one-out**, and `G` itself still saw the target, so that arm's guarantee is
   conditional on the graph. State that; never claim it away.
2. **The utility must be globally bounded in [0,1].** If a raw utility has range
   `delta_u != 1`, the calibration is `eta_tau = arctanh(tau)/delta_u`. An
   unnormalized utility does not weaken the bound, it removes it.
3. **Sample; never argmax.** Viterbi/MAP decoding is not a draw from the release
   distribution and carries no guarantee. The canonical decoder is `ffbs`.
4. **A failed attack is not evidence the bound holds.** Measured attacker accuracy
   is a lower bound on distinguishability. Only the proof gives the upper bound.
5. **Never mix two privacy levels, two utility functions, or two cohorts in one
   comparison.** Utility is comparable only at matched privacy, privacy only at
   matched utility. Leave-one-out gives each Arm 1 target its own `beta_v`, so its
   releases carry different utility functions: pool measured attack outcomes, never
   utilities.
6. **The guarantee covers the released path only.** Not the raw target genome, the
   graph, or any other released statistic about the same person.
7. **A utility comparison across weightings uses `e_match_nonref`.** `e_q_u` is
   `beta`-weighted and therefore weighting-specific. Compare at matched
   `acc_ceiling`, never at matched `eta_tau`: one `eta_tau` is two different privacy
   levels at two different scopes.

## Reporting

Metrics are rows, configurations are columns. **Every number carries its arm, its
privacy level, its denominator and its `n`** — one missing stamp makes a figure
uninterpretable rather than imprecise. Row labels are exact `docs/terms.md`
identifiers. Search that file before naming a quantity; never invent a synonym for
something already defined.

## Testing

```
pytest tests/
```

Several tests are deliberate gates rather than coverage, and some were added after
a defect reached a result. In particular `tests/test_mechanism.py` parses the
calibration table out of `docs/terms.md` and checks the code against it, because a
hardcoded copy in the test is a third transcription rather than a check. Do not
replace a parsing gate with a literal.

A gate that cannot fail is decoration. When adding one, inject the error it exists
to catch and confirm it fires.

## Long runs

Anything over about a minute runs detached with a status sentinel:

```
tmux new-session -d -s <name> '<cmd> > logs/<name>.log 2>&1; echo $? > logs/<name>.status'
```

Poll the sentinel; never report a run finished without reading it. **Do not edit
source while a run is in flight** — each chromosome imports its modules in a fresh
subprocess, so a mid-run edit splits the build by process start time and the
failure appears as a partial result rather than an error.

Any fan-out must assert its own completeness and fail loudly. A script that filters
results cannot distinguish "found nothing" from "ran nothing".

## Naming

Directories and files are `<phase>_<TARGET>[_<variant>]`, with the varying
parameter — above all the privacy level and the arm — in the name. No two names
differing only in case. Write `tau_0p10/`, not brace or glob shorthand.
