"""The acc line, and the ERROR BAR a measured accuracy must carry against it.

⚠⚠ WHAT WENT WRONG, and why these gates are shaped the way they are. In CEILING mode
`scripts/arm2_attack_a3a4.py` stored the `acc_ceiling` column labels in its `taus` slot
and then handed them to `p_succ_bound`, so the line it tested against was
(1 + acc_ceiling)/2 -- 0.7500 at the 0.50 column up to 0.9000 at 0.80 -- instead of the
label itself. The line sat 0.10 to 0.25 ABOVE the real one for every A3 and A4 cell in
`logs/B_acc_arm2/a34_ceil.md`, under which a cell can exceed the true ceiling by 0.2 and
still print as comfortably within. `scripts/arm2_attack.py` was never affected: it got
the ceiling branch and this driver did not.

⚠⚠⚠ AND THE GATE THAT WAS SUPPOSED TO CATCH IT COULD NOT. `docs/memory.md` and
`docs/plan.md` both record that "the bound check was shown to FIRE in ceiling mode by
injecting a 0.96 accuracy at the 0.80 column". **0.96 clears BOTH lines** -- the correct
0.80 and the wrong 0.90 -- so that injection fires identically either way and
distinguishes nothing. `test_the_injection_discriminates_the_two_lines` is the same idea
done so that it can fail: it injects BETWEEN them.

The second rule these gates carry (Dylan, 2026-10-29): wherever a measured accuracy is
shown against the acc line, it carries its error bar, and where a point sits above the
line the bar must be shown REACHING BELOW it. A bar that clears the line is not a
violation on its own -- with 84 cells a few are expected above by chance -- but it must
be stated rather than absorbed by a corrected threshold.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO))

A34 = REPO / "logs" / "B_acc_arm2" / "a34_ceil.json"
A1 = REPO / "logs" / "B_acc_arm2" / "a1_ceil.json"

# A3's corrected threshold over 21 cells, and A4's over 15,603 block-cells. Both are
# the values the shipped run printed, used here as fixed stand-ins.
Z_A3, Z_A4 = 2.82, 4.51


@pytest.fixture()
def m():
    import arm2_attack_a3a4
    return arm2_attack_a3a4


# ---------------------------------------------------------------------------
# THE TWO LINES, and an injection that can tell them apart
# ---------------------------------------------------------------------------

def test_p_succ_bound_of_a_ceiling_is_the_wrong_line() -> None:
    """The arithmetic that caused it, pinned with its magnitude. An `acc_ceiling` is a
    legal `tau` value, which is exactly why nothing raised."""
    from src.mechanism import calibration as C
    assert C.p_succ_bound(0.50) == pytest.approx(0.75)
    assert C.p_succ_bound(0.80) == pytest.approx(0.90)
    # the gap, in the direction that hides an exceedance
    for ceil in (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80):
        assert C.p_succ_bound(ceil) > ceil
    assert C.p_succ_bound(0.50) - 0.50 == pytest.approx(0.25)
    assert C.p_succ_bound(0.80) - 0.80 == pytest.approx(0.10)


def test_the_injection_discriminates_the_two_lines(m) -> None:
    """⚠⚠ THE MUTATION TEST THE RECORDED ONE WAS NOT.

    0.85 at the 0.80 column sits ABOVE the correct line (0.80) and BELOW the wrong one
    (0.90). So it fires against the right bound and is silent against the wrong bound:
    the outcome tells you which line the code is using. That is the whole point.
    """
    se = 0.01
    right, _z, _lo, _hi = m.classify_excess(0.85, se, bound=0.80, z_crit=Z_A3)
    wrong, _z, _lo, _hi = m.classify_excess(0.85, se, bound=0.90, z_crit=Z_A3)
    assert right == "viol", "an 0.85 against a 0.80 line must be a violation"
    assert wrong == "within", "an 0.85 against a 0.90 line must be silent"


def test_the_recorded_0p96_injection_could_not_discriminate(m) -> None:
    """⚠ The gate on record is DECORATION, and this is the proof. 0.96 exceeds 0.80 and
    0.90 alike, so it fires whichever line is in force and cannot establish that ceiling
    mode tests against the ceiling. Kept as a test so the claim in `docs/memory.md` is
    checkable rather than asserted."""
    se = 0.01
    a = m.classify_excess(0.96, se, bound=0.80, z_crit=Z_A3)[0]
    b = m.classify_excess(0.96, se, bound=0.90, z_crit=Z_A3)[0]
    assert a == b == "viol", "0.96 fires against both lines, so it discriminates nothing"


# ---------------------------------------------------------------------------
# THE FOUR VERDICTS, each reached by construction
# ---------------------------------------------------------------------------

def test_a_genuine_excess_is_a_violation(m) -> None:
    v, z, lo, _hi = m.classify_excess(0.95, 0.01, bound=0.80, z_crit=Z_A3)
    assert v == "viol" and z > Z_A3 and lo > 0.80


def test_a_bar_that_clears_the_line_is_called_out(m) -> None:
    """⚠ THE NEW VERDICT. 0.88 against 0.80 at SE 0.03 is z = +2.67, under A3's
    corrected 2.82 -- so the claimable test stays quiet -- yet its 95% lower limit is
    0.8212, which does NOT reach below the line. Before 2026-10-29 this case was filed
    under "reported, not claimed" beside cells whose bars did reach below."""
    v, z, lo, _hi = m.classify_excess(0.88, 0.03, bound=0.80, z_crit=Z_A3)
    assert v == "bar_clears"
    assert z < Z_A3, "must not be corrected-significant, or it is just a violation"
    assert lo > 0.80, "the bar must genuinely fail to reach below the line"


def test_an_excess_whose_bar_reaches_below_is_only_reported(m) -> None:
    v, z, lo, _hi = m.classify_excess(0.81, 0.03, bound=0.80, z_crit=Z_A3)
    assert v == "reported" and z < Z_A3 and lo < 0.80


def test_a_cell_at_or_under_the_line_is_silent(m) -> None:
    assert m.classify_excess(0.80, 0.02, bound=0.80, z_crit=Z_A3)[0] == "within"
    assert m.classify_excess(0.60, 0.02, bound=0.80, z_crit=Z_A3)[0] == "within"


def test_a_zero_se_cell_above_the_line_is_a_violation_not_a_shrug(m) -> None:
    """⚠⚠⚠ THIS TEST ASSERTED THE OPPOSITE UNTIL 2026-10-29, AND THE OPPOSITE WAS WRONG.

    The reasoning that produced the bug: "se = 0 means a CONSTANT ANSWER, so it is
    neither a violation nor a null." The first half is only half the story. se = 0
    requires TPR and TNR BOTH in {0, 1}, which is four cases, not one:
        tpr 1 / tnr 0  -> balanced 0.5  -- the constant answer, no information
        tpr 0 / tnr 1  -> balanced 0.5  -- likewise
        tpr 1 / tnr 1  -> balanced 1.0  -- PERFECT separation, all the information
        tpr 0 / tnr 0  -> balanced 0.0  -- perfect, with the labels swapped
    Every line in this project is >= 0.5, so a constant answer can NEVER be above one.
    The only way to reach a zero-SE cell ABOVE the line is perfect separation -- and
    the old branch filed exactly that as "no evidence either way" and exited 0 on a
    total break of the mechanism.
    """
    v, z, _lo, _hi = m.classify_excess(0.85, 0.0, bound=0.80, z_crit=Z_A3)
    assert v == "viol", "a zero-SE cell above the line is PERFECT separation"
    assert z == float("inf")


def test_a_constant_answer_can_never_be_above_the_line(m) -> None:
    """The other half, and the reason the case above is unambiguous. A constant answer
    scores exactly 0.5; the lowest line in use is 0.50. So it is always 'within'."""
    import numpy as np
    bal, se, tpr, tnr = m._acc(np.ones(250), np.ones(250))   # always answers "p"
    assert (bal, se, tpr, tnr) == (0.5, 0.0, 1.0, 0.0)
    assert m.classify_excess(bal, se, bound=0.50, z_crit=Z_A3)[0] == "within"


def test_the_two_zero_se_cases_do_not_share_a_label(m) -> None:
    """⚠ A perfect attacker printed as "a constant answer" would be the single most
    misleading cell this suite could emit.

    ⚠⚠ THIS GATE ASSERTED THE WRONG THING UNTIL 2026-10-30. It required the words
    "constant" and "perfect separation" INSIDE the cell string, i.e. a cell reading
    `0.5000 DEGENERATE (a constant answer)`. What shipped is a SHORT marker in the
    cell with the reading given once in the table's preamble, which is the better
    convention -- spelling it out inside every cell would widen an already
    nine-column table -- and `test_the_markers_are_defined_where_a_reader_will_see_them`
    gates the preamble half. The invariant is that the two cases are DISTINGUISHABLE
    and RESOLVABLE, never that the cell is self-contained.
    """
    deg, exact = m._fmt(0.5, 0.0), m._fmt(1.0, 0.0)
    assert deg != exact, "the two opposite zero-SE cases print identically"
    assert "DEGENERATE" in deg and "EXACT" not in deg
    assert "EXACT" in exact and "DEGENERATE" not in exact
    # Neither may carry a bar: there is no interval around a zero-variance estimate,
    # and `[0.5000, 0.5000]` would claim infinite precision for the one cell that
    # carries no information at all.
    assert "[" not in deg and "[" not in exact


# ---------------------------------------------------------------------------
# THE BAR ITSELF
# ---------------------------------------------------------------------------

def test_the_interval_floor_is_zero_not_one_half(m) -> None:
    """⚠⚠ LOAD-BEARING FOR THE WHOLE RULE. A cell at 0.5000 with SE 0.0224 has a true
    lower limit of 0.4561. Clipping the bar at 0.5 -- plausible, since an attacker
    "cannot do worse than chance" -- would make "the bar reaches below 0.50" impossible
    to show at the one column where the line IS 0.50, and the six exceedances in
    `a34_ceil.json` at the 0.50 and 0.55 columns would all become unqualifiable."""
    lo, hi = m._interval(0.5, 0.0223606797749979)
    assert lo == pytest.approx(0.4561, abs=1e-4)
    assert lo < 0.5 and hi > 0.5


def test_the_interval_is_clipped_to_the_unit_range(m) -> None:
    assert m._interval(0.99, 0.2)[1] == 1.0
    assert m._interval(0.01, 0.2)[0] == 0.0


def test_a_degenerate_cell_is_never_printed_as_a_tight_bar(m) -> None:
    """A zero-width bracket reads as infinite precision for a cell with no information,
    which is the opposite of what the bar is for."""
    s = m._fmt(0.5, 0.0)
    assert "[" not in s and "DEGENERATE" in s
    assert "[" in m._fmt(0.5, 0.02)


def test_the_bar_half_width_is_the_two_sided_95_percent_z(m) -> None:
    import math
    assert m.Z95 == pytest.approx(1.959963984540054, abs=1e-12)
    # P(|Z| > Z95) = 0.05
    back = 2.0 * (0.5 * (1.0 - math.erf(m.Z95 / math.sqrt(2.0))))
    assert back == pytest.approx(0.05, abs=1e-12)


# ---------------------------------------------------------------------------
# THE SHIPPED RECORD. A regression gate on the two runs the paper quotes.
# ---------------------------------------------------------------------------

def _a34_cells(d: dict):
    tab, taus = d["table"], d["taus"]
    for cond, per in tab["a3"].items():
        for t, c in per.items():
            yield f"A3 {cond} @ {t}", c["balanced"], c["se"], float(t), Z_A3
    for size, byc in tab["a4"].items():
        for cond, per in byc.items():
            for t, c in per.items():
                yield (f"A4 b{size} {cond} @ {t}", c["max_balanced"],
                       c["se_at_max"], float(t), Z_A4)


def _a1_cells(d: dict):
    assert d["knob"] == "acc_ceiling", f"not a ceiling run: knob={d['knob']}"
    for cond, per in d["accuracy"].items():
        for t, c in per.items():
            yield f"A1 {cond} @ {t}", c["balanced"], c["se_balanced"], float(t), 2.82


@pytest.mark.skipif(not A34.exists() or not A1.exists(),
                    reason="the ceiling-mode runs are not on this machine")
@pytest.mark.parametrize("path,walk", [(A34, _a34_cells), (A1, _a1_cells)])
def test_no_shipped_cell_violates_or_clears_the_correct_line(path, walk) -> None:
    """⚠ Against the CORRECT line, not the one the table printed. Every cell above the
    ceiling must have a bar reaching below it, and none may be corrected-significant.
    This is the statement the paper rests on; it is asserted here against the stored
    per-cell data so a re-table cannot quietly change it."""
    import arm2_attack_a3a4 as m
    d = json.loads(path.read_text())
    bad, above = [], []
    for name, bal, se, bound, zc in walk(d):
        v, z, lo, _hi = m.classify_excess(bal, se, bound, zc)
        # ⚠ WHITELIST, not blacklist. Listing the bad verdicts lets any NEW verdict
        # fall through both branches and pass a test whose docstring promises every
        # cell above the line carries a bar that reaches below it.
        if v == "within":
            continue
        if v == "reported":
            above.append(name)
        else:
            bad.append(f"{name}: {bal:.4f} vs {bound:.4f}, z={z:+.2f}, lo={lo:.4f} [{v}]")
    assert not bad, "cells above the line without a bar reaching below:\n  " + \
                    "\n  ".join(bad)
    # ⚠ and the gate must be looking at something: a run with nothing above the line
    # would pass this vacuously, so the known exceedances are pinned by count.
    assert above, "no cell sits above the line at all -- is the gate reading the data?"


@pytest.mark.skipif(not A34.exists(), reason="the A3/A4 ceiling run is not here")
def test_the_six_known_exceedances_are_exactly_these(m) -> None:
    """The record, pinned. Six cells, all at the 0.50 and 0.55 columns -- which are
    also the only two columns the experiment's 0.0556 minimum detectable advantage
    cannot resolve below."""
    d = json.loads(A34.read_text())
    over = sorted(name for name, bal, se, bound, zc in _a34_cells(d)
                  if m.classify_excess(bal, se, bound, zc)[0] == "reported")
    assert len(over) == 6, f"expected 6 cells above the line, got {len(over)}: {over}"
    assert all(("@ 0.5" in n) or ("@ 0.55" in n) for n in over), over


@pytest.mark.skipif(not A34.exists(), reason="the A3/A4 ceiling run is not here")
def test_the_wrong_line_hid_every_one_of_them(m) -> None:
    """⚠⚠ THE CONSEQUENCE, measured. Under `p_succ_bound(label)` not one of the six
    exceedances is above the line, which is why the shipped table could close with "no
    condition exceeded" and be both true and vacuous."""
    from src.mechanism import calibration as C
    d = json.loads(A34.read_text())
    hidden = [name for name, bal, se, bound, zc in _a34_cells(d)
              if m.classify_excess(bal, se, bound, zc)[0] != "within"
              and m.classify_excess(bal, se, C.p_succ_bound(bound), zc)[0] == "within"]
    assert len(hidden) == 6, f"expected the wrong line to hide 6, hid {len(hidden)}"


# ---------------------------------------------------------------------------
# END TO END. The pure gates above pass a bound in explicitly, so none of them
# would notice if the reporter went back to choosing the wrong line. These do.
# ---------------------------------------------------------------------------

import subprocess  # noqa: E402

SCRIPT = REPO / "scripts" / "arm2_attack_a3a4.py"
GRID = REPO / "logs" / "B_acc_arm2" / "arm2_chr21_kl_verified.json"
WRONG_LINE = "0.7500 | 0.7750 | 0.8000 | 0.8250 | 0.8500 | 0.8750 | 0.9000"


def _retable(path: Path, *extra: str) -> tuple[int, str]:
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--from-json", str(path), *extra],
        capture_output=True, text=True, cwd=REPO,
        env={"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin",
             "HOME": str(Path.home()), "MPLBACKEND": "Agg"})
    return r.returncode, r.stdout + r.stderr


@pytest.mark.skipif(not A34.exists() or not GRID.exists(),
                    reason="the A3/A4 ceiling run is not on this machine")
def test_ceiling_mode_prints_the_ceiling_as_the_bound() -> None:
    """⚠⚠ THE GATE A REVERT TRIPS. If the reporter goes back to `C.p_succ_bound(t)`,
    the bound row becomes 0.7500 .. 0.9000 and this fails on both assertions."""
    rc, out = _retable(A34, "--ceiling-grid", str(GRID))
    assert rc == 0, out[-3000:]
    assert "`acc_ceiling` [THE BOUND, measured per target] * | 0.5000 | 0.5500 | " \
           "0.6000 | 0.6500 | 0.7000 | 0.7500 | 0.8000 |" in out, \
           "the ceiling is not printed as the bound"
    assert WRONG_LINE not in out, \
        "a bound row carrying p_succ_bound(acc_ceiling) is back"
    assert "One-sided on the excess over `acc_ceiling`" in out
    # `p_succ_bound` must still appear, named as vacuous at the tilted columns --
    # but NOT at the 0.50 control column, where eta is 0 and it is exactly 0.5.
    row = next(ln for ln in out.splitlines() if ln.startswith("| `p_succ_bound`"))
    assert "0.5000 (exact, eta = 0)" in row, row
    assert "~1.0 (vacuous)" in row, row


@pytest.mark.skipif(not A34.exists() or not GRID.exists(),
                    reason="the A3/A4 ceiling run is not on this machine")
def test_an_injected_cell_between_the_two_lines_fails_the_run(tmp_path: Path) -> None:
    """⚠⚠ THE DISCRIMINATING INJECTION, end to end.

    0.85 at the 0.80 column, SE 0.01: z = +5.0 against the correct line, over A4's
    corrected 4.51, so it MUST fail the run. And 0.85 < p_succ_bound(0.80) = 0.90, so
    under the wrong line it would not even be above the bound -- the run would exit 0
    and say nothing. That asymmetry is what makes this a real gate rather than the
    0.96 injection on record, which clears both lines.
    """
    d = json.loads(A34.read_text())
    cell = d["table"]["a4"]["512"]["D2_nearest"]["0.8"]
    cell["max_balanced"], cell["se_at_max"] = 0.85, 0.01
    from src.mechanism import calibration as C
    assert 0.80 < 0.85 < C.p_succ_bound(0.80), "the injection must sit BETWEEN the lines"
    f = tmp_path / "injected.json"
    f.write_text(json.dumps(d))

    rc, out = _retable(f, "--ceiling-grid", str(GRID))
    assert rc == 1, f"a 5-SE excess over the ceiling did not fail the run\n{out[-3000:]}"
    assert "falsifies the IMPLEMENTATION" in out
    assert "A4 block 512 D2_nearest at acc_ceiling=0.8" in out


@pytest.mark.skipif(not A34.exists(), reason="the A3/A4 ceiling run is not here")
def test_a_json_with_no_recorded_knob_is_refused_not_guessed() -> None:
    """⚠⚠ The labels 0.50 .. 0.80 are legal in BOTH modes, so a re-table that guesses
    picks a bound line that differs by up to 0.25 without raising anything. Refusing is
    the fix. `a34_ceil.json` predates the `knob` field, so this runs on the real file."""
    assert "knob" not in json.loads(A34.read_text()), \
        "this fixture must be a pre-2026-10-29 JSON for the test to mean anything"
    rc, out = _retable(A34)
    assert rc != 0, "a knob-less JSON was re-reported without being told its mode"
    assert "records no `knob`" in out and "will not be guessed" in out


def test_a_ceiling_run_now_records_its_knob() -> None:
    """The other half: new runs must persist the mode, so the refusal above becomes
    unnecessary rather than permanent. Checked on the writer, not a fixture."""
    src = (REPO / "scripts" / "arm2_attack_a3a4.py").read_text()
    assert '"knob": ("acc_ceiling" if getattr(a, "ceiling_etas", None)' in src
    assert '"ceiling_etas": getattr(a, "ceiling_etas", None),' in src


def test_the_driver_can_fail_at_all() -> None:
    """⚠ `arm2_attack.py` has always returned `1 if violations else 0`; this driver
    returned 0 unconditionally, so a falsifying cell exited clean and every queue
    reading `$?` logged `fail=0`. Pinned on the source so it cannot regress silently."""
    src = (REPO / "scripts" / "arm2_attack_a3a4.py").read_text()
    assert src.count("return 1 if n_viol else 0") == 2, \
        "both the run path and the --from-json path must propagate a violation"


# ---------------------------------------------------------------------------
# THE SAME DISCRIMINATION, WITHOUT THE DATASET.
#
# ⚠⚠ Every gate above that can tell the two lines apart is skipif'd on files under
# logs/, and logs/ is gitignored (.gitignore:13). On a fresh clone or in CI they all
# skip, the suite goes green, and reverting bound_at to C.p_succ_bound(t) is invisible.
# These call `_report` directly with synthetic cells and a stubbed dataset, so the
# discriminating check runs on any machine.
# ---------------------------------------------------------------------------

CEILS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80)


def _synthetic(monkeypatch, m, a3_acc: dict, a4_acc: dict):
    """Drive _report with hand-built cells. -> captured stdout, violation count.

    `a3_acc` / `a4_acc` map a ceiling label to the balanced accuracy wanted for
    D2_nearest; everything else sits at chance.
    """
    import numpy as np

    class _Beta:
        beta = np.zeros(1000)

    monkeypatch.setattr(m, "load_state", lambda *a, **k: {"util": _Beta()})
    monkeypatch.setattr(m, "a3_sites", lambda st: {"n_bi": 302547})

    def trials_for(tau):
        """n draws whose tie-split balanced accuracy is exactly `acc`."""
        acc = a3_acc.get(tau, 0.5)
        n = 250
        k = int(round(acc * n))
        out = []
        for i in range(n):
            right = i < k
            cell = {"a3_under_p": 1.0 if right else -1.0,
                    "a3_under_q": -1.0 if right else 1.0}
            out.append({"target": f"T{i % 10}", "tau": tau,
                        "conditions": {c: dict(cell) if c == "D2_nearest"
                                       else {"a3_under_p": 0.0, "a3_under_q": 0.0}
                                       for c in m.CONDS}})
        return out

    res = [t for tau in CEILS for t in trials_for(tau)]
    prior = {"a4": {str(sz): {c: {tau: {
        "n_blocks": 10,
        "max_balanced": (a4_acc.get(tau, 0.5) if c == "D2_nearest" else 0.5),
        "argmax_block": 0, "se_at_max": 0.01, "mean_balanced": 0.5,
        "n_per_arm": 250, "block_balanced": [0.5] * 10}
        for tau in CEILS} for c in m.CONDS} for sz in m.BLOCK_SIZES}}
    ceilings = {f"T{i}": {f"{c:.2f}": (0.0 if c == 0.50 else 100.0 * c)
                          for c in CEILS} for i in range(10)}
    ctl = {"a3": {"rate": 1.0, "n_correct": 30, "n": 30},
           "a4": {"rate": 1.0, "n_correct": 30, "n": 30}}
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        _table, n_viol = m._report(res, ctl, tuple(CEILS), "synthetic", "rev",
                                   [f"T{i}" for i in range(10)], 25, pw=None,
                                   ceilings=ceilings, prior_table=prior)
    return buf.getvalue(), n_viol


def test_without_logs_the_ceiling_is_still_the_bound(monkeypatch, m) -> None:
    """The display half, dataset-free."""
    out, n_viol = _synthetic(monkeypatch, m, {}, {})
    assert "`acc_ceiling` [THE BOUND, measured per target] * | 0.5000 | 0.5500 | " \
           "0.6000 | 0.6500 | 0.7000 | 0.7500 | 0.8000 |" in out
    assert WRONG_LINE not in out, "the p_succ_bound-of-a-ceiling row is back"
    assert n_viol == 0


def test_without_logs_a_cell_between_the_lines_still_fails(monkeypatch, m) -> None:
    """⚠⚠ THE GATE THAT SURVIVES A FRESH CLONE. 0.85 at the 0.80 column is above the
    correct line and below p_succ_bound(0.80) = 0.90, so a revert to the wrong line
    turns this violation into silence."""
    from src.mechanism import calibration as C
    assert 0.80 < 0.85 < C.p_succ_bound(0.80)
    out, n_viol = _synthetic(monkeypatch, m, {}, {0.80: 0.85})
    assert n_viol >= 1, f"a 5-SE excess over the ceiling did not register\n{out[-2500:]}"
    assert "falsifies the IMPLEMENTATION" in out


def test_without_logs_the_control_column_bound_is_exact_not_vacuous(monkeypatch, m) -> None:
    """At the 0.50 column eta is 0 by construction, so tau is 0 and p_succ_bound is
    exactly 0.5 — a real bound, and the tightest in the table. It was printed as
    "~1.0 (vacuous)" alongside the genuinely vacuous columns."""
    out, _ = _synthetic(monkeypatch, m, {}, {})
    row = next(ln for ln in out.splitlines() if ln.startswith("| `p_succ_bound`"))
    assert "0.5000 (exact, eta = 0)" in row, row
    assert row.count("~1.0 (vacuous)") == len(CEILS) - 1, row


def test_the_corrected_thresholds_are_computed_not_assumed(monkeypatch, m) -> None:
    """⚠ Z_A3 and Z_A4 are LITERALS in this file, so every pure gate above would pass
    if the driver computed its thresholds wrongly — e.g. correcting A3 over A4's
    15,603 block-cells, which would turn an A3 violation into a `bar_clears` and exit
    0. This reads the thresholds the driver actually printed and checks them against
    the cell counts it actually examined."""
    import re
    out, _ = _synthetic(monkeypatch, m, {}, {})
    mt = re.search(r"A3 over (\d+) cells \(z_crit = ([\d.]+)\), A4 over ([\d,]+) "
                   r"block-cells \(z_crit = ([\d.]+)\)", out)
    assert mt, "the bound-check header did not print its thresholds"
    a3_cells, z3, a4_cells, z4 = (int(mt.group(1)), float(mt.group(2)),
                                  int(mt.group(3).replace(",", "")), float(mt.group(4)))
    assert a3_cells == len(m.CONDS) * len(CEILS) == 21
    assert a4_cells == len(m.BLOCK_SIZES) * len(m.CONDS) * len(CEILS) * 10
    assert z3 == pytest.approx(m._z_for_one_sided(0.05 / a3_cells), abs=5e-3)
    assert z4 == pytest.approx(m._z_for_one_sided(0.05 / a4_cells), abs=5e-3)
    assert z4 > z3, "more cells examined must mean a stricter per-cell threshold"


# ---------------------------------------------------------------------------
# THE FIGURE'S OWN GATE. docs/data.md records it as mutation-tested; until now
# nothing in tests/ touched scripts/arm2_figures.py at all.
# ---------------------------------------------------------------------------

def _series(points):
    return {"A1 — optimal distinguisher": {"D1_extremal": points,
                                           "D2_nearest": [(c, 0.5, 0.02) for c in CEILS],
                                           "D3_random": [(c, 0.5, 0.02) for c in CEILS]}}


def test_the_figure_refuses_to_draw_an_unqualified_exceedance() -> None:
    """⚠⚠ A point above the diagonal whose 95% bar does not reach below it would be
    drawn as an unqualified exceedance. That must stop the run, not reach a PDF."""
    import arm2_figures as F
    pts = [(c, 0.5, 0.02) for c in CEILS[:-1]] + [(0.80, 0.88, 0.03)]
    with pytest.raises(SystemExit) as e:
        F._gate_bars_reach_below(_series(pts), 2.82, 4.51)
    assert "NOT reaching below" in str(e.value)


def test_the_figure_allows_an_exceedance_whose_bar_reaches_below() -> None:
    """The other direction: the gate must not simply refuse every point above the
    line, or it would be unusable on the real data, where seven sit above it."""
    import arm2_figures as F
    pts = [(c, 0.5, 0.02) for c in CEILS[:-1]] + [(0.80, 0.81, 0.03)]
    above = F._gate_bars_reach_below(_series(pts), 2.82, 4.51)
    assert len(above) == 1 and "0.80" in above[0]


def test_the_figure_gate_rejects_a_perfect_attacker() -> None:
    """A zero-SE cell above the line is perfect separation. It has no bar at all, so
    it can never be drawn as qualified."""
    import arm2_figures as F
    pts = [(c, 0.5, 0.02) for c in CEILS[:-1]] + [(0.80, 1.0, 0.0)]
    with pytest.raises(SystemExit):
        F._gate_bars_reach_below(_series(pts), 2.82, 4.51)


def test_the_markers_are_defined_where_a_reader_will_see_them(monkeypatch, m) -> None:
    """The other half of the short-marker convention, and the half that can rot.

    A short marker is only legitimate if the emitted table says what it means. This
    drives `_report` dataset-free and asserts the preamble resolves BOTH markers, so
    dropping that line fails here rather than silently shipping a table whose
    `EXACT` -- a perfect attacker -- is undefined anywhere in the file.
    """
    out, _ = _synthetic(monkeypatch, m, {}, {})
    low = out.lower()
    assert "degenerate" in low and "exact" in low, "the markers are not named at all"
    assert "constant answer" in low, (
        "the emitted table uses DEGENERATE without defining it as a constant answer")
    assert "perfect separation" in low, (
        "the emitted table uses EXACT without defining it as perfect separation")
