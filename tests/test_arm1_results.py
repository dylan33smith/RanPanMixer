"""Can the arm 1 scoring script table a release that is not one?

⚠⚠ WHY THIS EXISTS. `scripts/arm1_results.py` is the only thing between 308 release
records and a number in the paper. A scoring script that formats whatever it is handed
cannot tell a release from a mislabelled one, and every failure here is SILENT: a record
with `loo_id` missing is a release drawn against a panel containing its own target, and it
looks identical in a table. So the script verifies on read, and these tests inject each
defect and assert it is caught.

⚠ The completeness refusal is tested too. A partial run tabled as a finished one is the
failure mode that makes a fan-out dangerous: it cannot tell "wrote nothing" from "ran
nothing".
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

_spec = importlib.util.spec_from_file_location(
    "arm1_results", REPO / "scripts" / "arm1_results.py")
AR = importlib.util.module_from_spec(_spec)
sys.modules["arm1_results"] = AR
_spec.loader.exec_module(AR)

AUT = [f"chr{i}" for i in range(1, 23)]


def _rec(**over):
    d = {"target_id": "HG00438", "loo_id": "HG00438", "tau": 0.5,
         "eta_tau": 0.5493061443340548, "log_z_p": 0.4484, "arm": "arm1",
         "scope": "genome", "n_chromosomes": 22, "chromosomes": list(AUT),
         "donor_switches": 830258, "T": 22131331, "u_ceiling": 0.9727,
         "u_achieved": 0.8161, "seed": 1, "ploidy": 2}
    d.update(over)
    return d


def _write(root: Path, name: str, rec: dict) -> None:
    p = root / name
    p.parent.mkdir(parents=True, exist_ok=True)
    (p.parent / f"{name}_genome.json").write_text(json.dumps(rec))


def test_a_clean_record_passes(tmp_path):
    """CONTROL. Without this the tests below could pass by rejecting everything."""
    _write(tmp_path, "ok", _rec())
    got = AR.load_cells(tmp_path)
    assert not got["problems"], got["problems"]
    assert list(got["cells"]) == [("HG00438", 0.5)]


@pytest.mark.parametrize("over,needle", [
    ({"loo_id": None}, "must hold out its own target"),
    ({"loo_id": "HG00621"}, "must hold out its own target"),
    ({"arm": "arm2"}, "not arm1"),
    ({"scope": "chromosome"}, "not genome"),
    ({"n_chromosomes": 21}, "not 22"),
    ({"chromosomes": AUT[:21]}, "not the 22 autosomes"),
    ({"log_z_p": 0.9}, "outside [0, eta_tau"),
    ({"log_z_p": -0.1}, "outside [0, eta_tau"),
])
def test_each_defect_is_caught(tmp_path, over, needle):
    """THE GATES, one injected defect each. Every one of these is silent in a table."""
    _write(tmp_path, "bad", _rec(**over))
    got = AR.load_cells(tmp_path)
    assert any(needle in p for p in got["problems"]), (
        f"injecting {over} produced no matching problem; got {got['problems']}")


def test_a_second_record_for_one_cell_is_caught(tmp_path):
    """Two genome records for one (target, `tau`) means constraint 5 was violated:
    independent re-releases of one genome compose and degrade the guarantee."""
    _write(tmp_path / "a", "one", _rec(seed=1))
    _write(tmp_path / "b", "two", _rec(seed=2))
    got = AR.load_cells(tmp_path)
    assert any("one-release-per-genome" in p for p in got["problems"])


def test_an_incomplete_set_refuses_to_report(tmp_path, capsys):
    """A partial run must not be tabled as a finished one."""
    _write(tmp_path, "ok", _rec())
    with pytest.raises(SystemExit) as e:
        AR.main(["--root", str(tmp_path)])          # defaults to all 7 pinned tau
    assert "not reportable" in str(e.value).lower()
    out = capsys.readouterr().out
    assert "NOT REPORTABLE" in out and "missing" in out


def test_a_verification_failure_refuses_even_when_complete(tmp_path, capsys):
    """A complete set that does not verify must still refuse -- completeness is not
    correctness."""
    _write(tmp_path, "bad", _rec(loo_id=None))
    with pytest.raises(SystemExit):
        AR.main(["--root", str(tmp_path), "--taus", "0.5"])
    assert "must hold out its own target" in capsys.readouterr().out


def test_the_full_metric_set_appears_with_reasons(tmp_path, capsys):
    """The contract forbids an omitted row: `n/a` must carry a REASON."""
    _write(tmp_path, "ok", _rec())
    assert AR.main(["--root", str(tmp_path), "--taus", "0.5"]) == 0
    out = capsys.readouterr().out
    for label in ("`target_fidelity` *", "`log_z_p` *", "`k_states` *",
                  "`p_succ_bound` *", "`donor_switches`", "`u_ceiling`",
                  "`utility_retained`", "`utility_loss_panmixer` *",
                  "`attacker_accuracy` *", "`af_loss`"):
        assert label in out, f"the metric set is missing {label}"
    assert "NOT ADOPTED, decision H" in out, "af_loss must say WHY it is absent"
    assert "decision M" in out, "the K = 88 scoring decision must be stamped"
    # the three caveats that must travel with every arm 1 number
    assert "JOINT bound does not hold" in out
    assert "DIFFERENT utility" in out
    assert "CONDITIONAL" in out
