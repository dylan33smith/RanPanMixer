"""Does every check in `scripts/verify_cohorts.py` actually FIRE?

⚠ WHY THIS TEST IS A MUTATION TEST AND NOT A "DOES IT PASS" TEST. The verifier's
whole job is to refuse a cohort tree that is subtly wrong. Pointing it at a GOOD
tree and asserting it passes proves only that it is quiet, and a checker that is
always quiet is indistinguishable from no checker. So each case below copies a real
built chromosome, injects ONE specific corruption, and asserts the verifier reports
it. A case that stops failing when its check is deleted is the bug this file exists
to prevent.

⚠ These tests need a real built chromosome and SKIP without one. That is a genuine
coverage hole in a fresh checkout, not a pass — the skip reason says so, rather than
the suite reporting green on checks that never ran.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

# chr19 is the smallest fully built autosome, so the copy is cheap; any
# declared_rate_fill chromosome would do.
SRC = REPO / "data" / "cohort" / "chr19"
CHROM = 19

pytestmark = pytest.mark.skipif(
    not (SRC / "cm_source.npy").exists(),
    reason=f"{SRC} is not built under declared_rate_fill, so the verifier's checks "
           f"cannot be exercised against real data in this checkout")


@pytest.fixture()
def tree(tmp_path: Path) -> Path:
    """A writable copy of one real chromosome's arrays, for mutating."""
    d = tmp_path / f"chr{CHROM}"
    d.mkdir()
    for f in SRC.iterdir():
        if f.suffix in (".npy", ".json") and ".prov" not in f.name:
            shutil.copy2(f, d / f.name)
    return tmp_path


def _check(root: Path) -> list[str]:
    import verify_cohorts as V
    return V.check(CHROM, root)


def _manifest(root: Path, **kw) -> None:
    f = root / f"chr{CHROM}" / "manifest.json"
    m = json.loads(f.read_text())
    m.update(kw)
    f.write_text(json.dumps(m))


def test_control_unmutated_tree_passes(tree: Path) -> None:
    """⚠ THE CONTROL. Without it a verifier that reports a problem on EVERY input
    would pass all ten mutation cases below and be worthless."""
    assert _check(tree) == []


def test_truncated_per_site_array_is_caught(tree: Path) -> None:
    f = tree / f"chr{CHROM}" / "support.npy"
    np.save(f, np.load(f)[:-1])
    bad = _check(tree)
    assert any("rows, haplotypes has" in b for b in bad), bad


def test_length_mismatch_returns_before_positional_checks(tree: Path) -> None:
    """⚠ Regression. The positional checks RAISE on mismatched shapes rather than
    reporting, so a length disagreement must stop the run. The first version of
    this verifier appended the length problem and then crashed on the next check,
    which loses the other 21 chromosomes."""
    f = tree / f"chr{CHROM}" / "support.npy"
    np.save(f, np.load(f)[:-1])
    bad = _check(tree)
    assert len(bad) == 1 and "rows, haplotypes has" in bad[0], bad


def test_genotype_support_disagreement_is_caught(tree: Path) -> None:
    """support(v) is derived from genotypes AND read from INFO/AN — two independent
    routes to one number, so they can disagree, which is what makes this a test."""
    f = tree / f"chr{CHROM}" / "haplotypes.npy"
    h = np.load(f)
    h[0, 5] = -1 if h[0, 5] != -1 else 0
    np.save(f, h)
    assert any("called-genotype count" in b for b in _check(tree)), _check(tree)


def test_an_info_disagreement_is_caught(tree: Path) -> None:
    f = tree / f"chr{CHROM}" / "an_info.npy"
    a = np.load(f)
    a[7] += 1
    np.save(f, a)
    assert any("INFO/AN" in b for b in _check(tree)), _check(tree)


def test_decreasing_positions_are_caught(tree: Path) -> None:
    """A DECREASE makes a genetic distance negative, which puts the Li--Stephens
    switch probability outside [0,1] and raises nothing."""
    f = tree / f"chr{CHROM}" / "positions.npy"
    p = np.load(f)
    p[10], p[11] = p[11], p[10]
    np.save(f, p)
    assert any("positions DECREASE" in b for b in _check(tree)), _check(tree)


def test_a_shared_position_is_allowed(tree: Path) -> None:
    """⚠ THE FALSIFYING CASE. The first version of this check demanded STRICTLY
    increasing positions and failed chr9, where a parent bubble (lv=1) and a 1 bp
    deletion nested inside it (lv=2) genuinely share a start coordinate. Nesting is
    routine -- 10-18% of records per autosome carry lv > 0 -- so a tie is legal
    `vg deconstruct` output and must NOT be reported as a defect.

    This test is the reason the check is "non-decreasing" rather than "increasing",
    and it fails if anyone tightens it back."""
    f = tree / f"chr{CHROM}" / "positions.npy"
    p = np.load(f)
    p[11] = p[10]                      # a tie, exactly as chr9 has at 37,752,316
    np.save(f, p)
    assert _check(tree) == [], (
        "a shared start coordinate was reported as a defect; two records at one "
        "position is legal nested-variant output, not corruption")


def test_nan_cm_inside_chain_is_caught(tree: Path) -> None:
    d = tree / f"chr{CHROM}"
    cm = np.load(d / "genetic_pos.npy")
    ic = np.load(d / "in_chain.npy").astype(bool)
    cm[np.flatnonzero(ic)[100]] = np.nan
    np.save(d / "genetic_pos.npy", cm)
    assert any("no cM" in b for b in _check(tree)), _check(tree)


def test_decreasing_cm_is_caught(tree: Path) -> None:
    """A negative genetic step puts the Li--Stephens switch probability outside
    [0,1], and nothing downstream would raise."""
    d = tree / f"chr{CHROM}"
    cm = np.load(d / "genetic_pos.npy")
    ic = np.load(d / "in_chain.npy").astype(bool)
    i = np.flatnonzero(ic)[200]
    cm[i] -= 1.0
    np.save(d / "genetic_pos.npy", cm)
    assert any("cM decreases" in b for b in _check(tree)), _check(tree)


def test_usable_position_dropped_from_chain_is_caught(tree: Path) -> None:
    """Decision C retired the coordinate cut, so `in_chain` must be exactly
    `support > 0`; any other exclusion is an accidentally reintroduced span."""
    d = tree / f"chr{CHROM}"
    ic = np.load(d / "in_chain.npy")
    ic[np.flatnonzero(ic)[300]] = 0
    np.save(d / "in_chain.npy", ic)
    assert any("not exactly (support > 0)" in b for b in _check(tree)), _check(tree)


def test_interior_declared_rate_position_is_caught(tree: Path) -> None:
    """`numpy.interp` cannot leave an interior gap, so an interior declared-rate
    position means the end-fill was applied where it should not have been."""
    d = tree / f"chr{CHROM}"
    src = np.load(d / "cm_source.npy")
    ic = np.flatnonzero(np.load(d / "in_chain.npy").astype(bool))
    src[ic[len(ic) // 2]] = 1
    np.save(d / "cm_source.npy", src)
    assert any("INSIDE the mapped" in b for b in _check(tree)), _check(tree)


def test_axis_inconsistent_manifest_is_caught(tree: Path) -> None:
    _manifest(tree, axis_consistent=False)
    assert any("axis_consistent=False" in b for b in _check(tree)), _check(tree)


def test_pre_decision_c_chain_rule_is_caught(tree: Path) -> None:
    _manifest(tree, chain_rule="coordinate_cut")
    assert any("coordinate_cut" in b for b in _check(tree)), _check(tree)


def test_missing_artifact_is_reported_not_raised(tree: Path) -> None:
    """⚠ A chromosome built before decision C has no `cm_source.npy`. That must be
    a FINDING; a traceback here stops every other chromosome being checked."""
    (tree / f"chr{CHROM}" / "cm_source.npy").unlink()
    _manifest(tree, chain_rule=None)
    bad = _check(tree)
    assert len(bad) == 1
    assert "missing cm_source.npy" in bad[0] and "PRE-DECISION-C" in bad[0], bad
