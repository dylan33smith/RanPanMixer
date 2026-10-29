"""ARM 2's per-target CEILING mode: the gates that keep it from renumbering a seed.

⚠⚠ WHY THESE EXIST. `--ceiling-grid` gives every target its OWN `eta_tau` at a shared
`acc_ceiling` column, because `c` = kl/eta^2 spreads 30.9274% across the ten arm 2
targets. The seed index of a cell can therefore no longer be the eta's POSITION in
`ETA_GRID` -- there are now 10 x 7 distinct etas -- so it became the CEILING's position
under a THIRD offset. `docs/bugs.md` records five occurrences of this project's
common-random-numbers defect, every one of them a seed index that moved. These assert
the namespaces cannot collide and that the pre-existing ones did not shift.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))


@pytest.fixture(scope="module")
def m():
    return importlib.import_module("arm2_attack")


def test_the_three_seed_namespaces_cannot_overlap(m):
    """⚠⚠ THE GATE THIS MODE EXISTS UNDER. A `tau` cell's index is its place in
    SEED_GRID, an `eta` cell's is ETA_SEED_OFFSET + its place in ETA_GRID, and a
    ceiling cell's is CEIL_SEED_OFFSET + its place in REPORTED_CEILINGS. If any two
    ranges met, two different cells would draw the SAME release stream and the table
    would silently report one experiment twice."""
    tau_hi = len(m.SEED_GRID) - 1
    eta_lo, eta_hi = m.ETA_SEED_OFFSET, m.ETA_SEED_OFFSET + len(m.ETA_GRID) - 1
    cl_lo, cl_hi = m.CEIL_SEED_OFFSET, m.CEIL_SEED_OFFSET + len(m.REPORTED_CEILINGS) - 1
    assert tau_hi < eta_lo, f"tau indices reach {tau_hi}, eta starts at {eta_lo}"
    assert eta_hi < cl_lo, f"eta indices reach {eta_hi}, ceiling starts at {cl_lo}"
    assert cl_lo > 0 and cl_hi > cl_lo


def test_the_pre_existing_eta_indices_did_not_move(m):
    """Every cell recorded on 2026-10-25 and 2026-10-28 keyed on its eta's position in
    ETA_GRID. Adding a mode must not renumber any of them, or those runs stop being
    reproducible from their key."""
    assert m.ETA_GRID[:6] == (0.0, 34.33, 68.66, 173.3, 387.39, 774.78)
    assert m.ETA_GRID[6:] == (75.7, 151.4, 227.1, 302.8, 378.5, 454.19)
    assert m.ETA_SEED_OFFSET == 1000
    assert m.SEED_BASE == 50_000


def test_reported_ceilings_match_the_grid_builder(m):
    """The axis is defined in `arm1_chr21_kl` and consumed here; two copies that drift
    would label a column with a ceiling nobody measured."""
    a1 = importlib.import_module("arm1_chr21_kl")
    assert tuple(m.REPORTED_CEILINGS) == tuple(a1.CEILINGS)


def _verified(tmp_path: Path, etas: dict, *, arm="arm2", scope="chr21") -> Path:
    per = {t: {"target": t, "first": {}, "refined":
               {lab: {"eta": e, "kl_exact": 1.0, "acc_ceiling": float(lab),
                      "error": 0.0} for lab, e in labs.items()}}
           for t, labs in etas.items()}
    p = tmp_path / "v.json"
    p.write_text(json.dumps({"arm": arm, "scope": scope, "k_states": 88,
                             "weights": "uniform", "worst_abs_error": 3e-8,
                             "per_target": per}))
    return p


ABOVE = ("0.55", "0.60", "0.65", "0.70", "0.75", "0.80")


def test_ceiling_grid_supplies_eta_zero_at_the_floor(m, tmp_path):
    """0.50 is `eta` = 0 EXACTLY -- the untilted control -- and `--verify` records only
    labels above 0.5, so the loader must supply it rather than leave a hole."""
    g = m.ceiling_grid(_verified(tmp_path, {
        "HG00096": {l: 100.0 * (i + 1) for i, l in enumerate(ABOVE)}}), ["HG00096"])
    assert g["HG00096"]["0.50"] == 0.0
    assert sorted(g["HG00096"], key=float) == ["0.50", *ABOVE]


def test_each_target_keeps_its_own_eta(m, tmp_path):
    """The whole point: a shared column, a per-target eta."""
    etas = {"HG00096": {l: 100.0 * (i + 1) for i, l in enumerate(ABOVE)},
            "HG00097": {l: 130.0 * (i + 1) for i, l in enumerate(ABOVE)}}
    g = m.ceiling_grid(_verified(tmp_path, etas), ["HG00096", "HG00097"])
    assert g["HG00096"]["0.80"] == 600.0
    assert g["HG00097"]["0.80"] == 780.0
    assert g["HG00096"]["0.80"] != g["HG00097"]["0.80"]


def test_a_fit_json_is_refused(m, tmp_path):
    """⚠ A grid quoted from a FIT is a derived number wearing a measured number's
    clothes, and in this mode the label IS the bound the run is checked against."""
    p = tmp_path / "fit.json"
    p.write_text(json.dumps({"arm": "arm2", "scope": "chr21", "weights": "uniform",
                             "per_target": {"HG00096": {"c": 1e-6, "grid": {}}}}))
    with pytest.raises(SystemExit, match="FIT and not a VERIFIED grid"):
        m.ceiling_grid(p, ["HG00096"])


def test_a_target_with_no_verified_ceiling_is_refused(m, tmp_path):
    """Running it under another target's eta is the bug this mode exists to remove."""
    g = _verified(tmp_path, {"HG00096": {l: 1.0 for l in ABOVE}})
    with pytest.raises(SystemExit, match="no verified ceiling"):
        m.ceiling_grid(g, ["HG00096", "HG00097"])


def test_the_wrong_arm_or_scope_is_refused(m, tmp_path):
    g = _verified(tmp_path, {"HG00096": {l: 1.0 for l in ABOVE}}, arm="arm1")
    with pytest.raises(SystemExit, match="not arm2"):
        m.ceiling_grid(g, ["HG00096"])


def test_a_partial_axis_is_refused(m, tmp_path):
    """Six columns where seven are reported would table a grid that stops early."""
    g = _verified(tmp_path, {"HG00096": {l: 1.0 for l in ABOVE[:-1]}})
    with pytest.raises(SystemExit, match="not the reported axis"):
        m.ceiling_grid(g, ["HG00096"])


# ---------------------------------------------------------------------------
# THE SEED-UNIQUENESS GATE. This is the gate occurrence SIX got past.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def m34():
    return importlib.import_module("arm2_attack_a3a4")


def _all_seeds(fn, n_targets: int, point_indices, n_draws: int) -> list:
    return [fn(ti, pi, k)
            for ti in range(n_targets)
            for pi in point_indices
            for k in range(n_draws)]


@pytest.mark.parametrize("mode", ["tau", "eta", "ceiling"])
def test_no_two_cells_share_a_release_stream(m, m34, mode):
    """⚠⚠ THE GATE OCCURRENCE SIX GOT PAST. `arm2_attack_a3a4.py` keyed its seed as
    `i * 10_000 + index * 1_000 + k`, which aliased once ETA_GRID reached 12 entries:
    (target i, eta 378.5, draw k) landed on (target i+1, eta 0.0, draw k). Measured on
    the shipped configuration: 450 of 1,750 cells, 25.71%, sharing a stream.

    Run over the REAL key space of each mode, for BOTH drivers, against the REAL
    `_cell_seed` -- not a copy of it."""
    if mode == "tau":
        idx = [m.SEED_GRID.index(t) for t in m.TAUS]
    elif mode == "eta":
        idx = [m.ETA_SEED_OFFSET + m.ETA_GRID.index(e) for e in m.ETA_GRID]
    else:
        idx = [m.CEIL_SEED_OFFSET + i for i in range(len(m.REPORTED_CEILINGS))]
    for name, fn in (("arm2_attack", m._cell_seed), ("arm2_attack_a3a4", m34._cell_seed)):
        seeds = _all_seeds(fn, 10, idx, 25)
        assert len(set(seeds)) == len(seeds), (
            f"{name} in {mode} mode: {len(seeds) - len(set(seeds))} of {len(seeds)} "
            f"cells share a release stream with another cell")


def test_the_two_drivers_do_not_share_streams(m, m34):
    """A1 and A3/A4 use the same SEED_BASE and the same (target, point, draw) key. If
    their streams coincided the two attacks would be scoring the SAME releases while
    reporting as independent evidence."""
    idx = [m.CEIL_SEED_OFFSET + i for i in range(len(m.REPORTED_CEILINGS))]
    a = set(_all_seeds(m._cell_seed, 10, idx, 25))
    b = set(_all_seeds(m34._cell_seed, 10, idx, 25))
    assert not (a & b), f"{len(a & b)} seeds are shared between the two drivers"


def test_the_old_arithmetic_form_really_did_collide(m):
    """⚠ The gate above must be about a REAL failure, not a hypothetical one. This
    reconstructs the pre-fix arithmetic and asserts it collides on the configuration
    that shipped on 2026-10-28 -- so if anyone restores that form, the gate above is
    known to catch something that actually happened."""
    etas = [0.0, 75.7, 151.4, 227.1, 302.8, 378.5, 454.19]
    seeds = [m.SEED_BASE + i * 10_000
             + (m.ETA_SEED_OFFSET + m.ETA_GRID.index(e)) * 1_000 + k
             for i in range(10) for e in etas for k in range(25)]
    dup = len(seeds) - len(set(seeds))
    assert dup == 225, f"expected the recorded 225 collisions, got {dup}"
