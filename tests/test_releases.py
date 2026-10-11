"""Does the release writer refuse every wrong write?

⚠ Each case below injects ONE specific error. A test that only checks the happy path
would pass on a writer that is a bare `np.save`, which is what this module replaced.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import releases as R                                    # noqa: E402
import store                                            # noqa: E402
from mechanism import ffbs                              # noqa: E402

DIGEST = "a" * 64
INPUTS = {"haplotypes.npy": "b" * 64}


def mk(chrom="chr1", *, tau=0.10, seed=7, share=1.0, loo_id=None, target="HG00096", T=5):
    return ffbs.Release(
        target_id=target, tau=tau, eta_tau=0.100335, utility_id="indicator/inv_support",
        graph_id="hprc-v1.0", seed=seed,
        state_path=np.arange(T, dtype=np.int32),
        released_path=np.arange(T, dtype=np.int16) % 2,
        # ⚠ BOTH SCALE WITH THE SHARE. The first version of this fixture used a flat
        # log_z_p = 0.05, so 22 parts summed to 1.05 against an eta_tau of 0.100335 --
        # and `GenomeRelease` caught it, which is the type doing exactly its job. A
        # per-chromosome log_z_p that does NOT scale means the budget was spent 22
        # times over, and the aggregate is the only thing that can see it.
        log_z_p=0.05 * share, u_achieved=0.4 * share, u_ceiling=0.8155,
        ploidy=1, loo_id=loo_id, chrom=chrom, share=share)


# ---------------------------------------------------------------------------
# the round trip, and the control for everything below
# ---------------------------------------------------------------------------

def test_a_release_round_trips_with_verified_provenance(tmp_path: Path) -> None:
    rel = mk()
    d = R.save_release(rel, arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    assert (d / "released_path.npy").exists() and (d / "state_path.npy").exists()
    # ⚠ store.verify RAISES on a digest mismatch, so calling it is the check
    rec = R.load_release(d)
    assert rec["release_id"] == rel.release_id
    assert np.array_equal(np.load(d / "released_path.npy"), rel.released_path)
    assert np.array_equal(np.load(d / "state_path.npy"), rel.state_path)


def test_the_axis_digest_reaches_every_sidecar(tmp_path: Path) -> None:
    """⚠ The digest is what lets a reader refuse to join a release to a different
    cohort build. It is useless if it is accepted and then dropped."""
    d = R.save_release(mk(), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    for f in ("released_path.npy", "state_path.npy", "release.json"):
        prov = json.loads((d / (f + store.PROV_SUFFIX)).read_text())
        assert prov["axis_digest"] == DIGEST, f
        assert prov["inputs"] == INPUTS, f
        assert prov["git_rev"], f"{f} records no git_rev"


# ---------------------------------------------------------------------------
# the collision
# ---------------------------------------------------------------------------

def test_a_second_write_to_one_key_raises(tmp_path: Path) -> None:
    """⚠ `Release`'s own docstring demanded this and nothing enforced it. Overwriting
    is the worst failure available: the second write looks exactly like the first, and
    the guarantee is about the path that was RELEASED, not the one left on disk."""
    rel = mk()
    R.save_release(rel, arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    with pytest.raises(FileExistsError, match="already exists"):
        R.save_release(rel, arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)


def test_a_different_seed_is_a_different_release_not_a_collision(tmp_path: Path) -> None:
    R.save_release(mk(seed=1), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    R.save_release(mk(seed=2), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    assert R.existing_seeds("HG00096", 0.10, "arm2", tmp_path) == [1, 2]


# ---------------------------------------------------------------------------
# the arm, which is never inferred
# ---------------------------------------------------------------------------

def test_arm1_without_leave_one_out_is_refused(tmp_path: Path) -> None:
    """Arm 1's guarantee holds ONLY under leave-one-out (constraint 1): the target IS
    a cohort member, so without it the release was drawn against a panel containing
    its own target."""
    with pytest.raises(ValueError, match="no loo_id"):
        R.save_release(mk(loo_id=None), arm="arm1", axis_digest=DIGEST,
                       inputs=INPUTS, root=tmp_path)


def test_arm2_with_leave_one_out_is_refused(tmp_path: Path) -> None:
    """Arm 2's target is external by construction, so there is nothing to hold out —
    a loo_id means an arm 1 release mislabelled."""
    with pytest.raises(ValueError, match="mislabelled"):
        R.save_release(mk(loo_id="HG00096"), arm="arm2", axis_digest=DIGEST,
                       inputs=INPUTS, root=tmp_path)


def test_arm1_with_leave_one_out_is_accepted(tmp_path: Path) -> None:
    """⚠ THE CONTROL for the two above: without it, a writer that refused EVERY arm 1
    release would pass them both."""
    d = R.save_release(mk(loo_id="HG00096"), arm="arm1", axis_digest=DIGEST,
                       inputs=INPUTS, root=tmp_path)
    assert d.is_dir() and "arm1" in str(d)


def test_an_unknown_arm_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="arm must be one of"):
        R.save_release(mk(), arm="arm3", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)


# ---------------------------------------------------------------------------
# provenance that cannot be skipped
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kwargs,match", [
    ({"axis_digest": ""}, "axis_digest is required"),
    ({"inputs": {}}, "inputs is empty"),
])
def test_unprovenanced_writes_are_refused(tmp_path: Path, kwargs, match) -> None:
    call = {"arm": "arm2", "axis_digest": DIGEST, "inputs": INPUTS, "root": tmp_path}
    call.update(kwargs)
    with pytest.raises(ValueError, match=match):
        R.save_release(mk(), **call)


# ---------------------------------------------------------------------------
# the tau directory name
# ---------------------------------------------------------------------------

def test_tau_directory_names_follow_the_naming_rule() -> None:
    assert R.tau_dir(0.0) == "tau_0p00"
    assert R.tau_dir(0.05) == "tau_0p05"
    assert R.tau_dir(0.10) == "tau_0p10"
    assert R.tau_dir(0.25) == "tau_0p25"
    assert R.tau_dir(0.99) == "tau_0p99"


def test_a_tau_that_would_collide_in_its_directory_name_is_refused() -> None:
    """⚠ 0.995 and 0.99 both render as `tau_0p99` at two decimals, so the second
    release would land in the first one's directory — a silent collision between two
    different PRIVACY LEVELS, which no reader could spot afterwards."""
    with pytest.raises(ValueError, match="does not survive"):
        R.tau_dir(0.995)
    assert R.tau_dir(0.99) == "tau_0p99"       # and the grid value still works


# ---------------------------------------------------------------------------
# genome scope -- the gap the in-memory type cannot see
# ---------------------------------------------------------------------------

def test_one_chromosome_claiming_genome_scope_is_refused(tmp_path: Path) -> None:
    """⚠⚠ THE GAP THIS MODULE WAS WRITTEN FOR. A `GenomeRelease` holding ONE
    chromosome with share 1.0 passes every check it makes — shares sum to 1, no
    duplicates, both bounds hold — so nothing in memory can tell it from a release of
    a whole genome. `tau` bounds the distinguishability of the WHOLE released path."""
    gr = ffbs.GenomeRelease(parts=(mk("chr1", share=1.0),))
    gr.record()                                  # the in-memory type is content
    with pytest.raises(AssertionError, match="not the 22 autosomes"):
        R.save_genome_release(gr, arm="arm2", axis_digest=DIGEST, inputs=INPUTS,
                              root=tmp_path, scope="genome")


def test_the_same_parts_are_fine_at_chromosome_scope(tmp_path: Path) -> None:
    """⚠ THE CONTROL: the refusal above must be about the SCOPE CLAIM, not about
    single-chromosome releases, which are legitimate and are most of our evidence."""
    gr = ffbs.GenomeRelease(parts=(mk("chr1", share=1.0),))
    out = R.save_genome_release(gr, arm="arm2", axis_digest=DIGEST, inputs=INPUTS,
                                root=tmp_path, scope="chromosome")
    assert json.loads(out.read_text())["scope"] == "chromosome"


def test_all_22_autosomes_at_genome_scope_are_accepted(tmp_path: Path) -> None:
    parts = tuple(mk(c, share=1.0 / 22) for c in R.AUTOSOMES)
    gr = ffbs.GenomeRelease(parts=parts)
    out = R.save_genome_release(gr, arm="arm2", axis_digest=DIGEST, inputs=INPUTS,
                                root=tmp_path, scope="genome")
    rec = json.loads(out.read_text())
    assert rec["scope"] == "genome" and rec["n_chromosomes"] == 22
    assert len(rec["parts"]) == 22


def test_a_missing_autosome_at_genome_scope_names_it(tmp_path: Path) -> None:
    """21 of 22, with the shares renormalised so every in-memory check still passes —
    which is exactly how this would reach disk unnoticed."""
    chroms = [c for c in R.AUTOSOMES if c != "chr7"]
    parts = tuple(mk(c, share=1.0 / len(chroms)) for c in chroms)
    gr = ffbs.GenomeRelease(parts=parts)
    with pytest.raises(AssertionError, match="chr7"):
        R.save_genome_release(gr, arm="arm2", axis_digest=DIGEST, inputs=INPUTS,
                              root=tmp_path, scope="genome")


# ---------------------------------------------------------------------------
# standing constraint 5, made visible
# ---------------------------------------------------------------------------

def test_the_sibling_seed_count_is_recorded_on_disk(tmp_path: Path) -> None:
    """Constraint 5 permits ONE release per genome. More seeds is a declared
    evaluation device, and the fact must be legible on disk rather than inferred by
    counting directories later."""
    d1 = R.save_release(mk(seed=1), arm="arm2", axis_digest=DIGEST, inputs=INPUTS,
                        root=tmp_path)
    d2 = R.save_release(mk(seed=2), arm="arm2", axis_digest=DIGEST, inputs=INPUTS,
                        root=tmp_path)
    p1 = json.loads((d1 / ("release.json" + store.PROV_SUFFIX)).read_text())
    p2 = json.loads((d2 / ("release.json" + store.PROV_SUFFIX)).read_text())
    assert p1["params"]["is_the_only_seed_at_this_tau"] is True
    assert p2["params"]["is_the_only_seed_at_this_tau"] is False
    assert p2["params"]["sibling_seeds_before_this_write"] == [1]


def test_two_seeds_cannot_be_quoted_as_a_deployed_release(tmp_path: Path) -> None:
    R.save_release(mk(seed=1), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    R.assert_one_release_per_genome("HG00096", 0.10, "arm2", tmp_path)   # one: fine
    R.save_release(mk(seed=2), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    with pytest.raises(AssertionError, match="COMPOSE"):
        R.assert_one_release_per_genome("HG00096", 0.10, "arm2", tmp_path)


def test_releases_at_different_tau_are_not_a_composition_problem(tmp_path: Path) -> None:
    """⚠ THE CONTROL: a `tau` grid is a declared evaluation device, so two releases at
    DIFFERENT tau must not trip the constraint-5 check — only two at the SAME tau do."""
    R.save_release(mk(tau=0.10), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    R.save_release(mk(tau=0.25), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    R.assert_one_release_per_genome("HG00096", 0.10, "arm2", tmp_path)
    R.assert_one_release_per_genome("HG00096", 0.25, "arm2", tmp_path)


def test_a_tampered_release_fails_verification(tmp_path: Path) -> None:
    """⚠ The sidecar is only worth something if reading goes THROUGH it."""
    d = R.save_release(mk(), arm="arm2", axis_digest=DIGEST, inputs=INPUTS, root=tmp_path)
    (d / "release.json").write_text(json.dumps({"release_id": "forged"}) + "\n")
    with pytest.raises(Exception):
        R.load_release(d)
