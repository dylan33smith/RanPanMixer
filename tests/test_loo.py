"""Does the leave-one-out chain actually hold the target out, and re-apply decision K?

⚠⚠ THE ONE TEST THAT MATTERS HERE is `test_a_block_only_the_target_could_span_is_dropped`.
`docs/plan.md` recorded the risk before chunk 8 opened: removing both of the target's
haplotypes shrinks every donor-intersection, so a zero-genetic-distance block that some
donor spans at K = 88 may be spanned by NONE at K = 86, and such a block DEAD-ENDS the
forward pass. The failure mode is specifically a block the TARGET was the only donor able
to span, which is unreachable with real data in a test, so it is CONSTRUCTED here: four
samples, one zero-distance block, and only the held-out individual called across the whole
block. The full cohort keeps it; leave-one-out must drop it.

⚠ The second real gate is the hold-out itself. Releasing one person's obfuscated path
while labelling it another's is silent and unrecoverable, so the rows held out are checked
against the ids rather than inferred from adjacency, and `support_loo <= support` is
asserted because removing donors cannot raise a called count.

⚠ `loo_chain`'s `n_dropped != PLOIDY` assertion is deliberately NOT tested: `sample_index`
already rejects any sample without exactly two haplotypes, so that branch is unreachable
through the public path. It is defence in depth against a future caller, and a contrived
test for it would test the contrivance. Named here so its absence is a decision.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from cohort import genmap, loo  # noqa: E402


def _write(d: Path, hap: np.ndarray, cm: np.ndarray, pos: np.ndarray, ids: list[str],
           support: np.ndarray | None = None) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "haplotypes.npy", hap.astype(np.int8))
    np.save(d / "genetic_pos.npy", cm.astype(np.float64))
    np.save(d / "positions.npy", pos.astype(np.int64))
    np.save(d / "haplotype_ids.npy", np.array(ids, dtype="<U32"))
    if support is not None:
        np.save(d / "support.npy", support.astype(np.int64))
    return d


def _fixture():
    """4 samples, 5 positions, ONE zero-distance block only the target can span.

    cm = [0, 1, 2, 2, 3] puts a zero step between positions 2 and 3, so (2, 3) is a
    zero-distance block. Across that block: A and C are called at 2 only, B at 3 only,
    and the TARGET at both. So the block is traversable at K = 8 and untraversable at
    K = 6.
    """
    ids = ["A:0", "A:1", "B:0", "B:1", "C:0", "C:1", "T:0", "T:1"]
    hap = np.zeros((8, 5), dtype=np.int8)
    hap[2:4, 2] = -1            # B missing at position 2
    hap[0:2, 3] = -1            # A missing at position 3
    hap[4:6, 3] = -1            # C missing at position 3
    cm = np.array([0.0, 1.0, 2.0, 2.0, 3.0])
    pos = np.array([100, 200, 300, 400, 500])
    return ids, hap, cm, pos


def test_a_block_only_the_target_could_span_is_dropped(tmp_path):
    """THE GATE. Arm 1's chain is per-target BECAUSE of this case."""
    ids, hap, cm, pos = _fixture()
    d = _write(tmp_path / "chrT", hap, cm, pos, ids)

    # at the FULL cohort the target spans the block, so nothing is dropped
    called = hap != -1
    full = genmap.chain_mask(pos, cm, None, None, called.sum(axis=0))
    full, trav_full = genmap.drop_untraversable(full, cm, called, pos)
    assert trav_full["untraversable_blocks_dropped"] == 0, (
        "the fixture is wrong: the block must be TRAVERSABLE with the target present, "
        "or the test proves nothing about leave-one-out")
    assert full.sum() == 5

    # holding the target out must now drop it
    lc = loo.loo_chain(d, "T")
    assert lc.k_states == 6
    assert lc.traversal["untraversable_blocks_dropped"] == 1
    assert lc.traversal["positions_dropped"] == 2
    assert lc.t_chain == 3
    assert not lc.in_chain[2] and not lc.in_chain[3]
    # and the refined chain must have no dead end left -- the property the gate exists for
    assert not genmap.untraversable_blocks(cm[lc.in_chain], called[:, lc.in_chain][
        [i for i in range(8) if i not in (6, 7)]])


def test_holding_out_a_different_member_keeps_the_block(tmp_path):
    """CONTROL. The drop must be caused by the HOLD-OUT, not by running the code."""
    ids, hap, cm, pos = _fixture()
    d = _write(tmp_path / "chrT", hap, cm, pos, ids)
    lc = loo.loo_chain(d, "A")
    assert lc.k_states == 6
    assert lc.traversal["untraversable_blocks_dropped"] == 0, (
        "holding out A leaves the target's spanning haplotypes in the panel, so the "
        "block is still traversable and nothing may be dropped")
    assert lc.t_chain == 5


def test_exactly_two_rows_leave_and_they_are_the_targets(tmp_path):
    ids, hap, cm, pos = _fixture()
    d = _write(tmp_path / "chrT", hap, cm, pos, ids)
    lc = loo.loo_chain(d, "T")
    assert lc.rows_dropped == (6, 7)
    assert lc.k_states == len(ids) - 2
    # support must be recomputed over the REMAINING donors, never reused
    assert lc.support[3] == 2, "only B is called at position 3 once T is held out"


def test_support_exceeding_the_full_cohort_raises(tmp_path):
    """A stored support LOWER than the recomputed one means the wrong rows left."""
    ids, hap, cm, pos = _fixture()
    bogus = np.zeros(5, dtype=np.int64)        # claims nothing is called anywhere
    d = _write(tmp_path / "chrT", hap, cm, pos, ids, support=bogus)
    with pytest.raises(AssertionError, match="exceeds the full-cohort"):
        loo.loo_chain(d, "T")


def test_an_arm2_target_is_refused(tmp_path):
    """Arm 2's target is outside D by construction; holding it out is meaningless."""
    ids, hap, cm, pos = _fixture()
    d = _write(tmp_path / "chrT", hap, cm, pos, ids)
    with pytest.raises(KeyError, match="not a cohort member"):
        loo.loo_chain(d, "HG00096")


def test_sample_index_rejects_unpaired_and_unlabelled():
    with pytest.raises(ValueError, match="no ':' separator"):
        loo.sample_index(np.array(["A:0", "BROKEN"], dtype="<U32"))
    with pytest.raises(ValueError, match="exactly 2 haplotypes"):
        loo.sample_index(np.array(["A:0", "A:1", "B:0"], dtype="<U32"))


def test_sample_index_does_not_assume_adjacency():
    """Reordered rows must still pair correctly, or the wrong person is held out."""
    idx = loo.sample_index(np.array(["A:0", "B:0", "A:1", "B:1"], dtype="<U32"))
    assert idx == {"A": (0, 2), "B": (1, 3)}


def test_shares_need_all_22_autosomes():
    """A share computed over a subset silently makes `tau` per-subset (constraint 4)."""
    with pytest.raises(ValueError, match="all 22 autosomes"):
        loo.loo_shares({"chr1": 100, "chr2": 200})
    full = {c: 1000 + i for i, c in enumerate(loo.AUTOSOMES)}
    sh = loo.loo_shares(full)
    assert abs(sum(sh.values()) - 1.0) < 1e-12
    assert sh["chr1"] < sh["chr22"], "shares must track each chromosome's own T"


def test_the_chain_ignores_the_held_out_data(tmp_path):
    """THE GATE THE RENAME RESTS ON. The chain depends on WHICH rows leave, never on
    what is in them.

    `src/cohort/loo.py` takes an individual as a parameter, which
    `tests/test_preprocess.py:test_no_cohort_function_takes_a_target` rightly forbids for
    a `target`. The parameter is named `held_out` because this module only DELETES those
    two rows. Renaming alone would be gaming that gate, so this proves the property it
    stands for: garble the held-out individual's two rows COMPLETELY -- alleles and
    missingness alike -- and every output must come back bit-identical. If any of
    `in_chain`, `support` or `k_states` moved, the cohort stage would be reading the
    held-out individual, which is the standing-constraint-1 violation.
    """
    ids, hap, cm, pos = _fixture()
    a = loo.loo_chain(_write(tmp_path / "clean", hap, cm, pos, ids), "T")

    garbled = hap.copy()
    rng = np.random.default_rng(0)
    # rows 6 and 7 are T's; replace them with arbitrary content INCLUDING missingness
    garbled[6:8, :] = rng.choice(np.array([-1, 0, 1], dtype=np.int8), size=(2, 5))
    assert not np.array_equal(garbled[6:8], hap[6:8]), "the fixture did not change"
    b = loo.loo_chain(_write(tmp_path / "garbled", garbled, cm, pos, ids), "T")

    assert np.array_equal(a.in_chain, b.in_chain), "in_chain moved with the held-out data"
    assert np.array_equal(a.support, b.support), "support moved with the held-out data"
    assert a.k_states == b.k_states
    assert a.traversal["positions_dropped"] == b.traversal["positions_dropped"]
    assert a.rows_dropped == b.rows_dropped


def test_garbling_a_REMAINING_donor_does_change_the_chain(tmp_path):
    """CONTROL for the test above: the invariance must be specific to the held-out rows.

    If the chain were invariant to EVERY row the test above would pass trivially and
    prove nothing. Garbling a donor that stays must move something.
    """
    ids, hap, cm, pos = _fixture()
    a = loo.loo_chain(_write(tmp_path / "clean", hap, cm, pos, ids), "T")
    other = hap.copy()
    other[0:2, :] = -1              # A becomes uncalled everywhere; it is a DONOR here
    b = loo.loo_chain(_write(tmp_path / "other", other, cm, pos, ids), "T")
    assert not np.array_equal(a.support, b.support), (
        "garbling a remaining donor left support unchanged, so the invariance test "
        "above is vacuous")
