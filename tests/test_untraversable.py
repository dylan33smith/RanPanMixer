"""Decision K: drop the zero-genetic-distance blocks no donor spans.

⚠ WHY THE RULE EXISTS. Where two chain positions have ZERO genetic distance the
transition is exactly the identity (`s_t` = 0 in `A_t = (1-s_t) I + (s_t/K) J`), so no
path can switch donors inside a RUN of such steps. The donors able to traverse a run are
the INTERSECTION of `allowed` over all of it. If that intersection is empty the forward
pass dead-ends -- measured on the rebuilt chr21, which raised "zero total mass" at
position 10,543,014.

⚠ THE FIXTURES BELOW ARE HAND-CHECKABLE ON PURPOSE. The real chromosomes have 52,870
such runs of which only 20 are untraversable, so a bug that drops too much or too little
would be invisible against real data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from cohort import genmap as G                            # noqa: E402


# ---------------------------------------------------------------------------
# the run structure
# ---------------------------------------------------------------------------

def test_a_single_zero_step_is_a_two_position_block() -> None:
    """⚠ The +1 on the end matters: a zero STEP at t spans positions t AND t+1, so a
    block of one step is two positions. Off by one here under-drops by a position and
    leaves the dead end in place."""
    cm = np.array([0.0, 1.0, 1.0, 2.0])      # one zero step, between index 1 and 2
    assert G.zero_distance_blocks(cm) == [(1, 2)]


def test_consecutive_zero_steps_merge_into_one_run() -> None:
    cm = np.array([0.0, 1.0, 1.0, 1.0, 1.0, 2.0])
    assert G.zero_distance_blocks(cm) == [(1, 4)]


def test_separated_runs_stay_separate() -> None:
    cm = np.array([0.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 4.0])
    assert G.zero_distance_blocks(cm) == [(1, 2), (4, 6)]


@pytest.mark.parametrize("cm", [np.array([]), np.array([1.0]), np.array([0.0, 1.0, 2.0])])
def test_no_zero_steps_means_no_blocks(cm) -> None:
    assert G.zero_distance_blocks(cm) == []


# ---------------------------------------------------------------------------
# traversability -- the intersection, not pairwise
# ---------------------------------------------------------------------------

def test_a_block_one_donor_spans_is_traversable() -> None:
    """⚠ THE CONTROL, and the important one: 52,850 of the real 52,870 runs are
    traversable, so a rule that dropped every zero-distance run would remove 357,918
    positions instead of 172,142 and would look like it was working."""
    cm = np.array([0.0, 1.0, 1.0, 1.0, 2.0])
    allowed = np.array([[1, 1, 1, 1, 1],       # donor 0 spans everything
                        [1, 1, 0, 1, 1]], bool)
    assert G.untraversable_blocks(cm, allowed) == []


def test_a_block_no_single_donor_spans_is_untraversable() -> None:
    """⚠ PAIRWISE CHECKING WOULD MISS THIS. Every adjacent pair shares a donor --
    (0,1) share donor 0, (1,2) share donor 1 -- but no donor spans all three, so no
    path traverses the run."""
    cm = np.array([0.0, 1.0, 1.0, 1.0, 2.0])
    #                     block = positions 1,2,3
    allowed = np.array([[1, 1, 0, 1, 1],       # donor A: at 1 and 3, not 2
                        [1, 1, 1, 0, 1],       # donor B: at 1 and 2, not 3
                        [1, 0, 1, 1, 1]], bool)  # donor C: at 2 and 3, not 1
    # ⚠ THREE donors are the minimum that makes this point. With two, any fixture
    # where no donor spans a 3-position run also breaks one of its pairs -- which is
    # what the first version of this test got wrong.
    for a, b in ((1, 2), (2, 3)):
        assert allowed[:, a:b + 1].all(axis=1).any(), \
            f"pair ({a},{b}) must share a donor for this test to mean anything"
    assert not allowed[:, 1:4].all(axis=1).any(), "no donor may span the run"
    assert G.untraversable_blocks(cm, allowed) == [(1, 3)]


def test_a_donor_dropping_out_three_positions_later_is_still_caught() -> None:
    cm = np.array([0.0] + [5.0] * 6 + [9.0])
    allowed = np.ones((3, 8), bool)
    allowed[0, 3] = False
    allowed[1, 5] = False
    allowed[2, 6] = False
    assert G.untraversable_blocks(cm, allowed) == [(1, 6)]


# ---------------------------------------------------------------------------
# the drop
# ---------------------------------------------------------------------------

def _fixture():
    #            0    1    2    3    4    5     <- site axis
    cm = np.array([0.0, 1.0, 1.0, 1.0, 2.0, 3.0])
    allowed = np.array([[1, 1, 1, 0, 1, 1],
                        [1, 1, 0, 1, 1, 1]], bool)
    pos = np.array([100, 200, 300, 400, 500, 600])
    return cm, allowed, pos


def test_the_whole_block_is_dropped_and_nothing_else() -> None:
    cm, allowed, pos = _fixture()
    mask = np.ones(6, bool)
    out, info = G.drop_untraversable(mask, cm, allowed, pos)
    assert out.tolist() == [True, False, False, False, True, True]
    assert info["untraversable_blocks_dropped"] == 1
    assert info["positions_dropped"] == 3


def test_the_block_span_is_recorded_in_base_pairs() -> None:
    """A count cannot be checked against an assembly; a bp span can."""
    cm, allowed, pos = _fixture()
    _, info = G.drop_untraversable(np.ones(6, bool), cm, allowed, pos)
    assert info["blocks"] == [{"pos_start": 200, "pos_end": 400, "n_positions": 3}]


def test_the_refined_chain_has_no_untraversable_block_left() -> None:
    """The one-pass argument, as a test: `cm` is non-decreasing and the block was
    maximal, so the flanking positions are strictly separated and their new adjacency
    carries positive distance."""
    cm, allowed, pos = _fixture()
    out, _ = G.drop_untraversable(np.ones(6, bool), cm, allowed, pos)
    assert G.untraversable_blocks(cm[out], allowed[:, out]) == []
    assert (np.diff(cm[out]) > 0).all(), "a new zero-distance step was created"


def test_an_already_excluded_position_is_not_counted_twice() -> None:
    """The rule operates on the CHAIN, so positions already excluded (support == 0)
    must not re-enter the block arithmetic."""
    cm, allowed, pos = _fixture()
    mask = np.ones(6, bool)
    mask[3] = False                      # the position that broke the intersection
    out, info = G.drop_untraversable(mask, cm, allowed, pos)
    # with position 3 gone, donor 0 spans what remains of the run -> nothing to drop
    assert info["untraversable_blocks_dropped"] == 0
    assert out.tolist() == mask.tolist()


def test_a_fully_traversable_chain_is_returned_unchanged() -> None:
    cm = np.array([0.0, 1.0, 2.0, 3.0])
    allowed = np.ones((2, 4), bool)
    mask = np.ones(4, bool)
    out, info = G.drop_untraversable(mask, cm, allowed, np.arange(4))
    assert out.tolist() == mask.tolist()
    assert info["positions_dropped"] == 0


def test_the_rule_reads_no_target() -> None:
    """⚠ STANDING CONSTRAINT 1, asserted from the signature. A target-dependent
    exclusion would void Theorem 1; this one takes only the map, the allowed matrix
    and the positions, all functions of D and the genetic map."""
    import inspect
    for fn in (G.zero_distance_blocks, G.untraversable_blocks, G.drop_untraversable):
        params = set(inspect.signature(fn).parameters)
        assert not (params & {"target", "target_id", "target_path", "p", "g"}), (
            f"{fn.__name__} takes a target-shaped argument: {sorted(params)}")
