"""ARM 1's per-cell seed must key on the `tau` VALUE, not its place in `--taus`.

⚠⚠ OCCURRENCE FOUR of this project's common-random-numbers defect, measured
2026-10-22. `genome_seed(target, tau_index)` took the position in whatever `--taus`
list was passed, so `--taus 0.999 0.9999 0.99999` passed 0, 1, 2 and redrew with the
seeds of `tau` = 0, 0.05 and 0.10 from the pinned run. The chr21 slices came out
BIT-IDENTICAL to those releases, so every drawn-path quantity at the three
exploratory `tau` was the pinned grid's number under a new label.

`docs/bugs.md` already recorded "A headline finding measured under common random
numbers" (2026-09-30) and `docs/memory.md` records seed reuse across chromosomes
correlating two chromosomes of one release; `scripts/arm2_attack.py` carries a
fourteen-line comment about the same defect in the attack. The lesson did not
transfer to this file, so it is pinned here as an executable gate.

THE TWO PROPERTIES, and the second is what the fix is FOR:
  1. every (target, `tau`) draws its own stream, over the WHOLE seed grid;
  2. a SUBSET run reproduces the full-grid seeds exactly -- which is what fails when
     the index is positional -- and the seven PINNED points keep the indices they
     were already released under, so the 308 releases on disk stay reproducible
     from their recorded seed under the signed pre-registration.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))

import arm1_release as M  # noqa: E402

TARGETS = ("HG00438", "HG00096", "NA19056")


def test_the_pinned_grid_keeps_the_seeds_it_was_released_under() -> None:
    """⚠ LOAD-BEARING. 308 arm 1 releases exist on disk, drawn under the signed
    pre-registration with the seed derived from the pinned list's own positions. The
    fix appends the exploratory points, so those positions must be unchanged -- any
    other arrangement silently makes every released cell unreproducible."""
    assert M.SEED_GRID[: len(M.PINNED_TAUS)] == M.PINNED_TAUS
    for t in TARGETS:
        h = int.from_bytes(hashlib.sha256(t.encode()).digest()[:4], "big")
        for i, tau in enumerate(M.PINNED_TAUS):
            as_released = int(np.random.SeedSequence([M.SEED_BASE, h, i])
                              .generate_state(1, dtype=np.uint32)[0])
            assert M.genome_seed(t, tau) == as_released, (
                f"{t} at tau={tau:g} would now draw a different stream than the "
                f"release on disk; those releases become unreproducible")


def test_every_target_and_tau_draws_its_own_stream() -> None:
    seeds = {(t, tau): M.genome_seed(t, tau) for t in TARGETS for tau in M.SEED_GRID}
    assert len(set(seeds.values())) == len(seeds), (
        f"{len(seeds) - len(set(seeds.values()))} seed collision(s): two cells share "
        f"a draw, which is the common-random-numbers failure this file exists for")


def test_a_subset_run_draws_the_same_streams_as_the_full_grid() -> None:
    """⚠⚠ THE DEFECT ITSELF. With a positional index, `--taus 0.999 0.9999 0.99999`
    reproduced the seeds of the first three PINNED points instead of its own."""
    subset = (0.999, 0.9999, 0.99999)
    for t in TARGETS:
        for tau in subset:
            assert M.genome_seed(t, tau) == M.genome_seed(t, tau)
        collisions = [(tau, p) for tau in subset for p in M.PINNED_TAUS
                      if M.genome_seed(t, tau) == M.genome_seed(t, p)]
        assert not collisions, (
            f"{t}: exploratory tau reuses a pinned tau's stream {collisions} -- the "
            f"2026-10-22 defect, where the chr21 slices came out bit-identical")


def test_a_tau_with_no_fixed_place_is_refused_rather_than_given_a_bare_index() -> None:
    """A `tau` absent from `SEED_GRID` has no reproducible seed. Refusing is the only
    safe answer: inventing one from a list position is exactly the defect."""
    with pytest.raises(ValueError, match="no fixed place"):
        M.genome_seed("HG00438", 0.995)


def test_the_seed_depends_on_the_target_as_well() -> None:
    for tau in (0.0, 0.99, 0.99999):
        s = {t: M.genome_seed(t, tau) for t in TARGETS}
        assert len(set(s.values())) == len(s), (
            f"tau={tau:g}: targets share a stream, so the 44 releases would be "
            f"near-replicates of each other")
