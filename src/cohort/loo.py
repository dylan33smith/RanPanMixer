"""The per-target leave-one-out cohort, `loo_cohort` — ARM 1's core construction.

⚠⚠ ARM 1 ONLY. Arm 2's target is outside `G`, `D`, every panel and the attack database,
so nothing is held out and K = 88. Every quantity here belongs to arm 1 and must be
stamped with its arm and its held-out individual. `docs/terms.md:loo_cohort` is the
authority; decision E (2026-10-08) settled that BOTH of the target's haplotypes leave the
donor panel, that it is applied PER TARGET, and that it is done by MASKING THE TWO ROWS
IN MEMORY rather than materialising 44 cohort trees (~160-380 GB cheaper, identical
artifacts).

⚠⚠ WHY THE CHAIN ITSELF IS PER-TARGET, which is the part that bites. Removing two
haplotypes lowers `support(v)`, so positions where only the target was called leave the
chain. It also shrinks every donor-intersection, so a zero-genetic-distance block that
SOME donor could span with 88 may be spanned by NONE with 86 — and such a block dead-ends
the forward pass. That is exactly the failure decision K exists to remove, so decision K
must be RE-APPLIED per target. Reusing arm 2's `in_chain` here would reintroduce the dead
end: measured on the rebuilt chr21, the un-refined chain raised "zero total mass" at
position 10,543,014.

⚠ WHAT IS *NOT* A SUBSET. It is tempting to assert `in_chain_loo` is contained in arm 2's
`in_chain`, and it is NOT provable. Dropping positions SHORTENS the zero-distance runs
that remain, and a shorter run needs a donor to span fewer positions, so a block arm 2
found untraversable can become traversable under leave-one-out. `support_loo <= support`
IS provable and is asserted; containment of the chain is MEASURED and reported, never
assumed.

⚠ Standing constraint 1 still holds and is worth saying precisely, because this module
looks like it breaks it. The held-out individual IS the target, so `D` here is chosen with
reference to the target. That is the DEFINITION of leave-one-out and it is what makes arm
1's guarantee conditional rather than void: `R_D` is built from `D = HPRC \\ {target}`,
which contains no information from the target at all, so the baseline the bound is stated
against is still target-independent GIVEN the hold-out. What remains conditional is `G`
itself, which was built from all HPRC assemblies including the target's. State that; never
claim it away.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cohort import genmap

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
PLOIDY = 2


# ⚠⚠ WHY THE PARAMETER IS `held_out` AND NOT `target`.
# `tests/test_preprocess.py:test_no_cohort_function_takes_a_target` forbids any function
# under `src/cohort/` from taking a target-shaped parameter, and it is RIGHT to: standing
# constraint 1 says `Y_D` and `R_D` are built from `G` and `D` only. This module takes an
# individual, but it only ever DELETES that individual's two rows -- it never reads their
# alleles, their path, or anything under `data/targets/`. `held_out` is the accurate name
# for that role, and renaming alone would be gaming the gate, so it is paired with a
# STRICTLY STRONGER check: `tests/test_loo.py:test_the_chain_ignores_the_held_out_data`
# garbles the held-out individual's two rows and asserts the chain comes back
# bit-identical. That proves the dependence is on WHICH rows leave, never on what is in
# them -- which is the property constraint 1 actually cares about.
@dataclass(frozen=True)
class LooChain:
    """One held-out individual's chain on one chromosome. Everything is arm 1."""

    held_out: str
    chrom: str
    rows_dropped: tuple[int, ...]
    k_states: int
    support: np.ndarray           # (n_sites,) recomputed over the REMAINING donors
    in_chain: np.ndarray          # (n_sites,) bool, decision K re-applied
    traversal: dict
    n_sites: int

    @property
    def t_chain(self) -> int:
        return int(self.in_chain.sum())


def sample_index(haplotype_ids: np.ndarray) -> dict[str, tuple[int, ...]]:
    """`{sample: (row, row)}` from the `SAMPLE:HAP` ids, with the pairing CHECKED.

    ⚠ The pairing is not assumed from adjacency. A cohort whose rows were ever reordered
    would silently hold out the wrong individual's haplotypes, which is the one error in
    this module that would produce a plausible release of the wrong person.
    """
    idx: dict[str, list[int]] = {}
    for i, raw in enumerate(haplotype_ids):
        s = str(raw)
        if ":" not in s:
            raise ValueError(
                f"haplotype id {s!r} has no ':' separator, so the sample it belongs to "
                f"cannot be read. Expected 'SAMPLE:HAP'.")
        idx.setdefault(s.rsplit(":", 1)[0], []).append(i)
    bad = {s: v for s, v in idx.items() if len(v) != PLOIDY}
    if bad:
        raise ValueError(
            f"{len(bad)} sample(s) do not have exactly {PLOIDY} haplotypes: "
            f"{dict(list(bad.items())[:3])}. A diploid hold-out needs both rows.")
    return {s: tuple(v) for s, v in idx.items()}


def loo_chain(chrom_dir: str | Path, held_out: str, *,
              haplotypes: np.ndarray | None = None) -> LooChain:
    """Build `held_out`'s leave-one-out chain on one chromosome. Decision K re-applied.

    `haplotypes` may be passed in so a caller sweeping 44 targets loads each
    chromosome's array ONCE instead of 44 times; it must be the full K-row array.
    """
    d = Path(chrom_dir)
    ids = np.load(d / "haplotype_ids.npy", allow_pickle=True)
    hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r")) \
        if haplotypes is None else haplotypes
    if hap_all.shape[0] != len(ids):
        raise ValueError(
            f"{hap_all.shape[0]} haplotype rows against {len(ids)} ids in "
            f"{d.name}; the panel and its labels disagree.")

    index = sample_index(ids)
    if held_out not in index:
        raise KeyError(
            f"{held_out!r} is not a cohort member of {d.name}; there are "
            f"{len(index)} members. ⚠ An ARM 2 target must NOT be passed here — it is "
            f"outside D by construction and nothing is held out for it.")
    rows = index[held_out]

    keep = np.ones(hap_all.shape[0], dtype=bool)
    keep[list(rows)] = False
    # ⚠ assert the hold-out, do not trust the arithmetic: releasing the wrong person's
    # obfuscated path while labelling it the target is silent and unrecoverable.
    n_dropped = int((~keep).sum())
    if n_dropped != PLOIDY:
        raise AssertionError(
            f"held out {n_dropped} rows for {held_out}, not {PLOIDY}. K = 86 requires "
            f"BOTH haplotypes to leave the donor panel (decision E).")
    hap = hap_all[keep]

    support = (hap != -1).sum(axis=0).astype(np.int64)
    ref = d / "support.npy"
    if ref.exists():
        # provable and worth checking: removing donors cannot RAISE the called count
        full = np.load(ref).astype(np.int64)
        if (support > full).any():
            j = int(np.flatnonzero(support > full)[0])
            raise AssertionError(
                f"{d.name} position {j}: leave-one-out support {int(support[j])} "
                f"exceeds the full-cohort {int(full[j])}. Removing two donors cannot "
                f"raise a called count, so the rows held out are not this target's.")

    positions = np.load(d / "positions.npy")
    cm = np.load(d / "genetic_pos.npy")
    allowed = hap != -1

    mask = genmap.chain_mask(positions, cm, None, None, support)
    mask, traversal = genmap.drop_untraversable(mask, cm, allowed, positions)

    # the same assertion `genmap.build` makes: an argument is not a check.
    left = genmap.untraversable_blocks(cm[mask], allowed[:, mask])
    if left:
        raise AssertionError(
            f"{d.name}/{held_out}: {len(left)} untraversable block(s) remain after "
            f"decision K, so the forward pass would dead-end. First at chain index "
            f"{left[0]}.")
    if not mask.any():
        raise AssertionError(f"{d.name}/{held_out}: the leave-one-out chain is empty")

    return LooChain(held_out=held_out, chrom=d.name, rows_dropped=tuple(int(r) for r in rows),
                    k_states=int(keep.sum()), support=support, in_chain=mask,
                    traversal=traversal, n_sites=int(len(positions)))


def masked_haplotypes(chrom_dir: str | Path, held_out: str, *,
                      haplotypes: np.ndarray | None = None) -> np.ndarray:
    """The (K - 2, n_sites) donor panel with `held_out`'s two rows gone."""
    d = Path(chrom_dir)
    ids = np.load(d / "haplotype_ids.npy", allow_pickle=True)
    hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r")) \
        if haplotypes is None else haplotypes
    rows = sample_index(ids)[held_out]
    keep = np.ones(hap_all.shape[0], dtype=bool)
    keep[list(rows)] = False
    return hap_all[keep]


def loo_shares(t_by_chrom: dict[str, int]) -> dict[str, float]:
    """`beta_share` per chromosome for ONE held-out individual, by variant count.

    ⚠⚠ THE SHARES ARE PER TARGET AND THEY MOVE. Decision J sets `beta_share` to
    `T_c / T_total`, and leave-one-out changes every `T_c`, so each of the 44 targets has
    its own 22 shares. Using arm 2's shares would mean `sum_v beta_v != 1` over the scope
    `tau` covers, which standing constraint 2 says does not weaken the bound but removes
    it.
    """
    if set(t_by_chrom) != set(AUTOSOMES):
        missing = sorted(set(AUTOSOMES) - set(t_by_chrom))
        extra = sorted(set(t_by_chrom) - set(AUTOSOMES))
        raise ValueError(
            f"shares need all 22 autosomes; missing {missing}, unexpected {extra}. "
            f"A share computed over a subset silently makes `tau` per-subset.")
    total = sum(t_by_chrom.values())
    if total <= 0:
        raise ValueError(f"total chain length is {total}; shares are undefined")
    shares = {c: t_by_chrom[c] / total for c in AUTOSOMES}
    s = sum(shares.values())
    if abs(s - 1.0) > 1e-12:
        raise AssertionError(f"shares sum to {s!r}, not 1; the calibration is void")
    return shares
