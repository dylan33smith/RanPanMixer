"""`baseline_model`: the cohort haplotype chain, built from G and D and nothing else.

This is `R_D`, the target-independent distribution the whole guarantee rests on.
Li-Stephens copying: the chain walks the chromosome, at each position occupies one
donor haplotype, and emits whatever that donor carries. Genetic distance decides
how readily it switches donors.

⚠ NO FUNCTION IN THIS MODULE TAKES A TARGET, and a test asserts it by reading the
source. Standing constraint 1 calls target leakage into `baseline_model` "the
easiest way to buy utility and silently lose the guarantee"; `cohort/` enforces
that for a whole subpackage by having no target parameter anywhere, and the same
device is available here because the tilt lives in `utility.py` instead.

THE CONSTANTS ARE RE-DERIVED FROM THE PUBLISHED METHODS, NOT FROM PANMIXER'S CODE
(`paper/s41467-026-77591-0_reference.pdf`, PDF page 10):

    Ne = 10,000,  r = 1.26,  Delta_x in centiMorgans,  d = Delta_x * Ne * r,
    p_stay = e^(-d/n) + (1 - e^(-d/n)) / n

PanMixer's shipped sampler gets four separate things wrong here -- its
`EFFECTIVE_N` is 1/10,000 rather than 10,000, its "distance" is a difference of
VCF row indices rather than centiMorgans, its forward pass and its sampler use
`d` values a factor 4 apart, and its non-anchor fill compares row indices against
base-pair positions. Measured consequence: total switch mass 5.8e-12 where the
published constants give 0.13-0.99, so the released mosaic is one donor copied
verbatim. `docs/bugs.md` has the full entry. Nothing here is inherited from it.

MISSINGNESS IS AN INDICATOR POTENTIAL, NOT A RENORMALISED TRANSITION ROW, and
this is a correction to the obvious reading of `missing_policy`. "Drop that donor
from the state distribution at that position" invites implementing it as: mask
row `j` of the transition to the allowed set and renormalise. That is not well
defined. Across a step with `Delta_x = 0` the transition is the identity, so a
donor that is available at `t-1` and unavailable at `t` has a masked row summing
to exactly ZERO -- there is no chain to renormalise, and the division is by zero.
chr21 has 618 interior map intervals with `dcM = 0` holding about 1,650 kept
variants, against 4.26% mean missingness, so this is reachable on the real data
rather than a contrived case.

So `R_D` is the Li-Stephens chain CONDITIONED on never occupying an unavailable
donor, the conditioning enters as a 0/1 potential, and the forward pass
normalises it away. ⚠ THE CONSEQUENCE IS IN `ffbs.py`: the conditioning constant
is part of `sum_t log c_t`, so `log_z_p` is that sum MINUS the baseline's own,
and forgetting the subtraction silently breaks the `[0, eta_tau]` assertion --
the cheapest correctness check we have.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

# The published Methods. ⚠ `exact.py` states these again rather than importing
# them, so the two implementations can disagree about the formula.
NE = 10_000.0
R = 1.26


@dataclass(frozen=True)
class CohortHMM:
    """`R_D` as the three things a forward pass needs, and nothing more.

    `allowed` is the state mask -- `missing_policy` = RENORMALISE as data.
    `switch_mass` is `s_t`, the coefficient of the all-ones part of
    `A_t = (1 - s_t) I + (s_t / K) J`. ⚠ It is the TOTAL switch mass, not the
    per-state one: `stay_switch` in `docs/terms.md` warns that confusing them
    rescales every transition by K. The per-state value is `s_t / K` and the
    probability of actually leaving state `j` is `s_t (K-1) / K`.

    ⚠ `A` IS NEVER MATERIALISED. At chr21 scale a dense matrix per position is
    340,849 * 88 * 88 * 8 B = about 21 TB. The rank-one-plus-diagonal structure
    makes the forward step O(K) from one scalar, and `ffbs.py` is the only place
    that arithmetic lives.
    """

    allowed: np.ndarray       # (K, T) bool     -- occupiable donors per CHAIN UNIT
    switch_mass: np.ndarray   # (T-1,) float64  -- s_t, between unit t and t+1
    genetic_pos: np.ndarray   # (T,) float64    -- cM at each unit's first position
    positions: np.ndarray     # (T,) int64      -- bp at each unit's first position
    cohort_id: str = ""
    # (T, 2) inclusive [start, end] chain-position span of each unit. MODEL B.
    # None means one position per unit, which IS model C -- so there is one code
    # path and model C is the identity case of it, not a separate branch.
    unit_spans: np.ndarray | None = None

    @property
    def K(self) -> int:
        return self.allowed.shape[0]

    @property
    def T(self) -> int:
        return self.allowed.shape[1]

    @property
    def n_positions(self) -> int:
        """Chain POSITIONS covered, which is T under model C and more under B."""
        if self.unit_spans is None:
            return self.T
        return int((self.unit_spans[:, 1] - self.unit_spans[:, 0] + 1).sum())

    @property
    def unit_sizes(self) -> np.ndarray:
        """(T,) positions per unit. All ones under model C."""
        if self.unit_spans is None:
            return np.ones(self.T, dtype=np.int64)
        return (self.unit_spans[:, 1] - self.unit_spans[:, 0] + 1).astype(np.int64)

    @property
    def model(self) -> str:
        return "C" if self.unit_spans is None else "B"

    @property
    def rho(self) -> np.ndarray:
        """The initial distribution: uniform over donors.

        A cohort-only choice, and the only defensible one without an external
        panel -- any other weighting would need frequencies we deliberately do
        not build (`docs/data.md`, the dropped-inputs table).
        """
        return np.full(self.K, 1.0 / self.K, dtype=np.float64)

    @property
    def p_stay(self) -> np.ndarray:
        """(T-1,) probability of remaining on the same donor. For reporting only."""
        return 1.0 - self.switch_mass + self.switch_mass / self.K

    @property
    def expected_switches(self) -> float:
        """E[number of donor changes] along the chain under the prior.

        The replacement for mis-specified exit gate 3, which demanded switch mass
        of order 0.1-1 between adjacent positions and would therefore have
        REJECTED a correct model: adjacent chr21 variants sit about 12 bp apart,
        where the biologically right value is about 0.0026. This quantity is
        comparable against the analytic segment length `n / (Ne * r)` instead.
        """
        return float((self.switch_mass * (self.K - 1) / self.K).sum())


def switch_mass_from_cm(genetic_pos: np.ndarray, n_states: int) -> np.ndarray:
    """`s_t = 1 - exp(-d_t / n)` with `d_t = Delta_x * Ne * r`, cross-checked.

    The decomposition `A = (1 - s) I + (s / K) J` and the published closed form
    `p_stay = e^(-d/n) + (1 - e^(-d/n))/n` are the same statement; computing both
    and comparing is a transcription check on the algebra that costs one array.
    It can genuinely fail -- which is the property `support.py` insists on and
    PanMixer's own equivalent assertion does not have.
    """
    if n_states < 2:
        raise ValueError(f"need at least 2 donor states, got {n_states}")
    if not np.isfinite(genetic_pos).all():
        bad = np.flatnonzero(~np.isfinite(genetic_pos))
        raise ValueError(
            f"{len(bad)} positions carry no interpolated cM (first at index "
            f"{bad[0]}). Inside chain_span every position must have a real cM; "
            f"fix the span, never clamp -- see genmap.chain_mask."
        )
    delta = np.diff(genetic_pos.astype(np.float64))
    if (delta < 0).any():
        i = np.flatnonzero(delta < 0)
        raise ValueError(f"cM decreases at index {i[0]}; the map is not monotonic")

    d = delta * NE * R
    decay = np.exp(-d / n_states)
    s = 1.0 - decay

    p_stay = 1.0 - s + s / n_states
    published = decay + (1.0 - decay) / n_states
    off = np.flatnonzero(~np.isclose(p_stay, published, rtol=0, atol=1e-15))
    if len(off):
        i = int(off[0])
        raise AssertionError(
            f"the rank-one decomposition disagrees with the published closed form "
            f"at index {i}: {p_stay[i]!r} vs {published[i]!r}"
        )
    return s


def from_arrays(haplotypes: np.ndarray, genetic_pos: np.ndarray,
                positions: np.ndarray, *, in_chain: np.ndarray | None = None,
                blocks: np.ndarray | None = None,
                cohort_id: str = "") -> CohortHMM:
    """Build `R_D` from cohort arrays. ⚠ There is no target parameter.

    `in_chain` is the `chain_span` mask over the FULL site axis. It is applied
    here rather than upstream because the cut is provisional: it is stored as a
    mask, never as a row filter, so moving it rewrites one small array and leaves
    every other artifact's row numbering -- and the site-axis digest -- untouched.

    `blocks` selects the CHAIN ARCHITECTURE and is the model-B fallback that
    `docs/plan.md` retains for closer parity with PanMixer:

    - `None` -> **MODEL C**, the decided architecture: one unit per variant.
    - a block id per chain position -> **MODEL B**: consecutive positions sharing
      an id become ONE ATOMIC unit. One donor is copied across the whole unit, so
      the chain can only switch donors BETWEEN units. Positions in no called block
      become their own singleton unit.

    ⚠ MODEL B'S THREE CONSEQUENCES, all of them measured and none of them hidden.
    (1) A donor may occupy a unit only if it is called at EVERY position in it --
    that is what atomic means -- so `allowed` is an AND across the unit rather than
    a per-position mask. Measured on chr21: the median unit still has all 88 donors
    callable, and only 6 units of 74,670 have none, which are excluded from the
    chain exactly as `support == 0` positions are. (2) The transition uses the GAP
    between units, from the last position of one to the first of the next, so the
    within-unit genetic length is discarded -- about half the recombination, which
    is the documented cost of model B. (3) The release carries a donor's real
    haplotype VERBATIM across each unit, up to 618 variants over 39.5 kb on chr21.
    That is departure (b), and `tau` bounds what the TARGET leaks, not the donors.

    ⚠ Where the blocks come from matters as much as the architecture. They must be
    a function of `G` and `D` or of a panel DECLARED IN ADVANCE that excludes the
    target -- see `scripts/build_blocks_chr21.sh`. PanMixer's own block set was
    computed on all 3,202 panel samples, which contains our test targets.
    """
    if haplotypes.ndim != 2:
        raise ValueError(f"haplotypes must be (K, n_sites), got {haplotypes.shape}")
    n_sites = haplotypes.shape[1]
    for name, arr in (("genetic_pos", genetic_pos), ("positions", positions)):
        if len(arr) != n_sites:
            raise ValueError(
                f"{name} has {len(arr)} entries but haplotypes has {n_sites} sites. "
                f"These were built against different VCFs."
            )
    if in_chain is None:
        in_chain = np.ones(n_sites, dtype=bool)
    elif len(in_chain) != n_sites:
        raise ValueError(f"in_chain has {len(in_chain)} entries, expected {n_sites}")

    hap = haplotypes[:, in_chain]
    gpos = genetic_pos[in_chain].astype(np.float64)
    pos = positions[in_chain].astype(np.int64)
    allowed = hap != -1
    unit_spans = None

    if blocks is not None:
        blocks = np.asarray(blocks)
        if len(blocks) == len(in_chain):
            blocks = blocks[in_chain]
        if len(blocks) != hap.shape[1]:
            raise ValueError(
                f"blocks has {len(blocks)} entries; expected {hap.shape[1]} chain "
                f"positions or {len(in_chain)} site-axis rows")
        # consecutive runs of one id are one unit; a change of id starts a new one
        cut = np.flatnonzero(np.diff(blocks) != 0) + 1
        starts = np.concatenate(([0], cut))
        ends = np.concatenate((cut - 1, [len(blocks) - 1]))
        unit_spans = np.stack([starts, ends], axis=1).astype(np.int64)
        # atomic: a donor may occupy a unit only if called at EVERY position in it
        allowed = np.stack([allowed[:, a:b + 1].all(axis=1) for a, b in unit_spans],
                           axis=1)
        keep = allowed.any(axis=0)
        if not keep.all():
            unit_spans = unit_spans[keep]
            allowed = allowed[:, keep]
        # ⚠ THE TRANSITION USES THE GAP, last position of one unit to first of the
        # next -- NOT start-to-start. The first version of this code used
        # start-to-start, which silently attributes a unit's own internal genetic
        # length to the transition out of it and contradicted the docstring above.
        # The gap form is what `docs/plan.md` specifies for model B and is why the
        # architecture "discards the within-block recombination".
        gaps = gpos[unit_spans[1:, 0]] - gpos[unit_spans[:-1, 1]]
        if (gaps < 0).any():
            raise ValueError("units overlap in genetic position; blocks are not disjoint")
        unit_cm_for_transitions = np.concatenate(([0.0], np.cumsum(gaps)))
        gpos = gpos[unit_spans[:, 0]]
        pos = pos[unit_spans[:, 0]]

    dead = np.flatnonzero(~allowed.any(axis=0))
    if len(dead):
        raise ValueError(
            f"{len(dead)} chain positions have NO called donor (first at chain "
            f"index {dead[0]}, pos {pos[dead[0]]}). RENORMALISE has an empty pool "
            f"there and the sampler has nothing to emit. `genmap.chain_mask` "
            f"excludes support == 0 upstream, so this means the mask and the "
            f"matrix have drifted apart -- rebuild, do not fall back."
        )

    # genetic_pos carries each unit's REAL first-position cM for provenance; the
    # transitions come from the inter-unit gaps, which are a different array.
    cm_for_transitions = gpos if unit_spans is None else unit_cm_for_transitions
    return CohortHMM(allowed=allowed,
                     switch_mass=switch_mass_from_cm(cm_for_transitions, hap.shape[0]),
                     genetic_pos=gpos, positions=pos, cohort_id=cohort_id,
                     unit_spans=unit_spans)


def from_cohort_dir(cohort_dir: str | Path, *, verify: bool = True) -> CohortHMM:
    """Load `R_D` from a cohort artifact directory. ⚠ Still no target parameter.

    `verify` re-checks each artifact against its recorded digest before use. The
    sidecar's sha256 was write-only until 2026-09-29 -- nothing compared it -- so
    a silently mutated artifact passed every guard; a digest nobody checks is
    decoration. It also refuses artifacts built against different site axes,
    which is the quiet failure: a re-sorted VCF relabels which variant is which
    and every downstream number is then about the wrong variants.
    """
    import store   # local: keeps this module importable without the src root on sys.path

    d = Path(cohort_dir)
    names = ["haplotypes.npy", "genetic_pos.npy", "positions.npy", "in_chain.npy"]
    if verify:
        digest = store.load_prov(d / "haplotypes.npy").get("axis_digest")
        for n in names:
            store.verify(d / n)
            if digest:
                store.check_axis(d / n, digest)
    return from_arrays(np.load(d / "haplotypes.npy"),
                       np.load(d / "genetic_pos.npy"),
                       np.load(d / "positions.npy"),
                       in_chain=np.load(d / "in_chain.npy"),
                       cohort_id=str(d))
