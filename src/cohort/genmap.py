"""Step 5: centiMorgans per site, and the chain_span MASK.

The mask is the single most important structural choice in the pipeline. The
coordinate cut is PROVISIONAL -- Dylan is settling it with his PI -- so it is
stored as a boolean over the full axis, never applied as a row filter. Moving the
cut then rewrites one small array, leaves the site-axis digest byte-identical,
and does not invalidate the target paths. Filtering rows instead would renumber
every artifact in the tree.
"""

from __future__ import annotations

from pathlib import Path

import re

import numpy as np

import store


def read_gmap(gmap_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a `pos / chr / cM` genetic map. Returns (positions, centimorgans)."""
    pos, cm = [], []
    with open(gmap_path) as fh:
        header = fh.readline().split()
        if header[:1] != ["pos"] or header[-1:] != ["cM"]:
            raise ValueError(f"{gmap_path}: expected a 'pos chr cM' header, got {header}")
        for line in fh:
            if not line.strip():
                continue
            f = line.split()
            pos.append(int(f[0])); cm.append(float(f[-1]))
    p = np.asarray(pos, dtype=np.int64)
    c = np.asarray(cm, dtype=np.float64)
    if not np.all(np.diff(p) > 0):
        raise ValueError(f"{gmap_path}: positions are not strictly increasing")
    if not np.all(np.diff(c) >= 0):
        raise ValueError(f"{gmap_path}: cM is not monotonic non-decreasing")
    return p, c


def interpolate(positions: np.ndarray, gmap_path: Path) -> np.ndarray:
    """Linear interpolation onto the genetic map, with NaN outside its span.

    `left=nan, right=nan` deliberately, NOT numpy's default clamp. Clamping is
    what produced the defect the coordinate cut exists for: 34,345 consecutive
    chr21 positions sharing one cM value, hence P(switch) = 0, hence one donor's
    real haplotype emitted verbatim across 7.25 Mb. A NaN is loud; a clamp is
    silent and looks exactly like data.

    Linear is also what Beagle, SHAPEIT4, hap-ibd and GLIMPSE all do. A cubic
    spline would be actively wrong here -- overshoot can produce a negative
    genetic distance between two adjacent markers.
    """
    mp, mc = read_gmap(gmap_path)
    return np.interp(positions.astype(np.float64), mp.astype(np.float64), mc,
                     left=np.nan, right=np.nan)


def chain_mask(positions: np.ndarray, cm: np.ndarray, start: int | None = None,
               end: int | None = None,
               support: np.ndarray | None = None) -> np.ndarray:
    """Which positions are in the chain, asserting no NaN survives among them.

    ⚠ `start` and `end` ARE NOW OPTIONAL, AND LEAVING THEM OUT IS THE POST-DECISION-C
    BEHAVIOUR. The coordinate cut was a way of excluding the positions no genetic
    map reaches; decision C includes them under a declared rate instead, so there
    is nothing left for a cut to do. With both omitted, the only exclusion is
    `support == 0`, which is a property of the cohort and not of the map.
    Passing a span still works, and is how the pre-decision-C chr21 artifacts were
    built -- they are superseded, not deleted.

    `start` is exclusive and `end` inclusive, matching configs/inputs.chr21.yaml.
    The assertion is the payoff for interpolating with NaN: it makes "every chain
    position carries a genuine interpolated cM" a checked invariant rather than a
    claim in a document.

    ⚠ Positions with `support == 0` are also excluded. There, NO donor haplotype
    has a called allele, so the sampler has nothing to emit and RENORMALISE has an
    empty pool -- the position is unusable, not merely awkward. Excluding it also
    removes the `w_v = 1/support(v)` division by zero at its source, which is
    better than patching the weight: `1/(support+1)` would keep the position in
    the chain and give it a finite weight it cannot ever earn. The exclusion is a
    function of the cohort alone, so it is target-independent.
    """
    mask = np.ones(len(positions), dtype=bool)
    if start is not None:
        mask &= positions > start
    if end is not None:
        mask &= positions <= end
    if support is not None:
        mask &= support > 0
    bad = mask & ~np.isfinite(cm)
    if bad.any():
        i = np.flatnonzero(bad)
        raise ValueError(
            f"{bad.sum()} chain positions have no cM (first at row {i[0]}, pos "
            f"{positions[i[0]]}). Under decision C the map-silent positions are "
            f"FILLED by `fill_map_silent` with the declared rate, not excluded -- "
            f"so this means the fill was not applied. Do not clamp, and do not "
            f"re-introduce a coordinate cut."
        )
    return mask


def _chrom_of(gmap_path: Path) -> str:
    """`chr7` from `.../chr7.b38.gmap`. The gap table keys on this exact spelling.

    ⚠ Derived from the map's filename rather than passed in, so the gap table and
    the genetic map cannot disagree about which chromosome is being built -- a
    mismatch there would subtract one chromosome's gaps from another's distances
    and nothing downstream would notice.
    """
    name = Path(gmap_path).name
    m = re.match(r"(chr[0-9XYMT]+)\.", name)
    if not m:
        raise ValueError(
            f"cannot read a chromosome name from {name!r}; expected something like "
            f"'chr7.b38.gmap' so the gap table can be keyed on it")
    return m.group(1)


def build(positions_npy: Path, gmap_path: Path, outdir: Path,
          *, start: int | None = None, end: int | None = None,
          gaps_path: Path | None = None,
          declared_rate: float | None = None) -> dict:
    """Emit `genetic_pos.npy`, `in_chain.npy`, and (post-decision-C) `cm_source.npy`.

    TWO MODES, and the second is the current one.

    **Decision C (pass `gaps_path` and `declared_rate`, omit the span).** Positions
    the genetic map does not reach are FILLED from the declared constant rate over
    gap-free base pairs, and the chain covers the whole assembled chromosome less
    only `support == 0`. `cm_source.npy` records, per position, whether its cM came
    from the map or from the declared rate, so any release can report what share of
    its chain ran on a declared distance.

    **Pre-decision-C (pass `start` and `end`).** The coordinate cut, which excluded
    the map-silent positions instead of filling them. Retained because the chr21
    artifacts on disk were built this way and are superseded rather than deleted;
    it must not be used for a new build.
    """
    outdir = Path(outdir)
    positions = np.load(positions_npy)
    cm = interpolate(positions, gmap_path)

    source = None
    if (gaps_path is None) != (declared_rate is None):
        raise ValueError(
            "pass `gaps_path` AND `declared_rate` together or neither. A rate with "
            "no gap table would accumulate distance across placeholder base pairs, "
            "which is the thing the gap table was acquired to prevent.")
    if declared_rate is not None:
        if start is not None or end is not None:
            raise ValueError(
                "a coordinate cut and the declared-rate fill are alternatives, not "
                "a pair. Decision C replaces the cut; passing both would fill the "
                "map-silent positions and then throw them away.")
        gaps = read_gaps(Path(gaps_path), _chrom_of(gmap_path))
        cm, source = fill_map_silent(positions, cm, gaps, declared_rate)

    sup_path = Path(outdir) / "support.npy"
    support = np.load(sup_path) if sup_path.exists() else None
    mask = chain_mask(positions, cm, start, end, support)

    meta = store.load_prov(positions_npy)
    digest = meta.get("axis_digest")
    if support is not None and digest:
        store.check_axis(sup_path, digest)   # the guard must actually be invoked
    inputs = {"positions": meta["sha256"], "gmap": store.sha256_file(gmap_path)}
    params = {"chain_span_start_exclusive": start, "chain_span_end_inclusive": end,
              "interp": "numpy.interp, left=nan right=nan",
              "provisional": declared_rate is None,
              "declared_rate_cM_per_Mb": declared_rate,
              "gap_spans_skipped": gaps_path is not None,
              "excludes_zero_support": support is not None}

    store.save(outdir / "genetic_pos.npy", cm, inputs=inputs, params=params,
               axis_digest=digest,
               note="cM per site; NaN outside the map span, by design")
    store.save(outdir / "in_chain.npy", mask, inputs=inputs, params=params,
               axis_digest=digest,
               note="the chain MASK over the full axis. Under decision C the only "
                    "exclusion is support == 0; a coordinate cut is superseded.")
    if source is not None:
        store.save(outdir / "cm_source.npy", source, inputs=inputs, params=params,
                   axis_digest=digest,
                   note="0 = cM from the genetic map, 1 = from the declared constant "
                        "rate. Saved so a release can report the share of its chain "
                        "running on a declared rather than a measured distance.")

    inside = cm[mask]
    d = np.diff(inside)
    return {
        "n_sites": int(len(positions)),
        "T": int(mask.sum()),
        # ⚠ Reported against the span ACTUALLY applied. Under decision C there is
        # no span, so "dropped" is 0 by construction rather than unmeasured, and
        # the exclusion count is the whole-axis support == 0 count.
        "excluded_zero_support": (int((support == 0).sum()) if start is None and end is None
                                  else int(((positions > start) & (positions <= end)
                                            & (support == 0)).sum())) if support is not None else None,
        "dropped_below": 0 if start is None else int((positions <= start).sum()),
        "dropped_above": 0 if end is None else int((positions > end).sum()),
        "declared_rate_positions": int((source == 1).sum()) if source is not None else None,
        "declared_rate_positions_in_chain": (int((source[mask] == 1).sum())
                                             if source is not None else None),
        "cm_span": [float(inside.min()), float(inside.max())] if mask.any() else None,
        "zero_distance_steps_inside": int((d == 0).sum()),
        "max_step_cm": float(d.max()) if len(d) else None,
    }


# ---------------------------------------------------------------------------
# ASSEMBLY GAPS AND THE DECLARED RATE (decision C's value settled 2026-10-10;
# gap handling settled 2026-10-10 as SKIP)
# ---------------------------------------------------------------------------

def read_gaps(gap_path: Path, chrom: str) -> np.ndarray:
    """One chromosome's assembly gaps as a sorted (n, 2) array of half-open [s, e).

    Reads the UCSC `gap.txt.gz` schema: bin, chrom, chromStart, chromEnd, ix, n,
    size, type, bridge. `chromStart` is 0-based half-open, which is why the
    overlap arithmetic below is half-open too.

    ⚠ Why gaps are needed at all. The declared constant rate converts BASE PAIRS
    to centiMorgans, and in a gap the base pairs are placeholders -- a convention
    standing in for an unknown amount of missing sequence, not a measurement.
    GRCh38 chr21's p arm is 25 contigs separated by fabricated N-gaps, 22 of them
    exactly 50 kb, and "exactly 50 kb, 22 times" is the tell. Multiplying a rate
    by that span manufactures genetic distance out of a spacer.
    """
    import gzip

    out = []
    opener = gzip.open if str(gap_path).endswith(".gz") else open
    with opener(gap_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[1] != chrom:
                continue
            out.append((int(f[2]), int(f[3])))
    if not out:
        raise ValueError(
            f"{gap_path} lists no gaps for {chrom!r}. Every GRCh38 autosome has at "
            f"least its telomere gaps, so an empty result means the chromosome "
            f"name does not match the file's convention (it uses 'chr21', not '21')."
        )
    a = np.asarray(sorted(out), dtype=np.int64)
    if (a[:, 1] <= a[:, 0]).any():
        raise ValueError(f"{chrom}: a gap has non-positive length")
    if len(a) > 1 and (a[1:, 0] < a[:-1, 1]).any():
        raise ValueError(f"{chrom}: gaps overlap; the overlap arithmetic assumes "
                         f"they are disjoint")
    return a


def gap_free_span(lo: np.ndarray, hi: np.ndarray, gaps: np.ndarray) -> np.ndarray:
    """Base pairs between `lo` and `hi` that are NOT inside an assembly gap.

    Vectorised over pairs via a cumulative gap-length function, so it is O(n log n)
    rather than O(n * n_gaps): `G(x)` is the total gap base pairs below `x`, and the
    gap length inside [lo, hi) is `G(hi) - G(lo)`.

    ⚠ MEASURED 2026-10-10, and it is what makes skipping safe rather than merely
    defensible: **no variant of ours falls inside any gap** (chr21: 0 of 340,849
    records, and 0 of the 26,636 in the map-silent region). So gaps contain nothing
    we release, and subtracting them changes NO position's chain membership and no
    support count. The only thing that changes is `Delta_x` across a gap -- which is
    precisely the quantity that was wrong. Autosome-wide this takes the declared
    map-silent distance from 109.5 cM to 22.1 cM, of which 87.4 cM was placeholder.
    """
    if gaps.ndim != 2 or gaps.shape[1] != 2:
        raise ValueError(f"gaps must be (n, 2), got {gaps.shape}")
    edges = np.concatenate([[0], gaps[:, 0], gaps[:, 1]])
    order = np.argsort(edges, kind="stable")
    # cumulative gap length at each sorted edge
    cum = np.zeros(len(edges), dtype=np.int64)
    lens = np.concatenate([[0], np.zeros(len(gaps), dtype=np.int64),
                           gaps[:, 1] - gaps[:, 0]])
    cum[order] = np.cumsum(lens[order])

    def G(x: np.ndarray) -> np.ndarray:
        """Total gap bp strictly below `x`, with partial gaps counted up to `x`."""
        x = np.asarray(x, dtype=np.int64)
        starts, ends = gaps[:, 0], gaps[:, 1]
        # fully-passed gaps
        full = np.searchsorted(ends, x, side="right")
        total = np.concatenate([[0], np.cumsum(ends - starts)])[full]
        # the gap we may be standing inside
        idx = np.searchsorted(starts, x, side="right") - 1
        inside = (idx >= 0) & (x > starts[np.clip(idx, 0, None)]) & \
                 (x < ends[np.clip(idx, 0, None)])
        partial = np.where(inside, x - starts[np.clip(idx, 0, None)], 0)
        return total + partial

    lo = np.asarray(lo, dtype=np.int64)
    hi = np.asarray(hi, dtype=np.int64)
    if (hi < lo).any():
        raise ValueError("hi < lo in a span")
    span = hi - lo
    gap_bp = G(hi) - G(lo)
    free = span - gap_bp
    if (free < 0).any():
        raise AssertionError(
            "a gap-free span came out negative, so the cumulative gap function "
            "disagrees with the raw span -- the gap table is probably not disjoint")
    return free


def declared_cm(lo: np.ndarray, hi: np.ndarray, gaps: np.ndarray,
                rate_cm_per_mb: float) -> np.ndarray:
    """Genetic distance from the DECLARED constant rate, over gap-free base pairs.

    ⚠ This is used ONLY where the genetic map is silent. Inside the mapped span the
    cM comes from the map, whose markers sit at real positions, so interpolating
    between them already accounts for whatever is in between -- gaps included. The
    gap subtraction belongs exactly where distance is being manufactured from base
    pairs, and nowhere else.
    """
    if not (0.0 < rate_cm_per_mb < 100.0):
        raise ValueError(f"implausible declared rate {rate_cm_per_mb!r} cM/Mb")
    return gap_free_span(lo, hi, gaps) / 1e6 * float(rate_cm_per_mb)


def fill_map_silent(positions: np.ndarray, cm: np.ndarray, gaps: np.ndarray,
                    rate_cm_per_mb: float) -> tuple[np.ndarray, np.ndarray]:
    """Give the map-silent positions a cM from the declared rate.

    Returns `(cm_filled, source)` where `source` is 0 where the cM came from the
    genetic map and 1 where it came from the declared rate. The source array is
    saved beside the cM, so every release can report what share of its chain ran
    on a declared rather than a measured distance.

    ⚠ WHY A NaN IS THE RIGHT INPUT HERE. `interpolate` returns NaN outside the
    map's marker span deliberately, rather than clamping. A clamp gives thousands
    of positions one identical cM, which makes the switch probability zero and
    emits a donor's real haplotype verbatim; a NaN is loud. Decision C does not
    change that -- it changes what we DO with the loud signal. The NaN still marks
    exactly the positions no map speaks to, and this function is the only place
    allowed to replace one.

    Only the two ENDS can be NaN: `numpy.interp` interpolates everywhere between
    the first and last marker, so a long marker-free stretch inside the span gets
    a linear interpolation rather than a NaN. Those interior stretches are a
    separate question (the interior zero-cM plateaus) and decision C does not
    touch them.
    """
    cm = np.asarray(cm, dtype=np.float64).copy()
    positions = np.asarray(positions, dtype=np.int64)
    if cm.shape != positions.shape:
        raise ValueError(f"cm is {cm.shape}, positions {positions.shape}")
    finite = np.isfinite(cm)
    if not finite.any():
        raise ValueError(
            "no position has a mapped cM at all, so there is no anchor to measure "
            "the declared distance FROM. The map and this chromosome disagree.")

    source = np.where(finite, 0, 1).astype(np.int8)
    idx = np.flatnonzero(finite)
    lo_i, hi_i = int(idx[0]), int(idx[-1])

    # ⚠ Interior NaNs would break the two-sided fill below, so refuse rather than
    # guess. They cannot arise from numpy.interp and would mean the cM came from
    # somewhere else.
    if not finite[lo_i:hi_i + 1].all():
        n = int((~finite[lo_i:hi_i + 1]).sum())
        raise ValueError(
            f"{n} positions between the first and last mapped position have no cM. "
            f"numpy.interp cannot produce that, so this cM array did not come from "
            f"`interpolate` and the fill's assumptions do not hold.")

    below = np.arange(len(cm)) < lo_i
    above = np.arange(len(cm)) > hi_i
    if below.any():
        back = declared_cm(positions[below], np.full(below.sum(), positions[lo_i]),
                           gaps, rate_cm_per_mb)
        cm[below] = cm[lo_i] - back
    if above.any():
        fwd = declared_cm(np.full(above.sum(), positions[hi_i]), positions[above],
                          gaps, rate_cm_per_mb)
        cm[above] = cm[hi_i] + fwd

    if not np.isfinite(cm).all():
        raise AssertionError("a NaN survived the declared-rate fill")
    d = np.diff(cm)
    if (d < -1e-12).any():
        i = int(np.flatnonzero(d < -1e-12)[0])
        raise AssertionError(
            f"cM goes BACKWARDS at row {i} (pos {positions[i]} -> "
            f"{positions[i + 1]}: {cm[i]!r} -> {cm[i + 1]!r}). A negative genetic "
            f"distance makes the switch probability leave [0,1]. The positions are "
            f"probably not sorted, or two positions are equal.")
    return cm, source
