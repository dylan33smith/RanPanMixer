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


def chain_mask(positions: np.ndarray, cm: np.ndarray, start: int, end: int) -> np.ndarray:
    """The chain_span mask, asserting no NaN survives inside it.

    `start` is exclusive and `end` inclusive, matching configs/inputs.chr21.yaml.
    The assertion is the payoff for interpolating with NaN: it makes "every chain
    position carries a genuine interpolated cM" a checked invariant rather than a
    claim in a document.
    """
    mask = (positions > start) & (positions <= end)
    bad = mask & ~np.isfinite(cm)
    if bad.any():
        i = np.flatnonzero(bad)
        raise ValueError(
            f"{bad.sum()} positions inside chain_span have no interpolated cM "
            f"(first at row {i[0]}, pos {positions[i[0]]}). The span and the map "
            f"disagree -- fix the span, do not clamp."
        )
    return mask


def build(positions_npy: Path, gmap_path: Path, outdir: Path,
          *, start: int, end: int) -> dict:
    """Emit genetic_pos.npy and in_chain.npy."""
    outdir = Path(outdir)
    positions = np.load(positions_npy)
    cm = interpolate(positions, gmap_path)
    mask = chain_mask(positions, cm, start, end)

    meta = store.load_prov(positions_npy)
    digest = meta.get("axis_digest")
    inputs = {"positions": meta["sha256"], "gmap": store.sha256_file(gmap_path)}
    params = {"chain_span_start_exclusive": start, "chain_span_end_inclusive": end,
              "interp": "numpy.interp, left=nan right=nan", "provisional": True}

    store.save(outdir / "genetic_pos.npy", cm, inputs=inputs, params=params,
               axis_digest=digest,
               note="cM per site; NaN outside the map span, by design")
    store.save(outdir / "in_chain.npy", mask, inputs=inputs, params=params,
               axis_digest=digest,
               note="the chain_span MASK over the full axis. PROVISIONAL: moving the "
                    "cut rewrites this one array and nothing else.")

    inside = cm[mask]
    d = np.diff(inside)
    return {
        "n_sites": int(len(positions)),
        "T": int(mask.sum()),
        "dropped_below": int((positions <= start).sum()),
        "dropped_above": int((positions > end).sum()),
        "cm_span": [float(inside.min()), float(inside.max())] if mask.any() else None,
        "zero_distance_steps_inside": int((d == 0).sum()),
        "max_step_cm": float(d.max()) if len(d) else None,
    }
