"""Step 5: centiMorgans per site, and the chain_span MASK.

The mask is the single most important structural choice in the pipeline. The
coordinate cut is PROVISIONAL -- Dylan is settling it with his PI -- so it is
stored as a boolean over the full axis, never applied as a row filter. Moving the
cut then rewrites one small array, leaves `site_axis_digest` byte-identical, and
does not invalidate the ten target paths. Filtering rows instead would renumber
every artifact in the tree.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np


def interpolate(positions: np.ndarray, gmap_path: Path) -> np.ndarray:
    """Linear interpolation onto the genetic map, with NaN outside its span.

    `left=nan, right=nan` deliberately, NOT numpy's default clamp. Clamping is
    what produced the defect this whole cut exists for: 34,345 consecutive chr21
    positions sharing one cM value, hence P(switch) = 0, hence one donor's real
    haplotype emitted verbatim across 7.25 Mb. A NaN is loud; a clamp is silent.
    """
    raise NotImplementedError


def chain_mask(positions: np.ndarray, cm: np.ndarray, start: int, end: int) -> np.ndarray:
    """The chain_span mask, asserting no NaN survives inside it.

    That assertion is the payoff for interpolating with NaN: it makes "every
    chain position carries a genuine interpolated cM" a checked invariant rather
    than a claim.
    """
    raise NotImplementedError
