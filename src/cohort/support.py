"""Step 4: support(v), and the free row-alignment check that comes with it.

`support(v)` is the count of non-missing haplotype entries at site v. It is the
denominator of the per-variant utility weight w_v = 1/support(v).

It needs no computation at all: measured 2026-09-28, the graph VCF's own INFO/AN
field is bit-identical to `(haplotypes != -1).sum(axis=0)` across all 340,824
chr21 records. We compute it from the matrix anyway and compare -- because two
independently-derived numbers agreeing is a row-alignment check that costs
nothing, and a disagreement means the matrix and the VCF have drifted apart.

⚠ The check must stay genuinely independent. PanMixer's equivalent assertion is
vacuous -- both sides equal site_mask.sum() by construction -- so it can never
fire. Here one side comes from the parsed genotypes and the other from a field
the caller wrote; they can disagree, which is the whole point.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import store


def build(cohort_dir: Path, *, strict: bool = False) -> dict:
    cohort_dir = Path(cohort_dir)
    hap = np.load(cohort_dir / "haplotypes.npy")
    an = np.load(cohort_dir / "an_info.npy")
    support = (hap != -1).sum(axis=0).astype(np.int32)

    have_an = an >= 0
    mismatch = np.flatnonzero(have_an & (support != an))
    if strict and len(mismatch):
        raise ValueError(
            f"{len(mismatch)} sites where INFO/AN disagrees with the matrix "
            f"(first at row {mismatch[0]}). The arrays and the VCF have drifted."
        )

    meta = store.load_prov(cohort_dir / "haplotypes.npy")
    store.save(cohort_dir / "support.npy", support,
               inputs={"haplotypes": meta["sha256"]},
               params={"definition": "count of non-missing haplotype entries per site"},
               axis_digest=meta.get("axis_digest"),
               note="w_v = 1/support(v); cross-checked against INFO/AN")

    return {
        "n_sites": int(len(support)),
        "sites_with_AN": int(have_an.sum()),
        "an_mismatches": int(len(mismatch)),
        "an_mismatch_rows": [int(i) for i in mismatch[:20]],
        "support_min": int(support.min()), "support_max": int(support.max()),
        "support_mean": float(support.mean()),
        "sites_with_zero_support": int((support == 0).sum()),
    }
