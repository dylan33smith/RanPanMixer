#!/usr/bin/env python
"""B-DAT-map-coverage: how much of each autosome the genetic map actually speaks to.

⚠ WHY THIS HAS TO RUN BEFORE THE RATE IS DECLARED. Decision C (2026-10-08) stops
excluding the regions where the map gives no distance and samples them under a
DECLARED CONSTANT RECOMBINATION RATE instead. The rate's value is not chosen, and
it cannot be: nothing in this repository measured how much sequence the rule
applies to. Declaring a rate for a region whose size has never been measured is
the wrong order, so this script measures it first.

THREE QUANTITIES, REPORTED SEPARATELY AND NEVER ADDED (see `map_coverage` in
`docs/terms.md`):

  (a) MAPPED       -- inside the marker span. The map speaks here.
  (b) MAP-SILENT   -- before the first marker or after the last. The map says
                      nothing, `np.interp` would clamp, and this is the class
                      decision C's declared rate covers.
  (c) INTERIOR ZERO-cM -- inside the span, where the map speaks and says ZERO.
                      Decision C deliberately does NOT override these: the repo's
                      position is that P(switch) = 0 between two variants with no
                      genetic separation is the CORRECT model, and
                      `tests/test_toy_enumeration.py` enforces it as a Theorem-1
                      support property. Reported so the deferred question (open
                      question 2) has a size.

⚠ (b) AND (c) ARE DIFFERENT KINDS OF THING and the whole value of this table is
keeping them apart. Adding them would overstate what the declared rate reaches by
roughly a factor of two, which is exactly the over-claim an earlier draft of the
`chain_span` glossary entry made.

⚠ WHAT THIS SCRIPT CANNOT MEASURE: the N-gap spans. Whether the declared rate
should accumulate distance across the fabricated N-gaps of GRCh38 placeholder
sequence needs the reference assembly's gap coordinates, and no FASTA, AGP or
assembly report is on disk. Reported as BLOCKED rather than skipped.

Base-pair figures need only the maps. POSITION counts need each chromosome's
deconstructed VCF, so they are emitted only for chromosomes whose cohort stage has
run -- chr21 today -- and the rest arrive with preprocessing.

    python scripts/map_coverage.py [--chrom 21] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
MAPS = REPO / "external" / "genetic_maps"
COHORT = REPO / "data" / "cohort"

# GRCh38 primary assembly chromosome lengths. Hard-coded because no FASTA, AGP or
# assembly report is on disk; chr21's value is cross-checked below against the
# length the project already recorded from the HPRC VCF header (46,709,983), which
# is the one independent check available on this table.
GRCH38_LEN = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}
ACROCENTRIC = (13, 14, 15, 21, 22)

# The superseded coordinate cut, kept only to reconcile against what is on record.
CHR21_CUT_LO, CHR21_CUT_HI = 12_968_320, 46_680_243


def read_map(chrom: int) -> tuple[np.ndarray, np.ndarray]:
    """`(pos, cM)` for one autosome, with the map's own invariants asserted.

    A map that is not sorted, or whose cM goes backwards, would silently produce
    negative `Delta_x` downstream and a transition probability outside [0,1]. The
    real pipeline asserts monotonicity too; asserting it here as well means a bad
    map is caught while measuring rather than while releasing.
    """
    f = MAPS / f"chr{chrom}.b38.gmap"
    if not f.exists():
        raise FileNotFoundError(
            f"{f} is not extracted. The SHAPEIT4 archive on disk holds all 22 "
            f"autosomes plus chrX -- extract the autosomes only (decision D)."
        )
    pos, cm = [], []
    with open(f) as fh:
        header = fh.readline().split()
        if header[:1] != ["pos"] or header[-1] != "cM":
            raise ValueError(f"{f}: unexpected header {header!r}")
        for line in fh:
            p = line.split()
            pos.append(int(p[0]))
            cm.append(float(p[-1]))
    pos_a, cm_a = np.asarray(pos, dtype=np.int64), np.asarray(cm, dtype=np.float64)
    if pos_a.size < 2:
        raise ValueError(f"chr{chrom}: {pos_a.size} markers is not a map")
    if not np.all(np.diff(pos_a) > 0):
        raise ValueError(f"chr{chrom}: marker positions are not strictly increasing")
    if not np.all(np.diff(cm_a) >= 0):
        raise ValueError(f"chr{chrom}: cM is not monotonic non-decreasing")
    return pos_a, cm_a


def measure(chrom: int) -> dict:
    pos, cm = read_map(chrom)
    length = GRCH38_LEN[chrom]
    if pos[-1] > length:
        raise ValueError(
            f"chr{chrom}: last marker {pos[-1]:,} is beyond the recorded "
            f"chromosome length {length:,}; the map and the assembly disagree")

    lead, trail = int(pos[0]), int(length - pos[-1])
    mapped_bp = int(pos[-1] - pos[0])

    flat = np.diff(cm) == 0.0
    interior_zero_bp = int(np.diff(pos)[flat].sum())
    interior_zero_n = int(flat.sum())

    total_cm = float(cm[-1] - cm[0])
    rate = total_cm / (mapped_bp / 1e6) if mapped_bp else float("nan")

    return {
        "chrom": chrom, "markers": int(pos.size), "length_bp": length,
        "map_first": int(pos[0]), "map_last": int(pos[-1]),
        "mapped_bp": mapped_bp,
        "lead_bp": lead, "trail_bp": trail, "map_silent_bp": lead + trail,
        "interior_zero_bp": interior_zero_bp,
        "interior_zero_intervals": interior_zero_n,
        "total_cM": total_cm, "mean_rate_cM_per_Mb": rate,
        "acrocentric": chrom in ACROCENTRIC,
    }


def positions(chrom: int, m: dict) -> dict | None:
    """Position-level counts, for chromosomes whose cohort stage has run.

    ⚠ Base pairs are not positions, and the rate is applied per POSITION. A
    chromosome can be 20% map-silent by sequence and carry far fewer than 20% of
    its variants there, because the map-silent regions are heterochromatic and
    variant-poor. So this is the number the declared rate actually acts on, and
    the bp table alone must not be read as a substitute.
    """
    d = COHORT / f"chr{chrom}"
    if not (d / "positions.npy").exists():
        return None
    p = np.load(d / "positions.npy")
    sup = np.load(d / "support.npy")
    ic = np.load(d / "in_chain.npy").astype(bool)
    silent_lo = p < m["map_first"]
    silent_hi = p > m["map_last"]
    out = {
        "records": int(p.size), "in_chain": int(ic.sum()),
        "below_map_start": int(silent_lo.sum()),
        "below_map_start_support_gt0": int((silent_lo & (sup > 0)).sum()),
        "above_map_end": int(silent_hi.sum()),
        "above_map_end_support_gt0": int((silent_hi & (sup > 0)).sum()),
        "support_eq0_anywhere": int((sup == 0).sum()),
    }
    if chrom == 21:
        # Reconciliation against the figures already on record, which is this
        # script's own falsification gate: if the split does not add up to the
        # recorded below-cut count, either the measurement or the record is wrong.
        in_span_plateau = (p >= m["map_first"]) & (p <= CHR21_CUT_LO)
        out["in_span_plateau"] = int(in_span_plateau.sum())
        out["in_span_plateau_support_gt0"] = int(
            (in_span_plateau & (sup > 0)).sum())
        out["below_cut_total"] = int((p <= CHR21_CUT_LO).sum())
        out["above_cut_total"] = int((p > CHR21_CUT_HI).sum())
        out["split_reconciles"] = bool(
            out["below_map_start"] + out["in_span_plateau"]
            == out["below_cut_total"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chrom", type=int, default=None,
                    help="one autosome; default is all 22")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    chroms = [a.chrom] if a.chrom else list(range(1, 23))
    rows = [measure(c) for c in chroms]
    pos_rows = {c: positions(c, r) for c, r in zip(chroms, rows)}

    # chr21's recorded length is the one external check on GRCH38_LEN
    if 21 in chroms:
        assert GRCH38_LEN[21] == 46_709_983, "chr21's recorded length changed"

    tot_len = sum(r["length_bp"] for r in rows)
    tot_mapped = sum(r["mapped_bp"] for r in rows)
    tot_lead = sum(r["lead_bp"] for r in rows)
    tot_trail = sum(r["trail_bp"] for r in rows)
    tot_silent = tot_lead + tot_trail
    tot_izero = sum(r["interior_zero_bp"] for r in rows)
    tot_cm = math.fsum(r["total_cM"] for r in rows)
    acro_lead = sum(r["lead_bp"] for r in rows if r["acrocentric"])

    print(f"## Map coverage per autosome — {len(rows)} chromosomes\n")
    print("⚠ bp only. Position counts need each chromosome's VCF; see the second "
          "table.\n")
    hdr = (f"| chr | markers | length Mb | lead Mb | trail Mb | map-silent % | "
           f"interior 0cM Mb | 0cM intervals | total cM | cM/Mb |")
    print(hdr)
    print("|" + "---|" * 10)
    for r in rows:
        mark = " *" if r["acrocentric"] else ""
        print(f"| {r['chrom']}{mark} | {r['markers']:,} | "
              f"{r['length_bp']/1e6:.1f} | {r['lead_bp']/1e6:.2f} | "
              f"{r['trail_bp']/1e6:.2f} | "
              f"{100*r['map_silent_bp']/r['length_bp']:.2f}% | "
              f"{r['interior_zero_bp']/1e6:.2f} | "
              f"{r['interior_zero_intervals']:,} | {r['total_cM']:.1f} | "
              f"{r['mean_rate_cM_per_Mb']:.3f} |")
    print("\n`*` = acrocentric (unassembled p arm).\n")

    print("### Autosomal totals — metrics are rows\n")
    print("| metric | value | share of autosome |")
    print("|---|---|---|")
    for label, bp in (("total length", tot_len), ("mapped span", tot_mapped),
                      ("map-silent (no cM at all)", tot_silent),
                      ("  of which before the first marker", tot_lead),
                      ("  of which after the last marker", tot_trail),
                      ("  of which the 5 acrocentric p arms", acro_lead),
                      ("interior zero-cM (inside the span)", tot_izero)):
        print(f"| {label} | {bp/1e6:,.1f} Mb | {100*bp/tot_len:.2f}% |")
    print(f"| total genetic length | {tot_cm:,.1f} cM | — |")
    print(f"| **autosomal mean mapped rate** | "
          f"**{tot_cm/(tot_mapped/1e6):.4f} cM/Mb** | — |")
    rates = np.array([r["mean_rate_cM_per_Mb"] for r in rows])
    print(f"| per-chromosome rate, min / median / max | "
          f"{rates.min():.3f} / {np.median(rates):.3f} / {rates.max():.3f} cM/Mb "
          f"| spread {rates.max()/rates.min():.2f}x |")

    print("\n### Position counts, where the cohort stage has run\n")
    have = {c: v for c, v in pos_rows.items() if v}
    if not have:
        print("None. No chromosome's cohort stage has produced `positions.npy`.")
    for c, v in have.items():
        print(f"**chr{c}** — {v['records']:,} records, {v['in_chain']:,} in chain")
        print(f"- below the map start: **{v['below_map_start']:,}** "
              f"({v['below_map_start_support_gt0']:,} with support > 0) "
              f"— the declared rate APPLIES")
        print(f"- above the map end: **{v['above_map_end']:,}** "
              f"({v['above_map_end_support_gt0']:,} with support > 0) "
              f"— the declared rate APPLIES")
        if "in_span_plateau" in v:
            print(f"- in-span zero-cM plateau: **{v['in_span_plateau']:,}** "
                  f"({v['in_span_plateau_support_gt0']:,} with support > 0) "
                  f"— the rate does NOT apply (open question 2)")
            print(f"- reconciliation: {v['below_map_start']:,} + "
                  f"{v['in_span_plateau']:,} = {v['below_cut_total']:,} below the "
                  f"superseded cut — **{'RECONCILES' if v['split_reconciles'] else 'DOES NOT RECONCILE'}**")
            if not v["split_reconciles"]:
                raise AssertionError(
                    "chr21's map-silent + in-span split does not equal the "
                    "below-cut count. Either this measurement or the recorded "
                    "34,371 is wrong; do not quote either until it is resolved.")
        print(f"- positions with support == 0 anywhere: "
              f"{v['support_eq0_anywhere']:,}\n")

    print("### ⚠ BLOCKED: N-gap spans\n")
    print("Whether the declared rate should accumulate distance across the "
          "fabricated N-gaps of GRCh38 placeholder sequence cannot be measured "
          "here: no reference FASTA, AGP or assembly report is on disk. A constant "
          "cM/Mb multiplied by fictitious base pairs yields fictitious cM, so the "
          "policy must be declared before any rate is applied. Needs the GRCh38 "
          "primary assembly's gap coordinates.")

    if a.json:
        a.json.write_text(json.dumps(
            {"per_chrom": rows, "positions": {str(k): v for k, v in pos_rows.items()},
             "totals": {"length_bp": tot_len, "mapped_bp": tot_mapped,
                        "map_silent_bp": tot_silent, "lead_bp": tot_lead,
                        "trail_bp": tot_trail, "acrocentric_lead_bp": acro_lead,
                        "interior_zero_bp": tot_izero, "total_cM": tot_cm,
                        "mean_rate_cM_per_Mb": tot_cm / (tot_mapped / 1e6)}},
            indent=1))
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
