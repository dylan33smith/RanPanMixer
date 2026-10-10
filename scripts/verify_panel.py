#!/usr/bin/env python
"""Check the 1000G 30x panel is sound, and that it is on OUR reference build.

⚠ WHY THE BUILD CHECK IS THE POINT. Decision F is a commitment not to inherit
PanMixer's reference mismatch: its released pipeline wires in 1000G Phase 3, which is
GRCh37, against a GRCh38 pangenome, and the cost is measurable -- 0.30% of chr21
pangenome records match Phase 3 on (POS, REF, ALT) against 75.88% for this panel, a
252x difference attributable to the build alone. A wrong-build panel does not fail
loudly. It produces a LOW MATCH RATE, which reads as biology.

So this script does not compare the panel against reference lengths typed in by hand.
It compares the panel's `##contig` lengths against the SAME lengths declared by our
own `vg deconstruct` output. Two independent files on disk, both describing the
assembly they were built on, and they must agree. A GRCh37 panel would disagree on
every chromosome (GRCh37 chr1 is 249,250,621 against GRCh38's 248,956,422), and no
constant in this file would need to be right for that to be caught.

The other checks cover how a 29 GB download actually goes wrong: a truncated file,
a missing index, a missing checksum, a leftover `.part`, or one chromosome carrying a
different sample set from the rest.

    python scripts/verify_panel.py [--deep]     # --deep re-reads 29 GB to verify sha256
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PANEL = Path("/data/ds85/RanPanMixer/data/raw")
# chr21 is excluded: it predates this fetch and is on disk under a different name.
CHROMS = [c for c in range(1, 23) if c != 21]


def _sh(cmd: list[str]) -> str:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd[:3])}... exited {r.returncode}: "
                           f"{r.stderr.strip()[:200]}")
    return r.stdout


def panel_vcf(c: int) -> Path:
    return PANEL / f"1000g_30x_chr{c}.vcf.gz"


def our_contig_length(c: int) -> int | None:
    """chrN's length as OUR OWN deconstruct output declares it.

    ⚠ Read from the header only -- `grep -m1` stops at the first match, so this does
    not stream the multi-GB body. The id is PanSN (`grch38#chrN`), not `chrN`.
    """
    f = REPO / "data" / "temp" / f"chr{c}" / "deconstruct.vcf"
    if not f.exists():
        return None
    with open(f, "r") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            m = re.match(rf"##contig=<ID=\S*chr{c},length=(\d+)>", line)
            if m:
                return int(m.group(1))
    return None


def panel_contig_length(c: int) -> int | None:
    """chrN's length as the index declares it -- index only, body untouched."""
    out = _sh(["bcftools", "index", "-s", str(panel_vcf(c))])
    for row in out.splitlines():
        parts = row.split("\t")
        if len(parts) >= 3 and parts[0] in (f"chr{c}", str(c)):
            return int(parts[1])
    return None


def check(c: int, *, deep: bool = False) -> tuple[list[str], dict]:
    bad: list[str] = []
    f = panel_vcf(c)
    info: dict = {"chrom": c}

    # 1. the three files a quotable artifact needs. A VCF with no checksum is not
    #    quotable and re-deriving one means re-reading 29 GB.
    for suffix, what in (("", "panel VCF"), (".tbi", "index"), (".sha256", "checksum")):
        p = Path(str(f) + suffix)
        if not p.exists() or p.stat().st_size == 0:
            bad.append(f"chr{c}: missing or empty {what} ({p.name})")
    if bad:
        return bad, info
    info["bytes"] = f.stat().st_size

    # 2. ⚠⚠ THE REFERENCE BUILD, against our own deconstruct output rather than a
    #    typed-in constant. This is the decision-F gate; a wrong-build panel produces
    #    a low match rate that reads as biology rather than an error.
    ours = our_contig_length(c)
    theirs = panel_contig_length(c)
    info["length"] = theirs
    if theirs is None:
        bad.append(f"chr{c}: the panel index declares no length for chr{c}")
    elif ours is None:
        bad.append(f"chr{c}: our own deconstruct.vcf is absent, so the panel's build "
                   f"CANNOT BE CHECKED -- do not treat this as a pass")
    elif ours != theirs:
        bad.append(f"chr{c}: REFERENCE BUILD MISMATCH. The panel declares chr{c} as "
                   f"{theirs:,} bp, our pangenome declares {ours:,} bp. These are "
                   f"different assemblies, and the symptom downstream would be a low "
                   f"(POS, REF, ALT) match rate that looks like biology, not an error.")

    # 3. record count, from the index.
    out = _sh(["bcftools", "index", "-n", str(f)]).strip()
    info["records"] = int(out) if out.isdigit() else None
    if not info["records"]:
        bad.append(f"chr{c}: the index reports {out!r} records")

    # 4. TRUNCATION, by SEEKING to the far end rather than streaming. A download cut
    #    short still opens, still indexes against its own partial content, and simply
    #    stops early -- so the last reachable record is the thing to look at.
    if theirs:
        lo = max(1, theirs - 3_000_000)
        body = _sh(["bcftools", "view", "-H", "-r", f"chr{c}:{lo}-{theirs}", str(f)])
        rows = [r for r in body.splitlines() if r.strip()]
        if not rows:
            bad.append(f"chr{c}: NO records in the last 3 Mb -- consistent with a "
                       f"truncated download")
        else:
            last = int(rows[-1].split("\t")[1])
            info["last_pos"] = last
            info["reach_pct"] = 100.0 * last / theirs
            if last > theirs:
                bad.append(f"chr{c}: a record at {last:,} lies BEYOND the declared "
                           f"contig length {theirs:,}")
            elif info["reach_pct"] < 99.0:
                bad.append(f"chr{c}: the last record reaches only "
                           f"{info['reach_pct']:.2f}% of the contig")

    # 5. the sample set. One chromosome fetched from a different release would carry a
    #    different panel, and a per-chromosome sample set is not one cohort.
    info["samples"] = len(_sh(["bcftools", "query", "-l", str(f)]).split())

    if deep:
        import hashlib
        h = hashlib.sha256()
        with open(f, "rb") as fh:
            while block := fh.read(1 << 22):
                h.update(block)
        want = Path(str(f) + ".sha256").read_text().split()[0]
        if h.hexdigest() != want:
            bad.append(f"chr{c}: sha256 MISMATCH -- recorded {want[:16]}..., "
                       f"file hashes to {h.hexdigest()[:16]}...")
    return bad, info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", type=int, nargs="+", default=CHROMS)
    ap.add_argument("--deep", action="store_true",
                    help="re-read all 29 GB and verify every sha256 (slow)")
    a = ap.parse_args(argv)

    problems, rows = [], []
    for c in a.chroms:
        bad, info = check(c, deep=a.deep)
        problems.extend(bad)
        info["ok"] = not bad
        rows.append(info)

    # a leftover .part means a download that never finished, whatever the sentinel says
    parts = sorted(PANEL.glob("*.part"))
    if parts:
        problems.append(f"{len(parts)} leftover .part file(s), so at least one "
                        f"download did not finish: {', '.join(p.name for p in parts[:4])}")

    samples = {r.get("samples") for r in rows if r.get("samples")}
    if len(samples) > 1:
        problems.append(f"the chromosomes do NOT share one sample set: {sorted(samples)}. "
                        f"A per-chromosome panel is not one cohort.")

    print("| chr | GB | records | contig bp | last record | samples | checks |")
    print("|" + "---|" * 7)
    for r in rows:
        if not r.get("records"):
            print(f"| {r['chrom']} | ? | ? | ? | ? | ? | **FAIL** |")
            continue
        print(f"| {r['chrom']} | {r['bytes']/1e9:.2f} | {r['records']:,} | "
              f"{r['length']:,} | {r.get('last_pos', 0):,} "
              f"({r.get('reach_pct', 0):.2f}%) | {r['samples']:,} | "
              f"{'pass' if r['ok'] else '**FAIL**'} |")
    tot = sum(r.get("bytes", 0) for r in rows)
    recs = sum(r.get("records") or 0 for r in rows)
    print(f"\n{sum(1 for r in rows if r['ok'])} of {len(a.chroms)} chromosome(s) pass; "
          f"{tot/1e9:.1f} GB, {recs:,} records, "
          f"{sorted(samples)[0] if len(samples) == 1 else '?'} samples")
    print(f"reference build checked against our own deconstruct output on "
          f"{sum(1 for r in rows if r.get('length'))} chromosome(s)"
          + ("" if a.deep else "; checksums NOT re-verified (use --deep)"))
    if problems:
        print("\nPROBLEMS — do not quote a number from this panel:")
        for p in problems:
            print(f"  {p}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
