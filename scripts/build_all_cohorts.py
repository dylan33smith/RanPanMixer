#!/usr/bin/env python
"""Build the cohort artifacts for all 22 autosomes. The 1..22 driver.

⚠ FAN-OUT IS BOUNDED BY MEMORY, NOT BY CORES, and the bound is measured rather than
guessed. `vg deconstruct` peaked at 9.48 GB RSS on chr21's 435 MB gzipped graph —
about 22 MB of RSS per MB of graph — so chr1's 1.2 GB graph implies roughly 26 GB.
Twenty-two at once would ask for several hundred GB against 251, and a fixed fan-out
either wastes the machine on small chromosomes or kills it on large ones.

So the chromosomes are BIN-PACKED largest-first: each batch's summed estimated RSS
fits the memory budget, and the batch runs at its own width. Small chromosomes end up
many-at-a-time, chr1 and chr2 end up nearly alone, and nothing has to be tuned by hand.

⚠ AND IT ASSERTS ITS OWN COMPLETENESS over all 22 chromosomes and every expected
artifact. A driver that filters results cannot tell "found nothing" from "ran nothing"
(`CLAUDE.md`), and a half-built genome looks exactly like a whole one to everything
downstream.

    python scripts/build_all_cohorts.py [--dry-run] [--chroms 1 2 21]
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import runner as R  # noqa: E402

# Measured on chr21: 9,476,712 KB of peak RSS for a 434,772,428 B gzipped GFA.
RSS_PER_GFA_BYTE = (9_476_712 * 1024) / 434_772_428        # ~22.3
EXPECTED = ("haplotypes.npy", "positions.npy", "support.npy", "genetic_pos.npy",
            "in_chain.npy", "cm_source.npy", "sites.tsv", "manifest.json")


def gfa_path(c: int) -> Path:
    return REPO / "data" / "raw" / f"chr{c}.hprc-v1.0-pggb.gfa.gz"


def estimate_rss(c: int) -> int:
    f = gfa_path(c)
    if not f.exists():
        raise FileNotFoundError(f"{f} is not downloaded")
    return int(f.stat().st_size * RSS_PER_GFA_BYTE)


def bin_pack(chroms: list[int], budget: int) -> list[list[int]]:
    """Largest-first batches whose summed estimated RSS fits `budget`.

    ⚠ A chromosome whose own estimate exceeds the budget still gets a batch of its
    own rather than being skipped. Dropping it would be the silent-incompleteness
    failure this driver exists to prevent; running it alone at least fails loudly
    if the machine cannot take it.
    """
    order = sorted(chroms, key=estimate_rss, reverse=True)
    batches: list[list[int]] = []
    cur: list[int] = []
    cur_sum = 0
    for c in order:
        need = estimate_rss(c)
        if cur and cur_sum + need > budget:
            batches.append(cur)
            cur, cur_sum = [], 0
        cur.append(c)
        cur_sum += need
    if cur:
        batches.append(cur)
    return batches


def _build_one(c: int) -> dict:
    """One chromosome, start to finish, through the real CLI."""
    cfg = REPO / "configs" / f"inputs.chr{c}.yaml"
    if not cfg.exists():
        raise FileNotFoundError(
            f"{cfg} does not exist. Run scripts/make_chrom_configs.py first; it "
            f"refuses to generate a partial set, so a missing config here means the "
            f"inputs are incomplete.")
    logdir = REPO / "logs" / "cohorts"
    logdir.mkdir(parents=True, exist_ok=True)
    log = logdir / f"chr{c}.log"
    t0 = time.time()
    # ⚠ Invoked as a MODULE with PYTHONPATH, not via the `prep-cohort` console
    # script. The script exists in pyproject.toml but only resolves if the package
    # has been pip-installed, which it is not in this environment -- and a driver
    # that depends on an install step nobody performed fails 20 chromosomes in.
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    with open(log, "w") as fh:
        p = subprocess.run([sys.executable, "-c",
                            "import cli, sys; sys.exit(cli.cohort_main())",
                            "all", "--chrom", f"chr{c}", "--config", str(cfg)],
                           stdout=fh, stderr=subprocess.STDOUT, cwd=REPO, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"chr{c}: cohort build exited {p.returncode}; see {log}")
    return {"chrom": c, "seconds": round(time.time() - t0, 1), "log": str(log)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", type=int, nargs="+", default=list(range(1, 23)))
    ap.add_argument("--dry-run", action="store_true",
                    help="print the schedule and the memory arithmetic, run nothing")
    ap.add_argument("--ram-fraction", type=float, default=0.7)
    a = ap.parse_args(argv)

    budget = int(R.available_ram_bytes() * a.ram_fraction)
    batches = bin_pack(a.chroms, budget)

    print(f"memory budget {budget / 1e9:.0f} GB "
          f"({a.ram_fraction:.0%} of {R.available_ram_bytes() / 1e9:.0f} GB available)")
    print(f"RSS estimate  {RSS_PER_GFA_BYTE:.1f} x gzipped GFA size, from chr21's "
          f"measured 9.48 GB peak\n")
    for i, b in enumerate(batches, 1):
        tot = sum(estimate_rss(c) for c in b)
        names = ", ".join(f"chr{c} ({estimate_rss(c) / 1e9:.0f}G)" for c in b)
        print(f"  batch {i}: {len(b)} wide, {tot / 1e9:.0f} GB est — {names}")
    print()
    if a.dry_run:
        print("dry run: nothing executed")
        return 0

    results = []
    for i, b in enumerate(batches, 1):
        plan = R.plan_workers(max(estimate_rss(c) for c in b), n_jobs=len(b),
                              ram_fraction=a.ram_fraction)
        print(f"--- batch {i}/{len(batches)}: {plan}", flush=True)
        results.extend(R.run_jobs(_build_one, b, workers=min(plan.workers, len(b)),
                                  label=f"batch-{i} chromosomes"))

    # --- completeness, asserted ---
    print("\n=== completeness ===")
    bad: list[str] = []
    for c in a.chroms:
        d = REPO / "data" / "cohort" / f"chr{c}"
        for f in EXPECTED:
            if not (d / f).exists() or (d / f).stat().st_size == 0:
                bad.append(f"chr{c}: missing or empty {f}")
    for r in sorted(results, key=lambda r: r["chrom"]):
        print(f"  chr{r['chrom']:<2} built in {r['seconds']:>7.1f} s")
    if bad:
        print(f"\nINCOMPLETE — {len(bad)} problem(s):")
        for b in bad[:20]:
            print(f"  {b}")
        print("\nDo NOT quote a genome-wide number from this tree.")
        return 1
    print(f"\nall {len(a.chroms)} chromosome(s) built and complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
