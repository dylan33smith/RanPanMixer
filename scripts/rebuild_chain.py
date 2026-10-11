#!/usr/bin/env python
"""Re-run the `genmap` and `manifest` steps on every autosome, for decision K.

⚠ WHY ONLY THESE TWO STEPS. Decision K changes `in_chain.npy` and nothing upstream of
it: `deconstruct`, `cohort-vcf`, `arrays` and `support` are untouched, and re-running
them would redo hours of work for no change. `genmap` is forced because its output
already exists; `manifest` is not in the resume map and always re-runs.

⚠ AND IT ASSERTS WHAT CHANGED. A silent rebuild of the chain is the most dangerous
edit available here -- every denominator in the project is a count of chain positions,
and `beta_share` is allocated by variant count under decision J, so a changed T changes
the budget split on every other chromosome too.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402
import runner as R                                        # noqa: E402


def _chain_size(c: int) -> int:
    f = REPO / "data" / "cohort" / f"chr{c}" / "in_chain.npy"
    return int(np.load(f).astype(bool).sum()) if f.exists() else -1


def _one(c: int) -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    cfg = REPO / "configs" / f"inputs.chr{c}.yaml"
    before = _chain_size(c)
    t0 = time.time()
    log = REPO / "logs" / "cohorts" / f"chr{c}.chain.log"
    with open(log, "w") as fh:
        for step, extra in (("genmap", ["--force"]), ("manifest", [])):
            p = subprocess.run(
                [sys.executable, "-c", "import cli, sys; sys.exit(cli.cohort_main())",
                 step, "--chrom", f"chr{c}", "--config", str(cfg), *extra],
                stdout=fh, stderr=subprocess.STDOUT, cwd=REPO, env=env)
            if p.returncode != 0:
                raise RuntimeError(f"chr{c}: {step} exited {p.returncode}; see {log}")
    man = json.loads((REPO / "data" / "cohort" / f"chr{c}" / "manifest.json").read_text())
    blocks = (man.get("artifacts", {}).get("in_chain.npy", {})
                 .get("params", {}).get("untraversable_blocks") or {})
    return {"chrom": c, "before": before, "after": _chain_size(c),
            "blocks": blocks.get("untraversable_blocks_dropped", 0),
            "dropped": blocks.get("positions_dropped", 0),
            "seconds": round(time.time() - t0, 1)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", type=int, nargs="+", default=list(range(1, 23)))
    a = ap.parse_args(argv)

    # the allowed matrix is the memory cost: K x T bools, plus positions and cm
    per_job = 88 * 2_000_000 + 3 * 8 * 2_000_000
    plan = R.plan_workers(per_job, n_jobs=len(a.chroms))
    print(f"{plan}\n")
    res = R.run_jobs(_one, a.chroms, workers=plan.workers, label="chromosomes",
                     progress_every=4)

    print("\n| chr | T before | T after | blocks dropped | positions dropped | % |")
    print("|---|---|---|---|---|---|")
    tb = ta = td = nb = 0
    for r in sorted(res, key=lambda r: r["chrom"]):
        pc = 100.0 * r["dropped"] / r["before"] if r["before"] > 0 else 0.0
        tb += r["before"]; ta += r["after"]; td += r["dropped"]; nb += r["blocks"]
        print(f"| {r['chrom']} | {r['before']:,} | {r['after']:,} | {r['blocks']} | "
              f"{r['dropped']:,} | {pc:.3f}% |")
    print(f"| **total** | **{tb:,}** | **{ta:,}** | **{nb}** | **{td:,}** | "
          f"**{100.0 * td / tb:.3f}%** |")

    # ⚠ ASSERT THE ARITHMETIC CLOSES. before - dropped must equal after on every
    # chromosome; if it does not, something other than decision K moved the chain.
    bad = [r for r in res if r["before"] - r["dropped"] != r["after"]]
    if bad:
        print(f"\n⚠⚠ {len(bad)} chromosome(s) where (before - dropped) != after, so "
              f"something OTHER than decision K changed the chain:")
        for r in bad:
            print(f"  chr{r['chrom']}: {r['before']:,} - {r['dropped']:,} "
                  f"!= {r['after']:,}")
        return 1
    print("\n(before - dropped) == after on every chromosome")
    return 0


if __name__ == "__main__":
    sys.exit(main())
