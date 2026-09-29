"""Two entry points, and deliberately two parsers.

`prep-cohort` and `prep-target` do NOT share an argparse parent. A shared
`parents=[common]` would put a target argument one refactor away from the cohort
stage, and standing constraint 1 is the rule this project says is most likely to
be broken by accident. The duplication is the point.

Every command takes `--root`, which relocates the entire output tree. That is how
the 2 Mb development slice is built and consumed by the identical code path: there
is no test mode to drift out of sync with the real one.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_ROOT = Path("data")


def _chrom_and_root(p: argparse.ArgumentParser) -> None:
    p.add_argument("--chrom", default="chr21", help="chromosome to process")
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                   help="output tree root; use data/dev for the slice fixture")
    p.add_argument("--config", type=Path, default=Path("configs/inputs.chr21.yaml"))
    p.add_argument("--force", action="store_true", help="rebuild even if outputs exist")


def cohort_main(argv=None) -> int:
    """Build the cohort artifacts. Takes no target, and must never learn how."""
    p = argparse.ArgumentParser(
        prog="prep-cohort",
        description="Graph + genetic map -> the cohort artifacts. No target input exists.",
    )
    _chrom_and_root(p)
    p.add_argument("step", choices=["deconstruct", "cohort-vcf", "arrays", "genmap", "all"],
                   help="run one step, or all of them in order")
    a = p.parse_args(argv)

    # The guard is structural, but assert it anyway: a future edit that adds a
    # target flag here should fail immediately and visibly.
    assert not any("target" in k for k in vars(a)), (
        "prep-cohort has acquired a target-shaped argument. That breaks standing "
        "constraint 1: the baseline must be a function of G and D only."
    )
    raise NotImplementedError(f"step {a.step!r} not implemented yet")


def target_main(argv=None) -> int:
    """Map one target onto the site axis. Consumes cohort artifacts read-only."""
    p = argparse.ArgumentParser(
        prog="prep-target",
        description="A phased GRCh38 target VCF -> its path over the site axis.",
    )
    _chrom_and_root(p)
    p.add_argument("--target", required=True,
                   help="sample id; must be in configs/clean_pool.chr21.txt")
    a = p.parse_args(argv)
    raise NotImplementedError(f"target {a.target!r} not implemented yet")


if __name__ == "__main__":
    sys.exit(cohort_main())
