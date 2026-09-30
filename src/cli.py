"""Two entry points, and deliberately two parsers.

`prep-cohort` and `prep-target` do NOT share an argparse parent. A shared
`parents=[common]` would put a target argument one refactor away from the cohort
stage, and standing constraint 1 is the rule this project says is most likely to
be broken by accident. The duplication is the point.

Every command takes `--root`, which relocates the entire output tree. That is how
the 2 Mb development slice is built and consumed by the identical code path:
there is no test mode to drift out of sync with the real one.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import store

DEFAULT_ROOT = Path("data")
COHORT_STEPS = ["deconstruct", "cohort-vcf", "arrays", "support", "genmap", "manifest"]


def _cfg(path: Path) -> dict:
    import yaml  # lazy: keeps the module importable where pyyaml is absent
    return yaml.safe_load(Path(path).read_text())


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                   help="output tree root; use data/dev for the slice fixture")
    p.add_argument("--config", type=Path, default=Path("configs/inputs.chr21.yaml"))
    p.add_argument("--force", action="store_true")


def cohort_main(argv=None) -> int:
    """Build the cohort artifacts. Takes no target, and must never learn how."""
    p = argparse.ArgumentParser(
        prog="prep-cohort",
        description="Graph + genetic map -> the cohort artifacts. No target input exists.")
    _common(p)
    p.add_argument("step", choices=COHORT_STEPS + ["all"])
    p.add_argument("--from-vcf", type=Path,
                   help="skip deconstruct and start from this cohort VCF")
    a = p.parse_args(argv)

    assert not any("target" in k for k in vars(a)), (
        "prep-cohort has acquired a target-shaped argument. That breaks standing "
        "constraint 1: the baseline must be a function of G and D only.")

    from cohort import arrays, genmap, support, vcf  # noqa: E402

    cfg = _cfg(a.config)
    temp = a.root / "temp" / a.chrom
    out = a.root / "cohort" / a.chrom
    out.mkdir(parents=True, exist_ok=True)
    steps = COHORT_STEPS if a.step == "all" else [a.step]
    summary: dict = {}

    for step in steps:
        done = {
            "arrays": out / "haplotypes.npy",
            "support": out / "support.npy",
            "genmap": out / "in_chain.npy",
            "cohort-vcf": temp / "cohort44.vcf.gz",
        }.get(step)
        if done and done.exists() and not a.force:
            print(f"  skip   {step}  ({done.name} exists; --force to rebuild)")
            continue

        print(f"  run    {step}")
        if step == "deconstruct":
            summary[step] = vcf.deconstruct(Path(cfg["inputs"]["graph"]["path"]),
                                            temp / "deconstruct.vcf", cfg=cfg)
        elif step == "cohort-vcf":
            src = a.from_vcf or (temp / "deconstruct.vcf")
            summary[step] = vcf.drop_chm13(src, temp / "cohort44.vcf.gz")
        elif step == "arrays":
            summary[step] = arrays.build(temp / "cohort44.vcf.gz", out, chrom=a.chrom)
        elif step == "support":
            summary[step] = support.build(out)
        elif step == "genmap":
            span = cfg["chain_span"]
            summary[step] = genmap.build(out / "positions.npy",
                                         Path(cfg["inputs"]["genetic_map"]["path"]), out,
                                         start=span["start"], end=span["end"])
        elif step == "manifest":
            summary[step] = _manifest(out, cfg, a.chrom)
        print(json.dumps(summary.get(step, {}), indent=2)[:1200])
    return 0


def _manifest(out: Path, cfg: dict, chrom: str) -> dict:
    """Roll every artifact's provenance into one file."""
    arts = {}
    for f in sorted(out.iterdir()):
        if f.name.endswith(store.PROV_SUFFIX) or f.name == "manifest.json":
            continue
        try:
            arts[f.name] = store.load_prov(f)
        except FileNotFoundError:
            arts[f.name] = {"ERROR": "no provenance sidecar -- not written by store.save()"}
    digests = {v.get("axis_digest") for v in arts.values() if v.get("axis_digest")}
    man = {
        "chrom": chrom,
        "chain_span": cfg["chain_span"],
        "tools": cfg.get("tools", {}),
        "artifacts": arts,
        "axis_digests_seen": sorted(digests),
        "axis_consistent": len(digests) <= 1,
    }
    (out / "manifest.json").write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
    return {"artifacts": len(arts), "axis_consistent": man["axis_consistent"],
            "axis_digests_seen": sorted(digests)}


def target_main(argv=None) -> int:
    """Map one target onto the site axis. Consumes cohort artifacts read-only."""
    p = argparse.ArgumentParser(
        prog="prep-target",
        description="A phased GRCh38 target VCF -> its path over the site axis.")
    _common(p)
    p.add_argument("--target", required=True, help="sample id")
    p.add_argument("--target-vcf", type=Path,
                   help="defaults to data/targets/<id>/<chrom>/target.<chrom>.vcf.gz")
    p.add_argument("--absent-policy", default="missing", choices=["missing", "reference"],
                   help="what a site absent from the target VCF means")
    p.add_argument("--no-split-alt-join", action="store_true",
                   help="disable recovery of multi-allelic records from the panel's "
                        "biallelic decomposition (on by default; disabling it discards "
                        "target information we hold)")
    a = p.parse_args(argv)

    from target import path as tpath  # noqa: E402

    tv = a.target_vcf or (a.root / "targets" / a.target / a.chrom /
                          f"target.{a.chrom}.vcf.gz")
    rep = tpath.build(tv, a.root / "cohort" / a.chrom,
                      a.root / "targets" / a.target / a.chrom,
                      sample=a.target, absent_policy=a.absent_policy,
                      split_alt_join=not a.no_split_alt_join)
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(cohort_main())
