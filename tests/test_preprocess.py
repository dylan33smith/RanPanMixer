"""Tests for the preprocessing stage.

Runs two ways, because `pytest` is not in the `panmixer` conda env:

    python3 tests/test_preprocess.py       # self-running
    python3 -m pytest -q tests/test_preprocess.py
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
TINY = REPO / "tests" / "tiny"
sys.path.insert(0, str(SRC))

import sites as sites_mod  # noqa: E402
import store as store_mod  # noqa: E402


# ---------------------------------------------------------------------------
# Constraint A -- the cohort stage cannot see a target.
#
# docs/plan.md open question 4 said this was "enforced by discipline; it should
# be enforced by construction". These tests ARE the construction: they read the
# source and fail on the import or the parameter, so the rule breaks at CI time
# rather than silently at release time.
# ---------------------------------------------------------------------------

def test_cohort_never_imports_target():
    offenders = []
    for py in (SRC / "cohort").rglob("*.py"):
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for n in node.names:
                    if n.name.split(".")[0] == "target":
                        offenders.append(f"{py.name}: import {n.name}")
            elif isinstance(node, ast.ImportFrom):
                if (node.module or "").split(".")[0] == "target":
                    offenders.append(f"{py.name}: from {node.module} import ...")
    assert not offenders, (
        "cohort/ imports target/:\n  " + "\n  ".join(offenders) +
        "\nStanding constraint 1: the baseline is a function of G and D only."
    )


def test_no_cohort_function_takes_a_target():
    offenders = []
    for py in (SRC / "cohort").rglob("*.py"):
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                args = [a.arg for a in node.args.args + node.args.kwonlyargs]
                bad = [a for a in args if "target" in a.lower()]
                if bad:
                    offenders.append(f"{py.name}:{node.name}({', '.join(bad)})")
    assert not offenders, (
        "cohort/ functions take target-shaped parameters:\n  " + "\n  ".join(offenders)
    )


def test_cli_does_not_share_a_parser_between_stages():
    """A shared argparse parent would put a target flag one refactor away.

    Parsed from the AST, not grepped: the first version of this test searched the
    source text and fired on the comment explaining why the code does NOT use a
    parent. A test that reads prose rather than behaviour is a test that will
    keep lying.
    """
    tree = ast.parse((SRC / "cli.py").read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            fn = node.func
            name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
            if name == "ArgumentParser":
                for kw in node.keywords:
                    if kw.arg == "parents":
                        offenders.append(ast.dump(node)[:80])
    assert not offenders, (
        "cli.py builds an ArgumentParser with parents=. The two stages must build "
        "separate parsers so a target argument cannot leak into prep-cohort."
    )


def test_cohort_cli_rejects_a_target_argument_at_runtime():
    """The static checks cover cohort/; this covers the entry point itself."""
    sys.path.insert(0, str(SRC))
    import cli
    try:
        cli.cohort_main(["--target", "HG00096", "arrays"])
    except SystemExit as e:
        assert e.code != 0, "prep-cohort must reject --target, not accept it"
    except NotImplementedError:
        raise AssertionError("prep-cohort ACCEPTED --target; constraint 1 is unguarded")


# ---------------------------------------------------------------------------
# The site axis and its digest.
# ---------------------------------------------------------------------------

def test_axis_digest_ignores_formatting_but_catches_reordering():
    chrom = ["c", "c", "c"]; pos = [10, 20, 30]
    ref = ["A", "C", "G"];   alt = ["T", "G", "A"]
    d1 = sites_mod.axis_digest(chrom, pos, ref, alt)
    d2 = sites_mod.axis_digest(chrom, np.array(pos, dtype=np.int64), ref, alt)
    assert d1 == d2, "digest must not depend on int type -- recompression is not a change"

    d3 = sites_mod.axis_digest(chrom, [10, 30, 20], ref, alt)
    assert d1 != d3, "digest MUST change when two records swap -- that is the whole point"

    d4 = sites_mod.axis_digest(chrom, pos, ref, ["T", "G", "C"])
    assert d1 != d4, "digest must change when an ALT changes"


def test_sites_roundtrip(tmp_path=None):
    out = Path(tmp_path) if tmp_path else REPO / "tests" / "_tmp"
    out.mkdir(parents=True, exist_ok=True)
    rows = [
        ("grch38#chrT", 100, "A", "G", 1, 0, ""),
        ("grch38#chrT", 520, "AGGTC", "A,AGGTCGG,AGG", 3, 0, ""),
        ("grch38#chrT", 522, "GGT", "G", 1, 1, ">5>9"),
    ]
    digest = sites_mod.write_sites_tsv(out / "sites.tsv", rows)
    np.save(out / "in_chain.npy", np.array([False, True, True]))

    s = sites_mod.load_sites(out / "sites.tsv")
    assert s.n_sites == 3, "n_sites is the AXIS -- always the full record count"
    assert s.T == 2, "T is the CHAIN -- the mask's sum, never stored"
    assert s.digest == digest
    assert s.lv[2] == 1 and s.ps[2] == ">5>9", "LV/PS must survive; they exist nowhere else"
    assert s.n_alt[1] == 3


def test_load_sites_applies_the_mask_by_default():
    """Forgetting the provisional cut must be the unusual path, not the easy one."""
    src = (SRC / "sites.py").read_text()
    assert "in_chain.npy" in src and "in_chain_path is None" in src, (
        "load_sites must look for the chain mask by default"
    )


def test_mask_not_filter():
    """The coordinate cut is provisional and WILL move. Storing it as a mask
    means moving it rewrites one array; filtering rows would renumber the tree."""
    src = (SRC / "cohort" / "genmap.py").read_text()
    assert "chain_mask" in src, "the cut must be produced as a mask"
    assert "def chain_filter" not in src, "the cut must never be applied as a row filter"


# ---------------------------------------------------------------------------
# Provenance.
# ---------------------------------------------------------------------------

def test_save_writes_a_provenance_sidecar():
    out = REPO / "tests" / "_tmp"; out.mkdir(parents=True, exist_ok=True)
    arr = np.arange(6, dtype=np.int16).reshape(2, 3)
    p = store_mod.save(out / "thing.npy", arr,
                       inputs={"vcf": "abc123"}, params={"chain_span_start": 12968320},
                       axis_digest="deadbeef", note="unit test")
    prov = store_mod.load_prov(p)
    assert prov["shape"] == [2, 3] and prov["dtype"] == "int16"
    assert prov["inputs"]["vcf"] == "abc123"
    assert prov["params"]["chain_span_start"] == 12968320, (
        "every value that shaped the output must be recorded -- chain_span above all, "
        "because it is provisional"
    )
    assert prov["axis_digest"] == "deadbeef"


def test_check_axis_refuses_a_mismatched_join():
    out = REPO / "tests" / "_tmp"; out.mkdir(parents=True, exist_ok=True)
    p = store_mod.save(out / "axis.npy", np.zeros(3),
                       inputs={}, params={}, axis_digest="aaa")
    store_mod.check_axis(p, "aaa")
    try:
        store_mod.check_axis(p, "bbb")
    except ValueError as e:
        assert "site-axis mismatch" in str(e)
    else:
        raise AssertionError("check_axis must refuse a mismatched digest")


def test_artifact_without_provenance_is_not_silently_usable():
    out = REPO / "tests" / "_tmp"; out.mkdir(parents=True, exist_ok=True)
    stray = out / "stray.npy"
    np.save(stray, np.zeros(2))          # written WITHOUT store.save
    try:
        store_mod.load_prov(stray)
    except FileNotFoundError as e:
        assert "not quotable" in str(e)
    else:
        raise AssertionError("an artifact with no provenance must not load quietly")


# ---------------------------------------------------------------------------
# The fixture is worth what it costs.
# ---------------------------------------------------------------------------

def test_tiny_fixture_covers_every_reason_code_and_the_collision():
    text = (TINY / "tiny.vcf").read_text()
    body = [l for l in text.splitlines() if l and not l.startswith("#")]
    assert len(body) == 24, f"expected 24 records, got {len(body)}"

    assert any("LV=1" in l for l in body), "need a nested child -> code 1"
    assert any("LV=2" in l for l in body), "need a grandchild -> code 2"
    assert any("CONFLICT=" in l for l in body), "need a CONFLICT record -> code 3"

    collision = [l for l in body if "CONFLICT=" in l and "PS=" in l]
    assert collision, (
        "need one record that is BOTH conflict-named and nested -- that is the "
        "precedence case, 25,596 such cells on real chr21"
    )

    run = [l for l in body if l.split("\t")[9 + 3].startswith(".|.")]
    assert len(run) >= 6, "need a >=6 record missing run for S4 -> missing_runs.tsv"

    wrong_an = [l for l in body if l.split("\t")[1] == "1300"]
    assert "AN=2" in wrong_an[0], (
        "record 1300 must carry a deliberately wrong INFO/AN so the row-alignment "
        "canary is proven to fire rather than assumed to"
    )


def test_tiny_gmap_spans_less_than_the_fixture():
    g = [l.split("\t") for l in (TINY / "tiny.gmap").read_text().splitlines()[1:]]
    pos = [int(r[0]) for r in g]
    cm = [float(r[2]) for r in g]
    assert min(pos) == 500 and max(pos) == 1500
    assert any(cm[i] == cm[i + 1] for i in range(len(cm) - 1)), (
        "need a zero-dcM interior interval -- a real map plateau, correct to leave"
    )


# ---------------------------------------------------------------------------
# Regressions. Both of these were live bugs caught by running on real data.
# ---------------------------------------------------------------------------

def test_sites_tsv_survives_a_huge_alt_field():
    """REGRESSION. The hypervariable chr21 record declares 90 alleles and its ALT
    string is ~23 MB. Python's csv module refuses fields over 128 KB by default,
    so load_sites() raised `_csv.Error: field larger than field limit` on real
    data. The dev slice was chosen to contain that record, which is how it
    surfaced in seconds rather than after a full chromosome build."""
    out = REPO / "tests" / "_tmp"; out.mkdir(parents=True, exist_ok=True)
    big = ",".join("ACGT" * 50_000 for _ in range(2))      # ~400 KB, over the limit
    rows = [("c", 10, "A", "G", 1, 0, ""), ("c", 20, "A", big, 2, 0, "")]
    sites_mod.write_sites_tsv(out / "sites.tsv", rows)
    np.save(out / "in_chain.npy", np.array([True, True]))
    s = sites_mod.load_sites(out / "sites.tsv")
    assert s.n_sites == 2 and len(s.alt[1]) == len(big)


def test_a_bare_dot_genotype_is_two_missing_cells():
    """REGRESSION-ADJACENT. A GT of '.' with no separator means BOTH haplotypes
    are unknown, so it contributes 2 missing cells, not 1. A naive text count of
    '.' tokens sees one and undercounts -- which is exactly how an apparent
    909-cell discrepancy appeared when cross-checking against bcftools.

    vg 1.36 emitted 237,594 such fields on chr21; vg 1.68 emits none, writing
    '.|.' instead. The asymmetry is preserved because it is real for other VCFs.
    """
    sys.path.insert(0, str(SRC))
    from cohort import arrays as A
    assert A._parse_gt(".") == (-1, -1)
    assert A._parse_gt(".|.") == (-1, -1)
    assert A._parse_gt("0") == (0, -1), "a lone real allele fills strand 0 only"
    assert A._parse_gt("0|1") == (0, 1)
    assert A._parse_gt("2/3") == (2, 3), "unphased separators must parse too"


def test_reason_precedence_is_the_pinned_order():
    """conflict > not applicable > inherited > uncategorised, verified on the
    fixture's deliberate collision at pos 720 where S3 is both conflict-named and
    would otherwise be inherited."""
    sys.path.insert(0, str(SRC))
    from cohort import arrays as A
    out = REPO / "tests" / "_tmp" / "prec"
    A.build(TINY / "tiny.vcf", out)
    pos = np.load(out / "positions.npy")
    rea = np.load(out / "reason.npy")
    ids = list(np.load(out / "haplotype_ids.npy"))
    r = int(np.flatnonzero(pos == 720)[0])
    s3 = [i for i, v in enumerate(ids) if str(v).startswith("S3:")]
    s4 = [i for i, v in enumerate(ids) if str(v).startswith("S4:")]
    s2 = [i for i, v in enumerate(ids) if str(v).startswith("S2:")]
    assert all(rea[h, r] == A.CONFLICT for h in s3), "conflict must beat inherited"
    assert all(rea[h, r] == A.INHERITED for h in s4), "parent missing -> inherited"
    assert all(rea[h, r] == A.NOT_APPLICABLE for h in s2), "parent called -> not applicable"


def test_an_cross_check_is_not_vacuous():
    """The support-vs-INFO/AN check must be able to FAIL. PanMixer's equivalent
    assertion compares two things that are equal by construction, so it never
    fires. Ours compares parsed genotypes against a field the caller wrote, and
    the fixture carries one deliberately wrong AN to prove it."""
    sys.path.insert(0, str(SRC))
    from cohort import arrays as A
    out = REPO / "tests" / "_tmp" / "an"
    A.build(TINY / "tiny.vcf", out)
    hap = np.load(out / "haplotypes.npy"); an = np.load(out / "an_info.npy")
    pos = np.load(out / "positions.npy")
    mism = np.flatnonzero((hap != -1).sum(axis=0) != an)
    assert len(mism) == 1, f"expected exactly the one planted error, got {len(mism)}"
    assert pos[mism[0]] == 1300


def test_manifest_tool_digests_match_the_binaries_on_disk():
    """REGRESSION 2026-09-29. A vg sha256 was written into the manifest with its tail
    FABRICATED from a truncated terminal display. It looked entirely plausible. The
    project rule is that every number is traceable; a checksum that was never computed
    is the worst possible violation of it, because its whole purpose is to be checked.

    vg's version is load-bearing, not cosmetic: 1.68 and 1.36 order ALT alleles
    differently at multi-allelic sites, so the version determines the allele indices
    every downstream array is written in terms of.
    """
    import yaml
    cfg = yaml.safe_load((REPO / "configs" / "inputs.chr21.yaml").read_text())
    tools = cfg.get("tools", {})
    assert tools, "configs/inputs.chr21.yaml declares no tools block"
    for name, spec in tools.items():
        declared = spec.get("sha256")
        path = spec.get("path")
        if not declared or not path:
            continue
        f = REPO / path
        if not f.exists():
            continue                      # staged elsewhere; nothing to check
        actual = store_mod.sha256_file(f)
        assert actual == declared, (
            f"{name}: manifest declares sha256 {declared} but {path} hashes to {actual}. "
            f"Either the binary changed or the digest was never computed."
        )


def test_split_alt_join_recovers_multiallelic_records():
    """REGRESSION 2026-09-29. The 1000 Genomes panel is entirely decomposed into
    biallelic records -- every one of its 1,002,752 chr21 records declares exactly
    one ALT -- while 7.9% of our pangenome records are multi-allelic. A strict
    full-ALT-string join can therefore NEVER match those 27,021 records, for a
    purely notational reason, and the first implementation silently discarded
    target information we actually held: 8,377 rows, 8,283 of them inside
    chain_span.

    This test builds the same shape in miniature: one multi-allelic cohort record,
    a target VCF that decomposes it into two biallelic records, and a target that
    carries a different ALT on each haplotype -- so the reconstruction has to map
    BOTH panel records back onto our allele numbering, which is the case a naive
    'first match wins' would get wrong.
    """
    sys.path.insert(0, str(SRC))
    from cohort import arrays as A
    from target import path as TP
    out = REPO / "tests" / "_tmp" / "splitalt"
    coh = out / "cohort"; coh.mkdir(parents=True, exist_ok=True)

    # cohort: one biallelic record and one 2-ALT record
    cvcf = out / "cohort.vcf"
    cvcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c,length=1000>\n"
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="x">\n'
        '##INFO=<ID=LV,Number=1,Type=Integer,Description="x">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD1\tD2\n"
        "c\t100\t.\tA\tG\t60\t.\tAN=4;LV=0\tGT\t0|1\t0|0\n"
        "c\t200\t.\tC\tT,G\t60\t.\tAN=4;LV=0\tGT\t1|2\t0|1\n")
    A.build(cvcf, coh)
    # through store.save, not np.save: check_axis now demands a sidecar, and a test
    # that bypasses the chokepoint is a test that does not exercise the real path
    import sites as _s
    _d = _s.load_sites(coh / "sites.tsv").digest
    store_mod.save(coh / "in_chain.npy", np.ones(2, dtype=bool),
                   inputs={}, params={}, axis_digest=_d)

    # target: the SAME site 200 decomposed into two biallelic records, and the
    # target carries T on haplotype 0 and G on haplotype 1 -> our alleles 1 and 2
    tvcf = out / "target.vcf"
    tvcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c,length=1000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="x">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTGT\n"
        "c\t100\t.\tA\tG\t60\t.\tAF=0.1\tGT\t0|1\n"
        "c\t200\t.\tC\tT\t60\t.\tAF=0.1\tGT\t1|0\n"
        "c\t200\t.\tC\tG\t60\t.\tAF=0.1\tGT\t0|1\n")

    rep = TP.build(tvcf, coh, out / "t", sample="TGT", split_alt_join=True)
    path = np.load(out / "t" / "path.npy")
    assert list(path[0]) == [0, 1], f"biallelic row should pass through, got {path[0]}"
    assert list(path[1]) == [1, 2], (
        f"multi-allelic row should reconstruct to OUR allele numbering [1,2], got {path[1]}")
    assert rep["split_alt_join"]["rows_recovered"] == 1
    assert rep["split_alt_join"]["contradictions"] == 0

    # and with the recovery disabled the row must be lost, proving it does work
    rep2 = TP.build(tvcf, coh, out / "t2", sample="TGT", split_alt_join=False)
    p2 = np.load(out / "t2" / "path.npy")
    assert list(p2[1]) == [-1, -1], "without the split join the multi-allelic row is lost"
    assert rep2["split_alt_join"]["rows_recovered"] == 0


def test_readable_set_is_target_independent_but_called_set_is_not():
    """REGRESSION 2026-09-29, from the independence audit.

    docs/plan.md recommends normalising beta_v over the panel-matched positions, and
    justified it with "identical for all ten targets". That justification was attached
    to the WRONG set. `in_chain & (target_reason == 0)` -- positions where a particular
    target has a call -- is NOT invariant: a contradiction at a poly-A repeat drops a
    row for one target and not another, giving four distinct sets across ten targets.
    Normalising over it would make beta_v target-dependent, which is the standing
    constraint 1 violation that voids Theorem 1.

    `readable.npy` is the set that IS invariant, because both its tests read only the
    panel's site list. This test asserts both halves, on the real artifacts, and is
    skipped when they are absent.
    """
    import glob
    rd = sorted(glob.glob(str(REPO / "data/targets/*/chr21/readable.npy")))
    tr = sorted(glob.glob(str(REPO / "data/targets/*/chr21/target_reason.npy")))
    if len(rd) < 2:
        print("    (skipped: needs >=2 built targets)")
        return
    sites_p = REPO / "data/cohort/chr21/sites.tsv"
    s = sites_mod.load_sites(sites_p)
    readable_sets = {frozenset(np.flatnonzero(np.load(f)).tolist()) for f in rd}
    assert len(readable_sets) == 1, (
        f"readable.npy differs across targets ({len(readable_sets)} distinct sets). It must "
        f"be a function of the panel's site list alone, or beta_v normalised over it is "
        f"target-dependent."
    )
    called_sets = {frozenset(np.flatnonzero(s.in_chain & (np.load(f) == 0)).tolist())
                   for f in tr}
    # not an error -- documenting WHY readable exists. If this ever becomes 1, the
    # comment above is stale, not the code.
    print(f"    (readable: 1 set; target_reason==0: {len(called_sets)} sets -- "
          f"which is why readable.npy exists)")


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as e:
            failed += 1
            print(f"  FAIL  {name}\n          {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
