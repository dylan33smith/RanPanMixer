"""Does every check in `scripts/verify_panel.py` actually FIRE?

⚠ The load-bearing one is the REFERENCE BUILD check, and it is load-bearing because a
wrong-build panel does not fail loudly -- it produces a low (POS, REF, ALT) match rate
that reads as biology. Decision F exists because PanMixer's released pipeline pairs a
GRCh37 panel with a GRCh38 pangenome and the cost is a 252x difference in match rate.
So the test below feeds the checker a GRCh37 chr22 length and requires it to refuse.

⚠ These tests need the real panel and SKIP without it, with a stated reason, rather
than reporting green on checks that never ran.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

CHROM = 22                      # the smallest panel file, so the seek is cheap

pytestmark = pytest.mark.skipif(
    not Path("/data/ds85/RanPanMixer/data/raw/1000g_30x_chr22.vcf.gz").exists()
    or not (REPO / "data" / "temp" / "chr22" / "deconstruct.vcf").exists(),
    reason="the 1000G panel and our own chr22 deconstruct output are both needed; "
           "the reference-build check compares them against each other")


@pytest.fixture()
def V():
    import verify_panel
    return verify_panel


def test_control_real_panel_chromosome_passes(V) -> None:
    """⚠ THE CONTROL. Without it, a checker that complained about every input would
    pass every mutation case below."""
    bad, info = V.check(CHROM)
    assert bad == [], bad
    assert info["samples"] == 3202
    assert info["records"] > 0


def test_grch37_length_is_refused(V, monkeypatch) -> None:
    """⚠ THE DECISION-F GATE. GRCh37 chr22 is 51,304,566 bp; GRCh38 is 50,818,468.
    Feeding the checker the GRCh37 value must produce a build mismatch."""
    monkeypatch.setattr(V, "our_contig_length", lambda c: 51_304_566)
    bad, _ = V.check(CHROM)
    assert any("REFERENCE BUILD MISMATCH" in b for b in bad), bad


def test_a_one_base_difference_is_refused(V, monkeypatch) -> None:
    """The check must be exact, not approximate. Two assemblies differing by a single
    base are still two assemblies, and a tolerance here would admit a patch release."""
    monkeypatch.setattr(V, "our_contig_length", lambda c: 50_818_469)
    bad, _ = V.check(CHROM)
    assert any("REFERENCE BUILD MISMATCH" in b for b in bad), bad


def test_absent_reference_is_reported_not_passed(V, monkeypatch) -> None:
    """⚠ If our own deconstruct output is missing, the build CANNOT be checked. That
    must be said out loud: silently passing is how an unchecked panel becomes a
    'verified' one."""
    monkeypatch.setattr(V, "our_contig_length", lambda c: None)
    bad, _ = V.check(CHROM)
    assert any("CANNOT BE CHECKED" in b for b in bad), bad


def test_truncation_is_caught(V, monkeypatch) -> None:
    """A download cut short still opens and still indexes against its own partial
    content; what gives it away is the last reachable record falling short of the
    contig. Simulated here by declaring a contig twice as long."""
    real = V.panel_contig_length(CHROM)
    monkeypatch.setattr(V, "panel_contig_length", lambda c: real * 2)
    bad, _ = V.check(CHROM)
    assert any("reaches only" in b or "NO records in the last" in b for b in bad), bad


def test_missing_checksum_is_caught(V, monkeypatch, tmp_path) -> None:
    """A panel VCF with no recorded checksum is not quotable, and re-deriving one
    means re-reading 29 GB."""
    monkeypatch.setattr(V, "panel_vcf", lambda c: tmp_path / f"1000g_30x_chr{c}.vcf.gz")
    (tmp_path / f"1000g_30x_chr{CHROM}.vcf.gz").write_bytes(b"x")
    (tmp_path / f"1000g_30x_chr{CHROM}.vcf.gz.tbi").write_bytes(b"x")
    bad, _ = V.check(CHROM)
    assert any("checksum" in b for b in bad), bad


def test_missing_index_is_caught(V, monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(V, "panel_vcf", lambda c: tmp_path / f"1000g_30x_chr{c}.vcf.gz")
    (tmp_path / f"1000g_30x_chr{CHROM}.vcf.gz").write_bytes(b"x")
    (tmp_path / f"1000g_30x_chr{CHROM}.vcf.gz.sha256").write_text("deadbeef\n")
    bad, _ = V.check(CHROM)
    assert any("index" in b for b in bad), bad


def test_leftover_part_file_fails_the_run(V, monkeypatch, tmp_path, capsys) -> None:
    """⚠ A leftover `.part` means a download that never finished, whatever the exit
    sentinel said. The per-chromosome checks cannot see it, so `main` must."""
    monkeypatch.setattr(V, "PANEL", tmp_path)
    (tmp_path / "1000g_30x_chr5.vcf.gz.part").write_bytes(b"partial")
    monkeypatch.setattr(V, "check", lambda c, deep=False: ([], {
        "chrom": c, "bytes": 1, "records": 1, "length": 1, "samples": 3202,
        "last_pos": 1, "reach_pct": 100.0}))
    rc = V.main(["--chroms", str(CHROM)])
    assert rc == 1
    assert "leftover .part" in capsys.readouterr().out


def test_mixed_sample_sets_fail_the_run(V, monkeypatch, capsys) -> None:
    """One chromosome fetched from a different release carries a different panel, and
    a per-chromosome sample set is not one cohort. Only `main` sees across them."""
    seen = {"n": 3202}

    def fake(c, deep=False):
        seen["n"] = 2504 if c == 2 else 3202      # Phase 3 had 2,504
        return [], {"chrom": c, "bytes": 1, "records": 1, "length": 1,
                    "samples": seen["n"], "last_pos": 1, "reach_pct": 100.0}

    monkeypatch.setattr(V, "check", fake)
    monkeypatch.setattr(V, "PANEL", Path("/nonexistent"))
    rc = V.main(["--chroms", "1", "2"])
    assert rc == 1
    assert "do NOT share one sample set" in capsys.readouterr().out
