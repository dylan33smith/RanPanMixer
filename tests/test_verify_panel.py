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


# ---------------------------------------------------------------------------
# ⚠⚠ THE CHECK THAT EVERY OTHER CHECK MISSED. On 2026-10-11 chr10 arrived with a CRC32
# mismatch in the bgzf block at offset 393,732,438. It passed the index check, the
# record count, the sample set, the reference-build comparison AND the truncation seek
# -- because the corruption is mid-file and a seek to the last 3 Mb never touches it.
#
# ⚠ The sha256 check passed too, and that is the deeper point: `get_panel.sh` computes
# the sidecar FROM the downloaded bytes, so it matches whatever arrived. A checksum you
# compute yourself after download verifies STORAGE, not TRANSFER. Only a structural
# check -- every bgzf block's own CRC -- can fail on its own terms.
# ---------------------------------------------------------------------------

def test_bgzf_integrity_is_checked_under_deep(V) -> None:
    """⚠ Structural: `--deep` must run `bgzip -t`, not only hash the file. Hashing
    against our own post-download sidecar cannot detect a corrupt transfer."""
    src = (REPO / "scripts" / "verify_panel.py").read_text()
    assert '"bgzip", "-t"' in src, (
        "--deep does not run `bgzip -t`. Without it a mid-file CRC failure passes "
        "every check in this script, which is how chr10 reached the pipeline.")
    assert "verifies STORAGE, not TRANSFER" in src, (
        "the reason the sha256 is insufficient must stay written down, or the next "
        "person will treat the checksum as transfer verification again")


def test_a_corrupt_bgzf_file_is_detected(tmp_path: Path) -> None:
    """⚠ THE POWER TEST, on a file corrupted the same way chr10 was: a byte flipped in
    the middle of a bgzf member, leaving length and structure intact."""
    import gzip
    import subprocess
    good = tmp_path / "good.vcf.gz"
    # a real multi-block bgzf file via bgzip, so the CRCs are genuine
    payload = b"".join(b"chr1\t%d\t.\tA\tG\t.\t.\t.\n" % i for i in range(60000))
    raw = tmp_path / "x.vcf"
    raw.write_bytes(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    + payload)
    subprocess.run(["bgzip", "-c", str(raw)], stdout=good.open("wb"), check=True)
    assert subprocess.run(["bgzip", "-t", str(good)]).returncode == 0, "control failed"

    bad = tmp_path / "bad.vcf.gz"
    data = bytearray(good.read_bytes())
    mid = len(data) // 2
    data[mid] ^= 0xFF                      # flip one byte, mid-file
    bad.write_bytes(bytes(data))
    assert subprocess.run(["bgzip", "-t", str(bad)],
                          capture_output=True).returncode != 0, (
        "bgzip -t did not detect a flipped mid-file byte, so it cannot be relied on "
        "as the transfer check")
    # ⚠ AND THE POINT: a sha256 of the corrupt file is perfectly self-consistent
    import hashlib
    h = hashlib.sha256(bad.read_bytes()).hexdigest()
    sidecar = tmp_path / "bad.vcf.gz.sha256"
    sidecar.write_text(h + "\n")           # exactly what get_panel.sh used to do
    assert hashlib.sha256(bad.read_bytes()).hexdigest() == sidecar.read_text().strip(), (
        "this assertion is the finding: hashing a corrupt download produces a sidecar "
        "that matches it, so the checksum can never reveal the corruption")


def test_the_fetcher_verifies_before_accepting(V) -> None:
    """A corrupt `.part` must be DISCARDED, not renamed into place and not resumed --
    `curl -C -` would see a complete file and skip it forever."""
    src = (REPO / "scripts" / "get_panel.sh").read_text()
    assert "bgzip -t" in src, "get_panel.sh accepts a download without testing it"
    i, j = src.index("bgzip -t"), src.index('mv "$f.part" "$f"')
    assert i < j, "the integrity test must run BEFORE the file is moved into place"
    assert 'rm -f "$f.part"' in src, (
        "a failed .part must be removed, or curl -C - will resume a corrupt file")
