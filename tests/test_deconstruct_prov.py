"""The deconstruct step must record the GFA it read BEFORE deleting it.

⚠ WHY THIS FILE EXISTS. On 2026-10-10 six chromosomes (chr10-14, chr17) failed
after vg had already spent one to two hours each, because the decompressed GFA was
unlinked three lines before its digest was taken. Nothing in the suite could see it:
every existing test either passed a plain `.gfa` (so `decompressed_here` was False
and the unlink never fired) or stubbed the provenance call away. The failure needed
the .gz path AND real provenance AND the cleanup, together.

⚠ AND THE TRACEBACK POINTED AT THE WRONG THING. It named `store.sha256_file`, so it
read as a provenance bug; the cause was three lines above, in cleanup code added for
an unrelated reason. That is why the test below asserts the ORDER by its observable
consequence -- a sidecar whose digest matches the decompressed bytes -- rather than
asserting that no exception was raised.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import store                                   # noqa: E402
from cohort import vcf as V                     # noqa: E402

GFA_TEXT = (b"H\tVN:Z:1.0\n"
            b"S\t1\tACGTACGTAC\n"
            b"S\t2\tGGTT\n"
            b"L\t1\t+\t2\t+\t0M\n"
            b"P\tgrch38#0#chr21\t1+,2+\t*\n"
            b"P\tHG00096#1#chr21\t1+\t*\n")

VCF_TEXT = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\n"
            "chr21\t5\t.\tA\tG\t.\t.\tAC=1;AN=2\tGT\t0|1\n")


@pytest.fixture()
def staged(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """A gzipped GFA and a stubbed vg that writes a VCF, so the real code path runs.

    ⚠ The GFA must be GZIPPED. With a plain `.gfa` the function sets
    `decompressed_here = False`, never unlinks, and the ordering bug is invisible --
    which is exactly why the pre-existing tests missed it.
    """
    gz = tmp_path / "chr21.hprc-v1.0-pggb.gfa.gz"
    with gzip.open(gz, "wb") as fh:
        fh.write(GFA_TEXT)

    fake = {"path": "/fake/vg", "version": "vg 1.68.0", "sha256": "f" * 64,
            "digest_checked_against_manifest": False}
    monkeypatch.setattr(V, "resolve_vg", lambda cfg, vg: Path("/fake/vg"))
    monkeypatch.setattr(V, "verify_vg", lambda resolved, cfg: dict(fake))

    def fake_run(cmd, stdout_path=None, **kw):
        # vg's contract as this code depends on it: the decompressed GFA must still
        # be on disk when vg runs, and the VCF is written before anything else.
        gfa_arg = Path(cmd[-1])
        assert gfa_arg.exists(), f"vg was handed {gfa_arg}, which does not exist"
        Path(stdout_path).write_text(VCF_TEXT)
        return {"returncode": 0}

    monkeypatch.setattr(V, "_run", fake_run)
    return gz, tmp_path


def test_decompressed_gfa_is_deleted(staged) -> None:
    """The control for the test below: cleanup must actually happen, or the
    ordering assertion would pass trivially on a function that never deletes."""
    gz, tmp = staged
    V.deconstruct(gz, tmp / "out" / "deconstruct.vcf")
    assert not (tmp / "out" / "chr21.hprc-v1.0-pggb.gfa").exists(), \
        "the decompressed copy was kept; the 4-9 GB per chromosome cleanup regressed"


def test_prov_sidecar_records_the_decompressed_gfa_digest(staged) -> None:
    """⚠ THE REGRESSION. If the unlink moves back before the digest, this raises
    FileNotFoundError inside `deconstruct` and the test fails. It also pins the
    VALUE, so computing the digest of the .gz by mistake is caught too."""
    gz, tmp = staged
    out = tmp / "out" / "deconstruct.vcf"
    V.deconstruct(gz, out)

    prov = json.loads((tmp / "out" / ("deconstruct.vcf" + store.PROV_SUFFIX)).read_text())
    import hashlib
    expect_plain = hashlib.sha256(GFA_TEXT).hexdigest()
    expect_gz = store.sha256_file(gz)
    assert prov["inputs"]["gfa"] == expect_plain, (
        "the sidecar does not carry the digest of the bytes vg actually read")
    assert prov["inputs"]["gfa"] != expect_gz, (
        "the sidecar carries the .gz digest, not the decompressed one — that silently "
        "changes the meaning of the field the first 11 chromosomes already recorded")


def test_sha256_gz_plain_equals_the_decompressed_digest(staged) -> None:
    """The property the whole recovery path rests on: the digest of the deleted
    decompressed file is recoverable from the .gz, because gzip is deterministic."""
    gz, _ = staged
    import hashlib
    assert store.sha256_gz_plain(gz) == hashlib.sha256(GFA_TEXT).hexdigest()


def test_recovery_attaches_a_sidecar_without_rerunning_vg(staged, monkeypatch) -> None:
    """A complete VCF with no sidecar is what the crash left behind; recovery must
    produce a usable sidecar WITHOUT invoking vg, which is the hours it saves."""
    gz, tmp = staged
    out = tmp / "out" / "deconstruct.vcf"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(VCF_TEXT)

    def must_not_run(*a, **k):
        raise AssertionError("recovery re-ran vg; the point is that it does not")
    monkeypatch.setattr(V, "_run", must_not_run)

    V.attach_deconstruct_prov(gz, out)
    prov = json.loads((tmp / "out" / ("deconstruct.vcf" + store.PROV_SUFFIX)).read_text())
    import hashlib
    assert prov["inputs"]["gfa"] == hashlib.sha256(GFA_TEXT).hexdigest()
    assert prov["params"]["prov_recovered_after_crash"] is True, (
        "a reconstructed sidecar must SAY it was reconstructed, or it is "
        "indistinguishable from provenance recorded by the run itself")


def test_recovery_refuses_to_overwrite_real_provenance(staged) -> None:
    """⚠ Replacing a genuine sidecar with an after-the-fact reconstruction would
    look identical and be weaker evidence. It must refuse."""
    gz, tmp = staged
    out = tmp / "out" / "deconstruct.vcf"
    V.deconstruct(gz, out)
    with pytest.raises(FileExistsError, match="real provenance"):
        V.attach_deconstruct_prov(gz, out)


def test_recovery_refuses_a_missing_or_empty_vcf(staged) -> None:
    gz, tmp = staged
    out = tmp / "out" / "deconstruct.vcf"
    with pytest.raises(FileNotFoundError, match="nothing to attach"):
        V.attach_deconstruct_prov(gz, out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("")
    with pytest.raises(FileNotFoundError, match="nothing to attach"):
        V.attach_deconstruct_prov(gz, out)
