"""Writing artifacts, with provenance attached whether you want it or not.

This module exists for one reason. `CLAUDE.md` says every number must be
traceable to an artifact, and `docs/data.md` says a file not registered there is
not quotable. That is a rule someone has to remember. Routing every write
through `save()` turns it into something the code does on your behalf: you
cannot produce an artifact without also producing the record of where it came
from, because they are the same function call.

It is deliberately small. If this file grows past a couple of hundred lines it
has stopped paying for itself -- the point is a chokepoint, not a framework.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

PROV_SUFFIX = ".prov.json"


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    """Content digest of a file on disk."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def git_rev(repo: str | Path | None = None) -> str:
    """The commit the code was at. 'unknown' rather than an exception --
    provenance should never be the thing that fails a run."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo) if repo else None,
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def save(
    path: str | Path,
    obj,
    *,
    inputs: dict[str, str],
    params: dict,
    axis_digest: str | None = None,
    note: str = "",
) -> Path:
    """Write one artifact and its provenance sidecar.

    inputs      name -> sha256 (or a bare identifier for a non-file input).
                This is what makes a stale artifact detectable.
    params      every value that shaped the output. A parameter missing here is
                a parameter nobody can later prove was used -- and several of
                ours (`chain_span` above all) are explicitly provisional, so
                this is the field that matters most.
    axis_digest the site-axis digest this artifact is indexed by, when it has
                one. Row order is the join key for everything and nothing
                upstream checks it; carrying the digest is what lets a reader
                refuse to join two artifacts built against different VCFs.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(obj, np.ndarray):
        np.save(path, obj, allow_pickle=False)
        written = path.with_suffix(".npy") if path.suffix != ".npy" else path
        shape, dtype = list(obj.shape), str(obj.dtype)
    elif isinstance(obj, (dict, list)):
        written = path
        written.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
        shape, dtype = None, "json"
    elif isinstance(obj, str):
        written = path
        written.write_text(obj)
        shape, dtype = None, "text"
    else:
        raise TypeError(f"save() does not handle {type(obj)!r}; add it deliberately")

    prov = {
        "artifact": written.name,
        "sha256": sha256_file(written),
        "shape": shape,
        "dtype": dtype,
        "inputs": inputs,
        "params": params,
        "axis_digest": axis_digest,
        "git_rev": git_rev(Path(__file__).resolve().parents[2]),
        "note": note,
    }
    Path(str(written) + PROV_SUFFIX).write_text(
        json.dumps(prov, indent=2, sort_keys=True) + "\n"
    )
    return written


def load_prov(path: str | Path) -> dict:
    """The provenance record for an artifact, or a clear failure."""
    p = Path(str(path) + PROV_SUFFIX)
    if not p.exists():
        raise FileNotFoundError(
            f"{path} has no provenance sidecar. It was not written by store.save(), "
            f"so nothing records what produced it -- treat it as not quotable."
        )
    return json.loads(p.read_text())


def check_axis(path: str | Path, expected_digest: str) -> None:
    """Refuse to use an artifact built against a different site axis.

    The failure this prevents is the quiet one: a re-sorted or re-filtered VCF
    silently re-labels which variant is which, every downstream number is then
    about the wrong variants, and nothing errors. PanMixer has exactly this
    hazard and its own guard is vacuous -- both sides of its shape assertion
    equal `site_mask.sum()` by construction, so it can never fire.
    """
    got = load_prov(path).get("axis_digest")
    if got != expected_digest:
        raise ValueError(
            f"site-axis mismatch for {path}:\n"
            f"  artifact was built against {got}\n"
            f"  you are joining it against  {expected_digest}\n"
            f"These are different variant orderings. Rebuild, do not override."
        )
