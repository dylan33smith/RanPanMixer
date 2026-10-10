"""Run independent jobs across cores, with the worker count bounded by MEMORY.

⚠ WHY THIS EXISTS. Nearly everything expensive in this project is embarrassingly
parallel — a release per (target, `tau`, seed), a chromosome per deconstruct, a draw
per trial — and it was all being run serially on a 48-core, 251 GB machine. Costs
were being quoted in hours that were really minutes, and a wrong cost estimate
changes which experiments look affordable, so it changes the design. Dylan's
instruction, 2026-10-10: parallelise future runs wherever they can be.

⚠ AND THE BOUND IS MEMORY, NOT CORES, which is the whole reason this is a helper
rather than a one-line `Pool(cpu_count())`. A chr21 forward pass stores an alpha
array of `T * K * 8` bytes — about 215 MB — and a release needs the potentials and
the backward message too, so three of those. `state_marginals` needs about five.
Forty-eight of those at once is real memory, and on the largest autosome each array
is five times bigger. A pool sized by cores will swap or be killed; a pool sized by
memory will not.

Use `plan_workers` to size the pool from a per-job estimate, and `run_jobs` to
execute and ASSERT COMPLETENESS. A parallel map that silently drops a failed job
cannot tell "found nothing" from "ran nothing" — the rule `CLAUDE.md` states for
fan-out, which applies no less when the fan-out is in-process.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass


def available_ram_bytes() -> int:
    """RAM we may actually use, from MemAvailable rather than MemTotal.

    MemTotal would count memory other work already holds. MemAvailable is the
    kernel's own estimate of what a new allocation can get without swapping, which
    is the question being asked.
    """
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return 8 << 30          # a deliberately small fallback: better to under-use


def array_bytes(t: int, k: int, n_arrays: int = 3, dtype_bytes: int = 8) -> int:
    """Peak bytes for a job holding `n_arrays` dense (T, K) float64 arrays.

    `n_arrays = 3` covers a release (potentials, alphas, backward); use 5 for
    `state_marginals`, which also materialises the product and the normalised copy.
    """
    return int(t) * int(k) * int(dtype_bytes) * int(n_arrays)


@dataclass(frozen=True)
class Plan:
    workers: int
    per_job_bytes: int
    headroom_bytes: int
    reason: str

    def __str__(self) -> str:
        return (f"{self.workers} workers "
                f"({self.per_job_bytes / 1e9:.2f} GB/job, "
                f"{self.headroom_bytes / 1e9:.0f} GB usable) — {self.reason}")


def plan_workers(per_job_bytes: int, *, n_jobs: int | None = None,
                 reserve_cores: int = 2, ram_fraction: float = 0.7) -> Plan:
    """How many workers to run, and why — bounded by cores, memory, and job count.

    `ram_fraction` leaves a margin, because a per-job estimate is an estimate:
    numpy temporaries, the interpreter and the page cache all cost more than the
    arrays being counted. 0.7 has no theory behind it; it is a margin, and the
    `reason` string says which constraint actually bound so an over-tight estimate
    is visible rather than mysterious.
    """
    if per_job_bytes <= 0:
        raise ValueError(f"per_job_bytes must be positive, got {per_job_bytes}")
    cores = max(1, (os.cpu_count() or 2) - reserve_cores)
    usable = int(available_ram_bytes() * ram_fraction)
    by_ram = max(1, usable // per_job_bytes)
    cap = cores if n_jobs is None else min(cores, max(1, n_jobs))
    workers = int(min(cap, by_ram))
    if workers == by_ram and by_ram < cap:
        why = f"MEMORY-bound: {usable / 1e9:.0f} GB usable / {per_job_bytes / 1e9:.2f} GB"
    elif n_jobs is not None and workers == n_jobs:
        why = f"job-count-bound: only {n_jobs} job(s)"
    else:
        why = f"CORE-bound: {os.cpu_count()} cores less {reserve_cores} reserved"
    return Plan(workers=workers, per_job_bytes=per_job_bytes,
                headroom_bytes=usable, reason=why)


def run_jobs(fn, jobs: list, *, workers: int, label: str = "jobs",
             progress_every: int = 0):
    """Map `fn` over `jobs` in parallel, in order, asserting every one finished.

    ⚠ A failure RAISES rather than returning a short list. The alternative is a
    result set that looks complete and is not, which is the failure mode this
    project has already been bitten by in a different guise: a filter cannot
    distinguish "nothing matched" from "nothing ran".

    `fn` must be a module-level callable and `jobs` must be picklable, because the
    workers are processes — threads would serialise on the GIL for the numpy-heavy
    work this exists for.
    """
    if not jobs:
        raise ValueError(f"no {label} to run")
    if workers < 1:
        raise ValueError(f"workers must be >= 1, got {workers}")
    if workers == 1:
        return [fn(j) for j in jobs]

    out: list = [None] * len(jobs)
    done = 0
    failures: list[str] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(fn, j): i for i, j in enumerate(jobs)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                out[i] = fut.result()
            except Exception as e:                      # noqa: BLE001
                failures.append(f"  {label}[{i}] = {jobs[i]!r}: "
                                f"{type(e).__name__}: {e}")
            done += 1
            if progress_every and done % progress_every == 0:
                print(f"  {done}/{len(jobs)} {label}", flush=True)

    if failures:
        raise RuntimeError(
            f"{len(failures)} of {len(jobs)} {label} FAILED, so this result set is "
            f"incomplete and must not be quoted:\n" + "\n".join(failures[:10]) +
            (f"\n  ... and {len(failures) - 10} more" if len(failures) > 10 else ""))
    if any(o is None for o in out):
        n = sum(o is None for o in out)
        raise RuntimeError(
            f"{n} of {len(jobs)} {label} returned nothing and did not raise. That "
            f"should be impossible; do not quote this result set.")
    return out
