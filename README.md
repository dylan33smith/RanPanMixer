# RanPanMixer

A randomized path-release mechanism for pangenome graphs. Given a public graph built
from a public cohort, it releases an external target's path as a single draw from a
tilted cohort HMM. Any two input paths produce output distributions within total
variation distance `tau` of each other.

## Layout

| path | what |
|---|---|
| `scripts/` | pipeline code: preprocessing first, then the sampler and evaluation |
| `configs/` | run parameters and the input manifest (URLs, checksums) |
| `logs/` | run logs and exit-status sentinels (contents gitignored) |
| `data/` | symlink to `/data/ds85/RanPanMixer/data` (gitignored, see below) |

## Data

The inputs are too large for git: the pangenome VCF alone is 5.7 GB, and the full
input set is about 43 GB. `data/` is a symlink to storage on `/data`. Anything
needed to rebuild it, such as source URLs and sha256 checksums, lives in `configs/`.

## Hard rule

The baseline model and the output support are built from the graph and the cohort
only, never from the target. Preprocessing therefore never takes a target as an
argument. Anything target-specific, including leave-one-out views of the cohort,
happens later, in a separate step.
