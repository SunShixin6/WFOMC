# Reproducibility Guide

This directory contains scripts for reproducing the experimental results
reported in the paper, including benchmark families for `C²` and `C²_mod`.

For solver usage and algorithmic background, see the top-level
[README](../README.md).

## Prerequisites

- Python `>=3.11`
- Dependencies installed from repository root (recommended: `uv sync`)
- Project virtual environment at `.venv`
- The run scripts use `.venv` when available and automatically fall back to `uv run python` otherwise.

## Quick Start

Run from repository root:

```bash
# 1) Validate the pipeline without launching heavy jobs
bash reproduce/run_all.sh --dry-run

# 2) Lightweight end-to-end sanity check
bash reproduce/run_all.sh --smoke --disable-python-logging

# 3) Full reproduction
bash reproduce/run_all.sh
```

## Pipeline Overview

`reproduce/run_all.sh` executes the following stages in order:

1. `reproduce.correctness.main`
2. `reproduce.correctness.odd_degree.main`
3. `reproduce.performance.main`
4. `reproduce.performance.odd_degree.main`
5. `reproduce/OEISsequence/run_ganak_odd_degree.sh`

## External Model Counters

Some stages invoke external model counters:

- `GANAK_BIN`: Ganak executable (exact counting)
- `APPROXMC_BIN`: ApproxMC executable (approximate counting)

For a full end-to-end run (`bash reproduce/run_all.sh`), both counters must be
available (via environment variables or `PATH`).

Lookup order:

1. Environment variable value, if set.
2. Command from `PATH` (`ganak`, `approxmc`).

Example (Linux/macOS):

```bash
export GANAK_BIN=/absolute/path/to/ganak
export APPROXMC_BIN=/absolute/path/to/approxmc
```

Example (PowerShell):

```powershell
$env:GANAK_BIN = "C:\\path\\to\\ganak.exe"
$env:APPROXMC_BIN = "C:\\path\\to\\approxmc.exe"
```

## Model File Resolution

Model files are resolved in this order:

1. `models/<model_filename>`
2. Any nested subdirectory under `models/` (recursive filename match)
3. `reproduce/models/<model_filename>`
4. Any nested subdirectory under `reproduce/models/` (recursive filename match)

## Logging and Quiet Mode

Each run creates timestamped logs under:

- `reproduce/logs/run_<timestamp>/` (for `run_all.sh`)
- `reproduce/logs/single_<timestamp>/` (for single-stage scripts)

Useful flags:

- `--disable-python-logging`
- `--enable-python-logging`

Equivalent environment variable:

- `REPRO_DISABLE_PYTHON_LOGGING=1` disables Python logging.
- `REPRO_DISABLE_PYTHON_LOGGING=0` keeps the default behavior.

To keep terminal output minimal while preserving full logs:

```bash
bash reproduce/run_all.sh --smoke --disable-python-logging > /tmp/repro_smoke.out 2>&1
```

## Smoke Mode

Smoke mode runs the same pipeline with reduced benchmark ranges for quick
validation:

```bash
bash reproduce/run_all.sh --smoke
```

Smoke mode can also be enabled with `REPRO_SMOKE=1`.

## Expected Outputs

After a successful full run, the following outputs should exist:

- `reproduce/correctness/correctness_results/raw_data/results_<timestamp>/`
- `reproduce/correctness/correctness_results/raw_data/odd_degree/results_<timestamp>/`
- `reproduce/performance/performance_results/raw_data/results_<timestamp>/`
- `reproduce/performance/performance_results/raw_data/odd_degree/results_<timestamp>/`
- `reproduce/results/Appendix.D/Table_2.csv`

## Mapping from Paper Figures and Tables to Scripts

- `Figure_5_to_9` -> `reproduce.performance.main` (Section 6 runtime plots)
- `Figure_10` -> `reproduce.performance.odd_degree.main`
- `Figure_11_to_12` -> `reproduce.correctness.main`
- `Figure_13` -> `reproduce.correctness.odd_degree.main`
- `Figure_14_to_16` -> `reproduce.performance.main` (Appendix B.2 memory plots)
- `Table_2` -> `reproduce/OEISsequence/ganak_odd_degree.py`

## Acceptance Checks

Minimal reviewer-facing checks:

1. `bash reproduce/run_all.sh` exits with code `0`.
2. One-command basic check passes:

```bash
bash reproduce/utils/check_results.sh
```

3. (Optional) Strict check passes for expected figure filenames and CSV coverage:

```bash
bash reproduce/utils/check_results.sh --strict
```

4. (Optional) Export a machine-readable report:

```bash
bash reproduce/utils/check_results.sh --strict --json reproduce/results/check_report.json
```

Quick checks:

```bash
test -s reproduce/results/Appendix.D/Table_2.csv
find reproduce/results/Section6 -name '*.pdf' | head -n 3
find reproduce/results/AppendixB.1 -name '*.pdf' | head -n 3
find reproduce/results/AppendixB.2 -name '*.pdf' | head -n 3
```

## Individual Entry Points

Run from repository root.

### Correctness

```bash
bash reproduce/correctness/run_correctness.sh
```

Defaults: timeout `10000`, epsilon `0.05`, delta `0.1`.

### Correctness (Odd-Degree, Figure 13)

```bash
bash reproduce/correctness/odd_degree/run_odd_degree.sh
```

Defaults: timeout `10000`, epsilon `0.05`, delta `0.1`.

### Performance

```bash
bash reproduce/performance/run_performance.sh
```

Default timeout: `10000`.

### Performance (Odd-Degree, Figure 10)

```bash
bash reproduce/performance/odd_degree/run_odd_degree.sh
```

Defaults: timeout `100`, epsilon `0.05`, delta `0.1`, with fixed `m ∈ {2,4,6}`.

### OEIS / Table 2 Generation

```bash
bash reproduce/OEISsequence/run_ganak_odd_degree.sh
```

Equivalent direct call:

```bash
PYTHONPATH=src:. python reproduce/OEISsequence/ganak_odd_degree.py
```

## Cleanup Generated Outputs

Use the cleanup helper to remove generated artifacts:

```bash
# List available figure/table keys
bash reproduce/utils/clean_results.sh --list

# Remove all generated outputs and related logs
bash reproduce/utils/clean_results.sh --all --yes

# Remove selected artifacts only
bash reproduce/utils/clean_results.sh --figure Figure_11_to_12 --figure Figure_10 --yes
```
