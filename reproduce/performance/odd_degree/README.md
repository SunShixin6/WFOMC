# Odd-Degree Benchmark Submodule

This package provides a standalone odd-degree performance runner under
`reproduce/performance/odd_degree/`.

## Included Files

- `benchmarks.py`: fixed runtime configuration and defaults
- `main.py`: entry point and CLI argument parsing
- `runner.py`: per-point execution, timeout/memory handling, CSV writing
- `render.py`: runtime plotting helper
- `transforms.py`: model transformation for odd-degree constraints
- `ganak_odd_degree.py`: CNF conversion + external counters (`ganak`/`approxmc`)

## Actual Default Setup

The authoritative defaults come from `benchmarks.py`:

- domain sizes: `n=5..9`
- directed edge count: `k=2n`
- odd-degree groups: `m=2,4,6`
- algorithms: `approxmc`, `ganak`, `incremental3`
- timeout: `100` seconds
- flush every: `1`

Note: the console banner text in `main.py` may still mention older values.
Use `benchmarks.py` values above as the actual runtime configuration.

## Run

From repository root:

```bash
PYTHONPATH=src:. python -m reproduce.performance.odd_degree.main
```

Or via helper script:

```bash
bash reproduce/performance/odd_degree/run_odd_degree.sh
```

Optional CLI arguments:

```bash
PYTHONPATH=src:. python -m reproduce.performance.odd_degree.main \
  --models-path /absolute/path/to/models \
  --timeout-seconds 100 \
  --flush-every 1
```

## External Counters

When using `ganak` or `approxmc`, binaries are resolved from:

1. environment variables: `GANAK_BIN`, `APPROXMC_BIN`
2. system `PATH`

Example:

```bash
export GANAK_BIN=/path/to/ganak
export APPROXMC_BIN=/path/to/approxmc
```

## Outputs

- CSV root: `reproduce/performance/performance_results/raw_data/odd_degree/results_<timestamp>/`
- CSV files: `odd-degree-n/m-odd-degree-graph-sc2/m-odd-degree-graph-sc2_m{2|4|6}.csv`
- CNF cache: `reproduce/performance/performance_results/raw_data/cnf/odd_degree/m-odd-degree-graph-sc2/`
- Runtime figures (full): `reproduce/performance/performance_results/Section6/`
- Runtime figures (publish copy): `reproduce/results/Section6/`
