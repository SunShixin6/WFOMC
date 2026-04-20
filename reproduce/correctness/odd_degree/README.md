# Figure13 Model-Count Reproduction

This package reproduces odd-degree model-count experiments under
`reproduce/correctness/odd_degree/`.

## Included Files

- `benchmarks.py`: Figure13 sweep definitions and defaults
- `main.py`: entry point and CLI parsing
- `runner.py`: isolated execution, timeout handling, CSV writing
- `render.py`: model-count plotting helper
- `run_odd_degree.sh`: convenience launcher

## Actual Sweep Definitions

The authoritative sweep settings are in `benchmarks.py`:

1. `vary_n`
- `n=2..8`
- `m=2`
- `k=1` (undirected)

2. `vary_k`
- `n=5`
- `m=2`
- `k=1..9` (undirected)

3. `vary_m`
- `n=8`
- `k=3` (undirected)
- `m=0,2,4,6`

Notes:

- `k` in the paper-facing interface is undirected edge count.
- Internal WFOMC path uses directed edge count and converts via `|E| = 2k`.
- Algorithms are fixed as `incremental3`, `ganak`, `approxmc`.

## Run

From repository root:

```bash
PYTHONPATH=src:. python -m reproduce.correctness.odd_degree.main
```

Or via shell script:

```bash
bash reproduce/correctness/odd_degree/run_odd_degree.sh
```

Optional CLI arguments (default values shown):

```bash
PYTHONPATH=src:. python -m reproduce.correctness.odd_degree.main \
  --models-path /absolute/path/to/models \
  --timeout-seconds 10000 \
  --epsilon 0.05 \
  --delta 0.1 \
  --flush-every 1
```

## Outputs

- CSV root:
  `reproduce/correctness/correctness_results/raw_data/odd_degree/results_<timestamp>/`
- CSV files by sweep:
  - `vary_n/m-odd-degree-graph-sc2/m-odd-degree-graph-sc2_vary_n_m2_k1.csv`
  - `vary_k/m-odd-degree-graph-sc2/m-odd-degree-graph-sc2_vary_k_n5_m2.csv`
  - `vary_m/m-odd-degree-graph-sc2/m-odd-degree-graph-sc2_vary_m_n8_k3.csv`
- CNF cache:
  `reproduce/correctness/correctness_results/raw_data/cnf/odd_degree/m-odd-degree-graph-sc2/`
- PDF plots:
  - `reproduce/correctness/correctness_results/AppendixB.1/`
  - `reproduce/results/AppendixB.1/`

CSV columns:

- `timestamp`
- `sweep`
- `variable_name`
- `variable_value`
- `formula`
- `domain_size`
- `k_value`
- `m_value`
- `algorithm`
- `model_count`
- `status`
