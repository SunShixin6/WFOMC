# WFOMC with Counting and Modulo Counting Quantifiers

This repository contains the code accompanying the paper:

> **A Fast Model Counting Algorithm for Two-Variable Logic with Counting and Modulo Counting Quantifiers**
>  Shixin Sun, Astrid Klipfel, Ondřej Kuželka, Yuanhong Wang, Yi Chang

The `modk` branch contains the implementation of `INCREMENTALWFOMC3`, a lifted algorithm for weighted first-order model counting (WFOMC) on the two-variable fragment **C²** and its modulo counting extension **C²_mod**. 

------
## Reproducibility

For the paper reproduction scripts, setup details, and expected outputs, please
refer to [reproduce/README.md](./reproduce/README.md).

Minimal entry point from the repository root:

```bash
bash reproduce/run_all.sh
```
------

## Quick Start

```bash
uv sync
uv run wfomc -i models/modk/0mod2-regular-graph.wfomcs -a incremental3
```

A successful run produces output ending with lines such as:

```text
WFOMC time: 0.12028004229068756
WFOMC (arbitrary precision): 64
WFOMC (round): 64 (exp(4.158883083359671856503392729))
```

> A successful run typically prints a parsing-time message, the exact WFOMC, its rounded exponential form, and several additional log lines.

## Requirements

- Python 3.11 or newer
- A working C/C++ build toolchain (required by the native dependency `pynauty`)

## Installation

We recommend [uv](https://github.com/astral-sh/uv) for dependency management. Install it via pip or follow the [official instructions](https://github.com/astral-sh/uv).

```bash
pip install uv
```

Sync the dependencies:

```bash
uv sync
```

Alternatively, install in editable mode with pip:

```bash
python -m pip install -e .
```

------

## Usage

```bash
uv run wfomc -i <input> [-o <output_dir>] [-a <algo>] [--debug]
```

| Option            | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| `-i <input>`      | Path to a `.wfomcs` input file                               |
| `-a <algo>`       | (Optional) Algorithm to use; defaults to fastv2 (see table below) |
| `-o <output_dir>` | (Optional) Output directory for logs; defaults to `./check-points`, and the solver writes logs to `OUTPUT_DIR/log.txt` |
| `--debug`         | (Optional) Enable debug logging                              |

**For the workflow used in the paper, use:**

```bash
uv run wfomc -i <input> -a incremental3
```

**Available algorithms:**

| Algorithm      | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `incremental3` | **Recommended.** The main algorithm introduced in this paper, supporting counting and modulo counting quantifiers directly. |
| `standard`     | Standard WFOMC (Beame et al., 2015)                          |
| `fast`         | Fast WFOMC (van Bremen & Kuželka, 2021)                      |
| `fastv2`       | Optimized fast WFOMC                                         |
| `incremental`  | IncrementalWFOMC with linear order axiom (Tóth & Kuželka, 2023) |
| `recursive`    | RecursiveWFOMC with linear order axiom (Meng et al., 2024)   |

> **Note:** 
>
> - Only `incremental3` supports modulo counting quantifiers (`\exists_{rmodk}`). For all modulo counting benchmarks, use `-a incremental3`.
>
> - If the input contains a linear order predicate (`LEQ`), use `incremental`, `recursive`, or `incremental3`.



------

## Input Format

See [Input-format.md](./Input-format.md) for the full input specification.

## Example Input Files

**2-regular graphs** — each vertex has degree exactly 2, using counting quantifiers:

```text
\forall X: (~E(X,X)) &
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (\exists_{=2} Y: (E(X,Y)))

V = 7
```

**2-colored graphs** — each vertex is red or black; adjacent vertices have different colors:

```text
\forall X: (\forall Y: ((E(X,Y) -> E(Y,X)) &
                        (R(X) | B(X)) &
                        (~R(X) | ~B(X)) &
                        (E(X,Y) -> ~(R(X) & R(Y)) & ~(B(X) & B(Y)))))

V = 7
```

**0mod2-regular graphs** — each vertex has even degree (modulo counting):

```text
\forall X: (~E(X,X)) &
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (\exists_{0mod2} Y: (E(X,Y)))

V = 5
```

**m-odd-degree graphs**

The file below (`m-odd-degree-graph-sc2.wfomcs`) counts undirected graphs on 4 vertices with exactly 0 odd-degree vertices and exactly 3 edges. 

```text
\forall X: (~E(X,X)) &
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (P(X) <-> (~Odd(X) & A(X) & C(X))) &
\forall X: (\forall Y: (P(X) & B(X,Y) -> U(Y))) &
\forall X: (\forall Y: (~P(X) -> (B(X,Y) <-> E(X,Y)))) &
\forall X: (Odd(X) | A(X)) &
\forall X: (A(X) | C(X)) &

\forall X: (\exists_{1mod2} Y: (B(X, Y))) &

\exists_{=0} X: (Odd(X)) &


\exists_{=1} X: (U(X))

n = 4
1 -1 C
|E| = 6 
```


------



## References

Please refer to [reference.bib](reference.bib) for the references of the algorithms.

## License

This project is released under the MIT License. See the [LICENSE](LICENSE) file for the full license text.