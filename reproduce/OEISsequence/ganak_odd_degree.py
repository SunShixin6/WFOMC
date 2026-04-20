"""Dedicated undirected model-counting script for m-odd-degree-graph.

Purpose and assumptions:
- Targets m-odd-degree-graph-origin.wfomcs with hardcoded undirected semantics:
    symmetry constraints + counting |E| only on c1 < c2.
- Builds CNF and invokes Ganak/ApproxMC for model counting; k is interpreted as
    the number of undirected edges (max n(n-1)/2).
- The input wfomcs file is only used as a naming placeholder; formula parsing is
    intentionally bypassed.

Quick usage:
- Run this script directly to generate odd_degree.csv in the script directory,
  and a second copy in reproduce/results/Appendix.D/Table_2.csv.
- Or call Fo2Counter(file_path, n, m, k, counter="ganak"/"approxmc").
- Ganak/ApproxMC binaries can be configured via GANAK_BIN/APPROXMC_BIN
    environment variables or CLI arguments.
"""

import csv
import os
import sys
import argparse
import copy
import sympy
import re
import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, Set, Iterable, Tuple
from itertools import product, combinations
from logzero import logger, loglevel, logfile

# Allow direct execution via:
# python reproduce/OEISsequence/ganak_odd_degree.py
if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[2]
    for candidate in (repo_root, repo_root / "src"):
        candidate_str = str(candidate)
        if candidate_str not in sys.path:
            sys.path.insert(0, candidate_str)

from wfomc import parse_input
from wfomc.fol.sc2 import SC2
from wfomc.fol.syntax import Const, X, Y, QFFormula, AtomicFormula, Pred, top

from pysat.formula import CNF
from pysat.solvers import Solver
from pysat.formula import CNF
from pysat.solvers import Solver
from pysat.card import CardEnc

from reproduce.utils.smoke_profile import resolve_oeis_max_n
from reproduce.utils.logging_control import env_disables_python_logging


GANAK_ENV_VAR = "GANAK_BIN"
APPROXMC_ENV_VAR = "APPROXMC_BIN"
APPROXMC_SEED = 42

# Resolve from environment first, then fall back to command names in PATH.
ganak_path = os.environ.get(GANAK_ENV_VAR, "ganak")
approxmc_path = os.environ.get(APPROXMC_ENV_VAR, "approxmc")


class CNFContext:
    """CNF conversion context for odd-degree undirected graph grounding."""

    def __init__(self, file_path, n, m, k):
        self.file_path = file_path  # Input model path (used for output naming).
        path = Path(self.file_path)
        self.file_name = path.name  # Input model filename.
        # Write CNF/clause artifacts into the script directory.
        self.file_dir = Path(__file__).resolve().parent
        #
        self.domain = {Const(str(i)) for i in range(
            n)}  # Domain generated from the provided n.
        self.k = k  # Number of undirected edges.
        self.m = m
        #
        self.expr = sympy.true  # Final symbolic expression placeholder.
        self.atom_to_id: Dict[AtomicFormula, int] = {}  # AtomicFormula -> CNF variable id.
        self.sym_to_id: Dict[sympy.Symbol, int] = {}  # Sympy symbol -> CNF variable id.
        self.next_var_id = 1  # Next available CNF variable id.
        #
        self.cnf_dir = os.path.join(self.file_dir, "cnf")
        cnf_file_name = f"{os.path.splitext(self.file_name)[0]}_n_{len(self.domain)}_m_{self.m}_k_{self.k}.cnf" # Output CNF filename.
        clause_file_name = f"{os.path.splitext(self.file_name)[0]}_n_{len(self.domain)}_m_{self.m}_k_{self.k}.txt" # Output clause filename.
        self.cnf_path = os.path.join(self.cnf_dir, cnf_file_name)
        self.clause_path = os.path.join(self.cnf_dir, clause_file_name)
        self.clauses: list[list[int]] = []  # Collected CNF clauses.

    def convert(self):
        """Run specialized m-odd-degree grounding and populate clauses."""
        os.makedirs(os.path.dirname(self.cnf_path),
                    exist_ok=True)  # Ensure output directory exists.

        self._ground_m_odd_degree_formulas()

    def _ground_m_odd_degree_formulas(self):
        """
        Specialized grounding routine for m-odd-degree-graph-origin.wfomcs.
        Encodes the following formulas:
        1. \forall X: (~E(X,X))
        2. \forall X: (\forall Y: (E(X,Y) -> E(Y,X)))
        3. \forall X: (Odd(X) <-> (\exists_{1 mod 2} Y: (E(X, Y))))
        4. \exists_{=m} X: (Odd(X))
        5. |E| = k
        """
        logger.info("Applying specialized grounding for m-odd-degree-graph-origin.")
        domain = self.domain
        e_pred = Pred('E', 2) # Edge predicate.
        odd_pred = Pred('Odd', 1) # Odd-degree predicate.

        ## Pre-register all E(c1, c2) and Odd(c) atoms with unique ids.
        for c1, c2 in product(domain, repeat=2):
            self._register_atom(AtomicFormula(e_pred, (c1, c2), True))
        for c in domain:
            self._register_atom(AtomicFormula(odd_pred, (c,), True))

        ## 1. \forall X: ~E(X,X) (irreflexive)
        for c in domain:
            atom = AtomicFormula(e_pred, (c, c), True)  # E(c, c)
            self.clauses.append([-self.atom_to_id[atom]])  # Add clause: ~E(c, c)

        ## 2. \forall X, Y: E(X,Y) -> E(Y,X) (symmetry)
        # Equivalent to: ~E(X,Y) v E(Y,X)
        for c1, c2 in combinations(domain, 2):
            atom1 = AtomicFormula(e_pred, (c1, c2), True)
            atom2 = AtomicFormula(e_pred, (c2, c1), True)
            var1 = self.atom_to_id[atom1]
            var2 = self.atom_to_id[atom2]
            # E(c1, c2) <-> E(c2, c1)
            self.clauses.append([-var1, var2])
            self.clauses.append([var1, -var2])

        ## 3. \forall X: (Odd(X) <-> (\exists_{1 mod 2} Y: E(X, Y)))
        # For each X, use an XOR chain to encode degree parity.
        for c1 in domain:
            odd_c1_var = self.atom_to_id[AtomicFormula(
                odd_pred, (c1,), True)]

            # Degree only counts edges to distinct neighbors.
            degree_vars = [self.atom_to_id[AtomicFormula(e_pred, (c1, c2), True)]
                           for c2 in domain if c1 != c2]

            # Handle the empty-neighbor case explicitly for completeness.
            if not degree_vars:
                # n=1 has degree 0, so Odd(c1) must be false.
                self.clauses.append([-odd_c1_var])
                continue
            # Build XOR chain when neighbors exist.
            # current_xor_out tracks accumulated parity.
            current_xor_out = degree_vars[0]

            # Build from the second degree variable onward.
            for i in range(1, len(degree_vars)):
                next_var = degree_vars[i]
                # Introduce an auxiliary variable for each XOR step.
                xor_out_new = self.next_var_id
                self.next_var_id += 1

                # Add clauses for xor_out_new <-> (current_xor_out XOR next_var)
                # (~xor_out_new V current_xor_out V next_var)
                self.clauses.append([-xor_out_new, current_xor_out, next_var])
                # (~xor_out_new V ~current_xor_out V ~next_var)
                self.clauses.append([-xor_out_new, -current_xor_out, -next_var])
                # (xor_out_new V ~current_xor_out V next_var)
                self.clauses.append([xor_out_new, -current_xor_out, next_var])
                # (xor_out_new V current_xor_out V ~next_var)
                self.clauses.append([xor_out_new, current_xor_out, -next_var])
                
                # Update chain state.
                current_xor_out = xor_out_new

            # Equate final XOR result with Odd(c1).
            # Odd(c1) <-> current_xor_out
            self.clauses.append([-odd_c1_var, current_xor_out])
            self.clauses.append([odd_c1_var, -current_xor_out])


        ## 4. \exists_{=m} X: (Odd(X))
        # Collect variables corresponding to Odd(c).
        odd_vars = [self.atom_to_id[AtomicFormula(
            odd_pred, (c,), True)] for c in domain]
        # Use CardEnc.equals to encode an exact-m cardinality constraint.
        card_clauses = CardEnc.equals(
            lits=odd_vars, bound=self.m, top_id=self.next_var_id - 1)
        self.clauses.extend(card_clauses.clauses)
        self.next_var_id = max(self.next_var_id, card_clauses.nv + 1)

        ## 5. |E| = k (cardinality constraint)
        # Because the graph is undirected (enforced by symmetry), only count c1 < c2.
        edge_vars = []
        # Sort domain constants by name for deterministic ordering.
        sorted_domain = sorted(list(domain), key=lambda c: c.name)
        for c1, c2 in combinations(sorted_domain, 2):
            # Use one directed atom to represent each undirected edge.
            atom = AtomicFormula(e_pred, (c1, c2), True)
            edge_vars.append(self.atom_to_id[atom])


        # Add the cardinality constraint |E| = k.
        logger.info(
            f"Adding cardinality constraint |E| = {self.k} on {len(edge_vars)} edge variables.")
        if self.k > len(edge_vars):
            self.clauses.append([])  # Impossible constraint, force UNSAT.
        elif edge_vars or self.k == 0:
            card_clauses_e = CardEnc.equals(
                lits=edge_vars, bound=self.k, top_id=self.next_var_id - 1)
            self.clauses.extend(card_clauses_e.clauses)
            self.next_var_id = max(self.next_var_id, card_clauses_e.nv + 1)

    def _register_atom(self, atom: AtomicFormula):
        """Register an atom and allocate a new id if needed."""
        if atom not in self.atom_to_id:
            self.atom_to_id[atom] = self.next_var_id  # Register atomic formula.
            self.sym_to_id[atom.expr] = self.next_var_id  # Register mapped sympy symbol.
            self.next_var_id += 1  # Advance id counter.

    def dump(self):
        """Write DIMACS CNF output to cnf_path."""
        num_vars = self.next_var_id - 1
        num_clauses = len(self.clauses)
        with open(self.cnf_path, 'w', encoding='utf-8') as f:
            f.write(f"p cnf {num_vars} {num_clauses}\n")
            for clause in self.clauses:
                f.write(" ".join(map(str, clause)) + " 0\n")

        logger.info(
            f"CNF file with {num_vars} vars and {num_clauses} clauses written to: {self.cnf_path}")


    @staticmethod
    def model_count_ganak(cnf_path: str) -> int:
        """Run Ganak and parse model count from its output."""
        result = subprocess.run([ganak_path, cnf_path],
                                capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Ganak execution failed: {result.stderr.strip()}")

        match = re.search(r"(?:s mc|c s exact arb int)\s+(\d+)", result.stdout)
        if match:
            return int(match.group(1))

        raise RuntimeError(f"Could not parse Ganak output: {result.stdout.strip()}")

    @staticmethod
    def model_count_approxmc(cnf_path: str, epsilon: float, delta: float) -> int:
        """Run ApproxMC and parse the integer estimate from 's mc' output."""
        cmd = [
            approxmc_path,
            f"--epsilon={epsilon}",
            f"--delta={delta}",
            f"--seed={APPROXMC_SEED}",
            cnf_path,
        ]
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"ApproxMC execution failed: {result.stderr.strip()}"
            )

        for line in result.stdout.splitlines():
            if line.startswith("s mc"):
                try:
                    return int(line.split()[-1])
                except (ValueError, IndexError):
                    break

        raise RuntimeError(f"Could not parse ApproxMC output: {result.stdout.strip()}")


def Fo2Counter(file_path, n, m, k, counter="ganak", epsilon=0.05, delta=0.1):
    """Single entry point: build odd-degree CNF and run the selected counter.

    Args:
    - file_path: Input path used for output naming.
    - n: Domain size.
    - m: Number of odd-degree nodes (exactly m).
    - k: Number of undirected edges (upper bound n(n-1)/2).
    - counter: "ganak" or "approxmc".
    - epsilon/delta: ApproxMC accuracy parameters.

    Returns: Integer model count.
    Raises: RuntimeError/FileNotFoundError/ValueError on failure.
    """
    if not os.path.exists(file_path):  # Validate input placeholder file exists.
        raise FileNotFoundError(f"Input file does not exist: {file_path}")


    context = CNFContext(file_path, n, m, k)  # Build CNF context.

    context.convert()  # Perform grounding/encoding.
    context.dump()  # Write CNF to disk.

    if counter == "ganak":
        count = CNFContext.model_count_ganak(context.cnf_path)
    elif counter == "approxmc":
        count = CNFContext.model_count_approxmc(context.cnf_path, epsilon, delta)
    else:
        raise ValueError(f"Unknown counter: {counter}")

    return count

    # logger.info(
    #     f"Result:\n InputFile: {file_path}\n Domain Size: {n}\n K:{k}\n M:{m}\n Model Count: {count}\n")


if __name__ == '__main__':
    quiet_python_logging = env_disables_python_logging()
    logger.setLevel(logging.ERROR if quiet_python_logging else logging.INFO)

    parser = argparse.ArgumentParser(
        description="m-odd-degree-graph counting script with CLI path overrides"
    )
    parser.add_argument(
        "--ganak-path",
        default=ganak_path,
        help=f"Path to Ganak executable (default: ${GANAK_ENV_VAR} or 'ganak')"
    )
    parser.add_argument(
        "--approxmc-path",
        default=approxmc_path,
        help=f"Path to ApproxMC executable (default: ${APPROXMC_ENV_VAR} or 'approxmc')"
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Enable smoke profile (sets max_n to 4 unless --max-n is provided).",
    )
    parser.add_argument(
        "--max-n",
        type=int,
        default=None,
        help="Maximum domain size n (default: 10; smoke default: 4)",
    )
    args = parser.parse_args()

    ganak_path = args.ganak_path
    approxmc_path = args.approxmc_path
    logger.info(f"Using Ganak path: {ganak_path}")
    logger.info(f"Using ApproxMC path: {approxmc_path}")

    script_dir = Path(__file__).resolve().parent
    file_path = str(script_dir / "m-odd-degree-graph.wfomcs") # Input path in script directory; used only for naming.
    max_n = resolve_oeis_max_n(args.max_n, args.smoke, default_max_n=10, smoke_max_n=4)
    if max_n < 1:
        raise ValueError("--max-n must be >= 1")
    logger.info(f"Using max_n={max_n}")
    script_output_path = script_dir / "odd_degree.csv"
    appendix_output_path = script_dir.parent / "results" / "Appendix.D" / "Table_2.csv"
    os.makedirs(appendix_output_path.parent, exist_ok=True)
    logger.info(f"Writing script-level CSV to: {script_output_path}")
    logger.info(f"Writing Appendix D CSV to: {appendix_output_path}")

    # Write both CSV outputs row by row with identical content.
    with open(script_output_path, "w", newline='') as script_file, open(appendix_output_path, "w", newline='') as appendix_file:
        script_writer = csv.writer(script_file)
        appendix_writer = csv.writer(appendix_file)
        # Sweep all n up to max_n.
        for n in range(1, max_n + 1):
            max_k = n * (n - 1) // 2
            # Write CSV header.
            if n == 1:
                header = ["n", "m"] + [f"k={k}" for k in range(max_k + 1)]
                script_writer.writerow(header)
                appendix_writer.writerow(header)
            for m in range(0, n + 1, 2):  # m must be even due to parity constraints.
                row = [n, m]
                for k in range(max_k + 1):
                    try:
                        count = Fo2Counter(file_path, n=n, m=m, k=k)
                    except Exception as exc:
                        raise RuntimeError(
                            f"Counting failed at n={n}, m={m}, k={k}"
                        ) from exc
                    row.append(count)
                    if not quiet_python_logging:
                        print(f"n={n}, m={m}, k={k} -> valid model count: {count}")
                script_writer.writerow(row)
                appendix_writer.writerow(row)
                script_file.flush()
                appendix_file.flush()

