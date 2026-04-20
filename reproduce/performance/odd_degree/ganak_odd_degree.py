from __future__ import annotations

import os
import re
import shutil
import subprocess
from itertools import combinations, product
from pathlib import Path
from typing import Dict

import sympy
from pysat.card import CardEnc

from wfomc.fol.syntax import AtomicFormula, Const, Pred

GANAK_ENV_VAR = "GANAK_BIN"
APPROXMC_ENV_VAR = "APPROXMC_BIN"
APPROXMC_SEED = 42


def _resolve_counter_binary(env_var: str, default_command: str) -> str:
    configured = os.environ.get(env_var, default_command)
    configured_path = Path(configured).expanduser()

    if configured_path.is_file():
        return str(configured_path)

    resolved = shutil.which(configured)
    if resolved:
        return resolved

    raise FileNotFoundError(
        f"Unable to locate counter binary '{configured}'. "
        f"Set {env_var} or install '{default_command}' in PATH."
    )


class CNFContext:
    """Build CNF constraints for the m-odd-degree benchmark model."""

    def __init__(
        self,
        file_path: str | Path,
        n: int,
        m: int,
        k: int,
        output_dir: str | Path | None = None,
    ) -> None:
        self.file_path = Path(file_path)
        self.file_name = self.file_path.name
        self.file_dir = Path(output_dir) if output_dir else self.file_path.parent
        self.file_dir.mkdir(parents=True, exist_ok=True)

        self.domain = {Const(str(i)) for i in range(n)}
        self.k = k
        self.m = m

        self.atom_to_id: Dict[AtomicFormula, int] = {}
        self.sym_to_id: Dict[sympy.Symbol, int] = {}
        self.next_var_id = 1

        stem = self.file_path.stem
        cnf_file_name = f"{stem}_n_{len(self.domain)}_m_{self.m}_k_{self.k}.cnf"
        self.cnf_path = self.file_dir / cnf_file_name
        self.clauses: list[list[int]] = []

    def convert(self) -> None:
        self._ground_m_odd_degree_formulas()

    def _register_atom(self, atom: AtomicFormula) -> None:
        if atom not in self.atom_to_id:
            self.atom_to_id[atom] = self.next_var_id
            self.sym_to_id[atom.expr] = self.next_var_id
            self.next_var_id += 1

    def _ground_m_odd_degree_formulas(self) -> None:
        """
        Ground the following formula family:
        1) forall X: ~E(X, X)
        2) forall X,Y: E(X,Y) <-> E(Y,X)
        3) forall X: Odd(X) <-> (exists_{1 mod 2} Y: E(X, Y))
        4) exists_{=m} X: Odd(X)
        5) |E| = k on undirected edges
        """
        domain = self.domain
        e_pred = Pred("E", 2)
        odd_pred = Pred("Odd", 1)

        for c1, c2 in product(domain, repeat=2):
            self._register_atom(AtomicFormula(e_pred, (c1, c2), True))
        for c in domain:
            self._register_atom(AtomicFormula(odd_pred, (c,), True))

        for c in domain:
            atom = AtomicFormula(e_pred, (c, c), True)
            self.clauses.append([-self.atom_to_id[atom]])

        for c1, c2 in combinations(domain, 2):
            atom1 = AtomicFormula(e_pred, (c1, c2), True)
            atom2 = AtomicFormula(e_pred, (c2, c1), True)
            var1 = self.atom_to_id[atom1]
            var2 = self.atom_to_id[atom2]
            self.clauses.append([-var1, var2])
            self.clauses.append([var1, -var2])

        for c1 in domain:
            odd_var = self.atom_to_id[AtomicFormula(odd_pred, (c1,), True)]
            degree_vars = [
                self.atom_to_id[AtomicFormula(e_pred, (c1, c2), True)]
                for c2 in domain
                if c1 != c2
            ]

            if not degree_vars:
                self.clauses.append([-odd_var])
                continue

            current_xor_out = degree_vars[0]
            for i in range(1, len(degree_vars)):
                next_var = degree_vars[i]
                xor_out_new = self.next_var_id
                self.next_var_id += 1

                self.clauses.append([-xor_out_new, current_xor_out, next_var])
                self.clauses.append([-xor_out_new, -current_xor_out, -next_var])
                self.clauses.append([xor_out_new, -current_xor_out, next_var])
                self.clauses.append([xor_out_new, current_xor_out, -next_var])

                current_xor_out = xor_out_new

            self.clauses.append([-odd_var, current_xor_out])
            self.clauses.append([odd_var, -current_xor_out])

        odd_vars = [
            self.atom_to_id[AtomicFormula(odd_pred, (c,), True)]
            for c in domain
        ]
        if self.m < 0 or self.m > len(odd_vars):
            self.clauses.append([])
        elif odd_vars or self.m == 0:
            odd_card_clauses = CardEnc.equals(
                lits=odd_vars,
                bound=self.m,
                top_id=self.next_var_id - 1,
            )
            self.clauses.extend(odd_card_clauses.clauses)
            self.next_var_id = max(self.next_var_id, odd_card_clauses.nv + 1)

        sorted_domain = sorted(domain, key=lambda c: c.name)
        edge_vars = [
            self.atom_to_id[AtomicFormula(e_pred, (c1, c2), True)]
            for c1, c2 in combinations(sorted_domain, 2)
        ]

        if self.k < 0 or self.k > len(edge_vars):
            self.clauses.append([])
        elif edge_vars or self.k == 0:
            edge_card_clauses = CardEnc.equals(
                lits=edge_vars,
                bound=self.k,
                top_id=self.next_var_id - 1,
            )
            self.clauses.extend(edge_card_clauses.clauses)
            self.next_var_id = max(self.next_var_id, edge_card_clauses.nv + 1)

    def dump(self) -> None:
        num_vars = self.next_var_id - 1
        num_clauses = len(self.clauses)
        with self.cnf_path.open("w", encoding="utf-8") as handle:
            handle.write(f"p cnf {num_vars} {num_clauses}\n")
            for clause in self.clauses:
                handle.write(" ".join(map(str, clause)) + " 0\n")

    @staticmethod
    def model_count_ganak(cnf_path: str | Path) -> int:
        binary = _resolve_counter_binary(GANAK_ENV_VAR, "ganak")
        result = subprocess.run(
            [binary, str(cnf_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Ganak execution failed: {result.stderr.strip()}")

        match = re.search(r"(?:s mc|c s exact arb int)\s+(\d+)", result.stdout)
        if not match:
            raise RuntimeError(f"Unable to parse Ganak output: {result.stdout.strip()}")
        return int(match.group(1))

    @staticmethod
    def model_count_approxmc(
        cnf_path: str | Path,
        epsilon: float,
        delta: float,
    ) -> int:
        binary = _resolve_counter_binary(APPROXMC_ENV_VAR, "approxmc")
        result = subprocess.run(
            [
                binary,
                f"--epsilon={epsilon}",
                f"--delta={delta}",
                f"--seed={APPROXMC_SEED}",
                str(cnf_path),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(f"ApproxMC execution failed: {result.stderr.strip()}")

        for line in result.stdout.splitlines():
            if line.startswith("s mc"):
                parts = line.split()
                if parts:
                    return int(parts[-1])

        raise RuntimeError(f"Unable to parse ApproxMC output: {result.stdout.strip()}")


def fo2_count_odd_degree(
    file_path: str | Path,
    *,
    n: int,
    m: int,
    k: int,
    counter: str = "ganak",
    epsilon: float = 0.01,
    delta: float = 0.01,
    output_dir: str | Path | None = None,
) -> int:
    """Convert odd-degree constraints to CNF and count models with external counters."""
    context = CNFContext(file_path, n=n, m=m, k=k, output_dir=output_dir)
    context.convert()
    context.dump()

    normalized_counter = counter.lower()
    if normalized_counter == "ganak":
        return CNFContext.model_count_ganak(context.cnf_path)
    if normalized_counter == "approxmc":
        return CNFContext.model_count_approxmc(context.cnf_path, epsilon, delta)

    raise ValueError(f"Unknown counter: {counter}")


def Fo2Counter(
    file_path: str | Path,
    n: int,
    m: int,
    k: int,
    counter: str = "ganak",
    epsilon: float = 0.01,
    delta: float = 0.01,
    output_dir: str | Path | None = None,
) -> int:
    """Backward-compatible alias for legacy odd-degree scripts."""
    return fo2_count_odd_degree(
        file_path=file_path,
        n=n,
        m=m,
        k=k,
        counter=counter,
        epsilon=epsilon,
        delta=delta,
        output_dir=output_dir,
    )
