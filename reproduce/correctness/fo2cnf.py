from __future__ import annotations

import copy
import os
import re
import shutil
import subprocess
from itertools import combinations, product
from pathlib import Path
from typing import Dict, Set

import sympy
from pysat.card import CardEnc
from pysat.formula import CNF
from pysat.solvers import Solver

from wfomc import parse_input
from wfomc.fol.syntax import AtomicFormula, Const, Pred, QFFormula, X, Y, top

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
    """Convert a WFOMC FO sentence into CNF and run model counting backends."""

    def __init__(
        self,
        file_path: str | Path,
        domain_size: int,
        output_dir: str | Path | None = None,
    ) -> None:
        self.file_path = Path(file_path)
        self.file_name = self.file_path.name
        self.file_dir = Path(output_dir) if output_dir else self.file_path.parent
        self.file_dir.mkdir(parents=True, exist_ok=True)

        self.problem = parse_input(str(self.file_path))
        self.problem.domain = {Const(str(i)) for i in range(domain_size)}
        self.domain = self.problem.domain

        self.expr = sympy.true
        self.atom_to_id: Dict[AtomicFormula, int] = {}
        self.sym_to_id: Dict[sympy.Symbol, int] = {}
        self.next_var_id = 1

        cnf_file_name = f"{self.file_path.stem}_domain_size_{domain_size}.cnf"
        self.cnf_path = self.file_dir / cnf_file_name
        self.clauses: list[list[int]] = []

    def convert(self) -> None:
        """Main conversion pipeline."""
        self._ground_universal()
        if self.expr is not sympy.false:
            self._ground_extension()
            self._ground_counting()
            if self._check_exists_linear_order():
                self._add_linear_order_axioms()

        self._add_clauses_from_expr()

        if (
            self.problem.cardinality_constraint
            and not self.problem.cardinality_constraint.empty()
        ):
            self._ground_cardinality_constraints()

    def _register_atom(self, atom: AtomicFormula) -> None:
        if atom not in self.atom_to_id:
            self.atom_to_id[atom] = self.next_var_id
            self.sym_to_id[atom.expr] = self.next_var_id
            self.next_var_id += 1

    def _add_clauses_from_expr(self) -> None:
        if self.expr is sympy.true:
            return
        if self.expr is sympy.false:
            self.clauses.append([])
            return

        cnf_expr = sympy.to_cnf(self.expr)
        expr_clauses = cnf_expr.args if isinstance(cnf_expr, sympy.And) else [cnf_expr]

        for cl in expr_clauses:
            lits: list[int] = []
            if isinstance(cl, (sympy.Symbol, sympy.Not)):
                atoms = [cl]
            elif isinstance(cl, sympy.Or):
                atoms = list(cl.args)
            else:
                continue

            for atom_expr in atoms:
                if isinstance(atom_expr, sympy.Symbol):
                    if atom_expr in self.sym_to_id:
                        lits.append(self.sym_to_id[atom_expr])
                elif isinstance(atom_expr, sympy.Not):
                    base = atom_expr.args[0]
                    if base in self.sym_to_id:
                        lits.append(-self.sym_to_id[base])

            if lits:
                self.clauses.append(lits)

    def _extract_qf(self, formula) -> QFFormula:
        while not isinstance(formula, QFFormula):
            formula = formula.quantified_formula
        return formula

    def _ground_universal(self) -> None:
        uni_qf = self._extract_qf(copy.deepcopy(self.problem.sentence.uni_formula))

        for e1, e2 in product(self.domain, repeat=2):
            grounded = uni_qf.substitute({X: e1, Y: e2})
            if grounded is top:
                continue

            if grounded.expr is None:
                self.expr = sympy.false
                return

            for atom in grounded.atoms():
                self._register_atom(atom)
            self.expr &= grounded.expr

    def _ground_extension(self) -> None:
        ext_qfs = [
            self._extract_qf(copy.deepcopy(formula.quantified_formula))
            for formula in self.problem.sentence.ext_formulas
        ]
        if not ext_qfs:
            return

        for e1 in self.domain:
            for ext_qf in ext_qfs:
                disjunction = sympy.false
                for e2 in self.domain:
                    grounded = ext_qf.substitute({X: e1, Y: e2})
                    for atom in grounded.atoms():
                        self._register_atom(atom)
                    if grounded.expr is not None:
                        disjunction |= grounded.expr
                self.expr &= disjunction

    def _ground_cardinality_constraints(self) -> None:
        for pred_map, op, bound in self.problem.cardinality_constraint.constraints:
            for pred, _ in pred_map.items():
                pred_name = str(pred)
                bound_int = int(bound)

                related_vars = [
                    var_id
                    for sym, var_id in self.sym_to_id.items()
                    if pred_name in str(sym)
                ]
                max_count = len(related_vars)

                # Convert out-of-range constraints into SAT/UNSAT directly,
                # instead of letting CardEnc raise "Wrong bound".
                if op == "=":
                    if bound_int < 0 or bound_int > max_count:
                        self.clauses.append([])
                        continue
                    if max_count == 0 and bound_int == 0:
                        continue

                elif op == "<=":
                    if bound_int < 0:
                        self.clauses.append([])
                        continue
                    if bound_int >= max_count:
                        continue
                    if max_count == 0:
                        continue

                elif op == ">=":
                    if bound_int <= 0:
                        continue
                    if bound_int > max_count:
                        self.clauses.append([])
                        continue
                    if max_count == 0:
                        self.clauses.append([])
                        continue

                else:
                    raise RuntimeError(f"Unknown operator: {op}")

                if op == "<=":
                    card_cnf = CardEnc.atmost(
                        lits=related_vars,
                        bound=bound_int,
                        top_id=self.next_var_id - 1,
                    )
                elif op == ">=":
                    card_cnf = CardEnc.atleast(
                        lits=related_vars,
                        bound=bound_int,
                        top_id=self.next_var_id - 1,
                    )
                elif op == "=":
                    card_cnf = CardEnc.equals(
                        lits=related_vars,
                        bound=bound_int,
                        top_id=self.next_var_id - 1,
                    )

                self.clauses.extend(card_cnf.clauses)
                self.next_var_id = card_cnf.nv + 1

    def _is_single_layer(self, formula) -> bool:
        return isinstance(formula.quantified_formula, QFFormula)

    def _ground_counting(self) -> None:
        cnt_formulas = self.problem.sentence.cnt_formulas
        single_layer_formulas: list = []
        double_layer_formulas: list = []

        for formula in cnt_formulas:
            if self._is_single_layer(formula):
                single_layer_formulas.append(formula)
            else:
                double_layer_formulas.append(formula)

        for formula in single_layer_formulas:
            scope = formula.quantifier_scope
            if scope.comparator not in ("=", "mod"):
                raise RuntimeError(
                    "Unsupported comparator in counting quantifier: "
                    f"{scope.comparator}"
                )

            if scope.comparator == "=":
                self.expr &= self._build_eq(
                    formula.quantified_formula,
                    None,
                    scope.quantified_var,
                    None,
                    scope.count_param,
                )
            else:
                self.expr &= self._build_mod(
                    formula.quantified_formula,
                    None,
                    scope.quantified_var,
                    None,
                    scope.count_param,
                )

        for e1 in self.domain:
            for formula in double_layer_formulas:
                inner_scope = formula.quantified_formula.quantifier_scope
                inner_qf = formula.quantified_formula.quantified_formula
                free_vars = inner_qf.vars() - {inner_scope.quantified_var}
                var_x = next(iter(free_vars)) if free_vars else None

                if inner_scope.comparator not in ("=", "mod"):
                    raise RuntimeError(
                        "Unsupported comparator in counting quantifier: "
                        f"{inner_scope.comparator}"
                    )

                if inner_scope.comparator == "=":
                    self.expr &= self._build_eq(
                        inner_qf,
                        var_x,
                        inner_scope.quantified_var,
                        e1,
                        inner_scope.count_param,
                    )
                else:
                    self.expr &= self._build_mod(
                        inner_qf,
                        var_x,
                        inner_scope.quantified_var,
                        e1,
                        inner_scope.count_param,
                    )

    def _build_eq(self, inner_qf, var_x, var_y, e1, k: int) -> sympy.Expr:
        clause_expr = sympy.false
        for selected_y in combinations(self.domain, k):
            sub_conj = sympy.true
            domain_set = set(self.domain)
            selected_y_set = set(selected_y)

            for y in selected_y_set:
                subst = {var_y: y, **({var_x: e1} if var_x else {})}
                grounded = inner_qf.substitute(subst)
                for atom in grounded.atoms():
                    self._register_atom(atom)
                sub_conj &= grounded.expr

            for y in domain_set - selected_y_set:
                subst = {var_y: y, **({var_x: e1} if var_x else {})}
                grounded = inner_qf.substitute(subst)
                for atom in grounded.atoms():
                    self._register_atom(atom)
                sub_conj &= ~grounded.expr

            clause_expr |= sub_conj

        return sympy.simplify_logic(clause_expr, form="cnf")

    def _build_mod(self, inner_qf, var_x, var_y, e1, rk) -> sympy.Expr:
        remainder, mod_base = rk
        clause_expr = sympy.false

        for count in range(len(self.domain) + 1):
            if count % mod_base != remainder:
                continue

            for selected_y in combinations(self.domain, count):
                sub_conj = sympy.true
                domain_set = set(self.domain)
                selected_y_set = set(selected_y)

                for y in selected_y_set:
                    subst = {var_y: y, **({var_x: e1} if var_x else {})}
                    grounded = inner_qf.substitute(subst)
                    for atom in grounded.atoms():
                        self._register_atom(atom)
                    sub_conj &= grounded.expr

                for y in domain_set - selected_y_set:
                    subst = {var_y: y, **({var_x: e1} if var_x else {})}
                    grounded = inner_qf.substitute(subst)
                    for atom in grounded.atoms():
                        self._register_atom(atom)
                    sub_conj &= ~grounded.expr

                clause_expr |= sub_conj

        return sympy.simplify_logic(clause_expr, form="cnf")

    def dump(self) -> None:
        num_vars = self.next_var_id - 1
        num_clauses = len(self.clauses)
        with self.cnf_path.open("w", encoding="utf-8") as handle:
            handle.write(f"p cnf {num_vars} {num_clauses}\n")
            for clause in self.clauses:
                handle.write(" ".join(map(str, clause)) + " 0\n")

    def _check_exists_linear_order(self) -> bool:
        return "LEQ" in [pred.name for pred in self._collect_all_predicates()]

    def _add_linear_order_axioms(self) -> None:
        axioms_expr = sympy.true
        leq_pred = Pred("LEQ", 2)

        for c1, c2 in product(self.domain, repeat=2):
            self._register_atom(AtomicFormula(leq_pred, (c1, c2), True))

        for x in self.domain:
            axioms_expr &= AtomicFormula(leq_pred, (x, x), True).expr

        for x, y in combinations(self.domain, 2):
            axioms_expr &= sympy.Or(
                AtomicFormula(leq_pred, (x, y), False).expr,
                AtomicFormula(leq_pred, (y, x), False).expr,
            )

        for x, y, z in product(self.domain, repeat=3):
            axioms_expr &= sympy.Or(
                AtomicFormula(leq_pred, (x, y), False).expr,
                AtomicFormula(leq_pred, (y, z), False).expr,
                AtomicFormula(leq_pred, (x, z), True).expr,
            )

        for x, y in product(self.domain, repeat=2):
            axioms_expr &= sympy.Or(
                AtomicFormula(leq_pred, (x, y), True).expr,
                AtomicFormula(leq_pred, (y, x), True).expr,
            )

        self.expr &= axioms_expr

    def _collect_all_predicates(self) -> Set[Pred]:
        all_preds: Set[Pred] = set()
        formulas_to_check = []

        if self.problem.sentence.uni_formula:
            formulas_to_check.append(self.problem.sentence.uni_formula)
        formulas_to_check.extend(self.problem.sentence.ext_formulas)
        formulas_to_check.extend(self.problem.sentence.cnt_formulas)

        for formula in formulas_to_check:
            qf_formula = self._extract_qf(copy.deepcopy(formula))
            for atom in qf_formula.atoms():
                all_preds.add(atom.pred)

        if (
            self.problem.cardinality_constraint
            and not self.problem.cardinality_constraint.empty()
        ):
            for pred_map, _, _ in self.problem.cardinality_constraint.constraints:
                for pred in pred_map.keys():
                    all_preds.add(pred)

        if self.problem.unary_evidence:
            for pred in self.problem.unary_evidence.keys():
                all_preds.add(pred)

        return all_preds

    @staticmethod
    def model_count_pysat(cnf_path: str | Path) -> int:
        cnf = CNF(from_file=str(cnf_path))
        count = 0
        with Solver(bootstrap_with=cnf) as solver:
            for _ in solver.enum_models():
                count += 1
        return count

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
            raise RuntimeError(
                f"ApproxMC execution failed: {result.stderr.strip()}"
            )

        for line in result.stdout.splitlines():
            if line.startswith("s mc"):
                return int(line.split()[-1])

        raise RuntimeError(f"Unable to parse ApproxMC output: {result.stdout.strip()}")

    def exists_cnf_file(self) -> Path | None:
        if self.cnf_path.exists():
            return self.cnf_path
        return None


def fo2_count(
    file_path: str | Path,
    *,
    domain_size: int,
    counter: str = "ganak",
    epsilon: float = 0.01,
    delta: float = 0.01,
    output_dir: str | Path | None = None,
) -> int:
    """Convert one model to CNF and count models with a selected backend."""
    context = CNFContext(file_path, domain_size=domain_size, output_dir=output_dir)

    if not context.exists_cnf_file():
        context.convert()
        context.dump()

    normalized_counter = counter.lower()
    if normalized_counter == "pysat":
        return CNFContext.model_count_pysat(context.cnf_path)
    if normalized_counter == "ganak":
        return CNFContext.model_count_ganak(context.cnf_path)
    if normalized_counter == "approxmc":
        return CNFContext.model_count_approxmc(context.cnf_path, epsilon, delta)

    raise ValueError(f"Unknown counter: {counter}")


def Fo2Counter(
    file_path: str | Path,
    domain_size: int,
    counter: str = "ganak",
    epsilon: float = 0.01,
    delta: float = 0.01,
    output_dir: str | Path | None = None,
) -> int:
    """Backward-compatible alias for legacy check scripts."""
    return fo2_count(
        file_path=file_path,
        domain_size=domain_size,
        counter=counter,
        epsilon=epsilon,
        delta=delta,
        output_dir=output_dir,
    )
