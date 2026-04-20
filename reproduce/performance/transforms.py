from __future__ import annotations


def apply_benchmark_transforms(problem, model_name: str, domain_size: int) -> None:
    """Apply benchmark-specific mutations to parsed problems."""
    # This special-case transform is only for the BA benchmark where |Eq| = n.
    if model_name == "BA_CC":
        if (
            problem.cardinality_constraint
            and problem.cardinality_constraint.constraints
        ):
            new_constraints = []
            for old_constraint in problem.cardinality_constraint.constraints:
                new_constraint = (
                    old_constraint[0],
                    old_constraint[1],
                    domain_size,  # Assume the cardinality bound equals domain_size.
                )
                new_constraints.append(new_constraint)
            problem.cardinality_constraint.constraints = new_constraints
