from __future__ import annotations

import argparse
import fnmatch
import shutil
from dataclasses import dataclass
from pathlib import Path


FIGURE_5_TO_8 = "Figure_5_to_8"
FIGURE_9 = "Figure_9"
FIGURE_10 = "Figure_10"
FIGURE_11_TO_12 = "Figure_11_to_12"
FIGURE_13 = "Figure_13"
FIGURE_14_TO_16 = "Figure_14_to_16"
TABLE_2 = "Table_2"

CANONICAL_TARGETS = [
    FIGURE_5_TO_8,
    FIGURE_9,
    FIGURE_10,
    FIGURE_11_TO_12,
    FIGURE_13,
    FIGURE_14_TO_16,
    TABLE_2,
]


@dataclass(frozen=True)
class GlobSpec:
    pattern: str
    exclude: tuple[str, ...] = ()


@dataclass(frozen=True)
class CleanupSpec:
    summary: str
    globs: tuple[GlobSpec, ...]


TARGET_CLEANUP_MAP: dict[str, CleanupSpec] = {
    FIGURE_5_TO_8: CleanupSpec(
        summary="performance main artifacts for Section6 Figure_5 to Figure_8",
        globs=(
            GlobSpec("reproduce/performance/performance_results/raw_data/results_*"),
            GlobSpec(
                "reproduce/performance/performance_results/Section6/*_time_comparison.*",
                exclude=(
                    "reproduce/performance/performance_results/Section6/m-odd-degree-graph-sc2_m*_time_comparison.*",
                    "reproduce/performance/performance_results/Section6/0mod2-regular-graph_time_comparison.*",
                    "reproduce/performance/performance_results/Section6/1mod2-regular-graph_time_comparison.*",
                    "reproduce/performance/performance_results/Section6/2mod4-regular-graph_time_comparison.*",
                ),
            ),
            GlobSpec("reproduce/results/Section6/Figure_5_*.pdf"),
            GlobSpec("reproduce/results/Section6/Figure_6_*.pdf"),
            GlobSpec("reproduce/results/Section6/Figure_7_*.pdf"),
            GlobSpec("reproduce/results/Section6/Figure_8_*.pdf"),
            GlobSpec(
                "reproduce/results/Section6/*_time_comparison.*",
                exclude=(
                    "reproduce/results/Section6/m-odd-degree-graph-sc2_m*_time_comparison.*",
                    "reproduce/results/Section6/0mod2-regular-graph_time_comparison.*",
                    "reproduce/results/Section6/1mod2-regular-graph_time_comparison.*",
                    "reproduce/results/Section6/2mod4-regular-graph_time_comparison.*",
                ),
            ),
        ),
    ),
    FIGURE_9: CleanupSpec(
        summary="performance rmodk artifacts for Section6 Figure_9",
        globs=(
            GlobSpec("reproduce/performance/performance_results/raw_data/rmodk/results_*"),
            GlobSpec("reproduce/performance/performance_results/raw_data/cnf/rmodk"),
            GlobSpec(
                "reproduce/performance/performance_results/Section6/0mod2-regular-graph_time_comparison.*",
            ),
            GlobSpec(
                "reproduce/performance/performance_results/Section6/1mod2-regular-graph_time_comparison.*",
            ),
            GlobSpec(
                "reproduce/performance/performance_results/Section6/2mod4-regular-graph_time_comparison.*",
            ),
            GlobSpec("reproduce/results/Section6/Figure_9_*.pdf"),
            GlobSpec("reproduce/results/Section6/0mod2-regular-graph_time_comparison.*"),
            GlobSpec("reproduce/results/Section6/1mod2-regular-graph_time_comparison.*"),
            GlobSpec("reproduce/results/Section6/2mod4-regular-graph_time_comparison.*"),
        ),
    ),
    FIGURE_10: CleanupSpec(
        summary="performance odd-degree artifacts for Section6 Figure_10",
        globs=(
            GlobSpec("reproduce/performance/performance_results/raw_data/odd_degree/results_*"),
            GlobSpec("reproduce/performance/performance_results/raw_data/cnf/odd_degree"),
            GlobSpec(
                "reproduce/performance/performance_results/Section6/m-odd-degree-graph-sc2_m*_time_comparison.*"
            ),
            GlobSpec("reproduce/results/Section6/Figure_10_*.pdf"),
            GlobSpec("reproduce/results/Section6/m-odd-degree-graph-sc2_m*_time_comparison.*"),
        ),
    ),
    FIGURE_11_TO_12: CleanupSpec(
        summary="correctness main artifacts for Figure_11 to Figure_12",
        globs=(
            GlobSpec("reproduce/correctness/correctness_results/raw_data/results_*"),
            GlobSpec(
                "reproduce/correctness/correctness_results/raw_data/cnf/*",
                exclude=(
                    "reproduce/correctness/correctness_results/raw_data/cnf/odd_degree",
                ),
            ),
            GlobSpec("reproduce/correctness/correctness_results/AppendixB.1/*_correctness.*"),
            GlobSpec("reproduce/results/AppendixB.1/Figure_11_*.pdf"),
            GlobSpec("reproduce/results/AppendixB.1/Figure_12_*.pdf"),
            GlobSpec("reproduce/results/AppendixB.1/*_correctness.*"),
        ),
    ),
    FIGURE_13: CleanupSpec(
        summary="correctness odd-degree artifacts for Figure_13",
        globs=(
            GlobSpec("reproduce/correctness/correctness_results/raw_data/odd_degree/results_*"),
            GlobSpec("reproduce/correctness/correctness_results/raw_data/cnf/odd_degree"),
            GlobSpec("reproduce/correctness/correctness_results/AppendixB.1/odd_degree_*.*"),
            GlobSpec("reproduce/results/AppendixB.1/Figure_13_*.pdf"),
            GlobSpec("reproduce/results/AppendixB.1/odd_degree_*.*"),
        ),
    ),
    FIGURE_14_TO_16: CleanupSpec(
        summary="performance main artifacts for AppendixB.2 Figure_14 to Figure_16",
        globs=(
            GlobSpec("reproduce/performance/performance_results/raw_data/results_*"),
            GlobSpec("reproduce/performance/performance_results/AppendixB.2/*_memory_comparison.*"),
            GlobSpec("reproduce/results/AppendixB.2/Figure_14_*.pdf"),
            GlobSpec("reproduce/results/AppendixB.2/Figure_15_*.pdf"),
            GlobSpec("reproduce/results/AppendixB.2/Figure_16_*.pdf"),
            GlobSpec("reproduce/results/AppendixB.2/*_memory_comparison.*"),
        ),
    ),
    TABLE_2: CleanupSpec(
        summary="OEISsequence outputs for Table_2",
        globs=(
            GlobSpec("reproduce/OEISsequence/odd_degree.csv"),
            GlobSpec("reproduce/results/Appendix.D/Table_2.csv"),
            GlobSpec("reproduce/OEISsequence/cnf/*.cnf"),
        ),
    ),
}

TARGET_LOG_GLOBS: dict[str, tuple[GlobSpec, ...]] = {
    FIGURE_5_TO_8: (
        GlobSpec("reproduce/logs/run_*/03_reproduce_performance_main.log"),
        GlobSpec("reproduce/logs/single_*/performance_main.log"),
        GlobSpec("reproduce/logs/single_*/performance_main.meta.json"),
    ),
    FIGURE_9: (
        GlobSpec("reproduce/logs/run_*/05_reproduce_performance_rmodk_main.log"),
        GlobSpec("reproduce/logs/single_*/performance_rmodk.log"),
        GlobSpec("reproduce/logs/single_*/performance_rmodk.meta.json"),
    ),
    FIGURE_10: (
        GlobSpec("reproduce/logs/run_*/04_reproduce_performance_odd_degree_main.log"),
        GlobSpec("reproduce/logs/single_*/performance_odd_degree.log"),
        GlobSpec("reproduce/logs/single_*/performance_odd_degree.meta.json"),
    ),
    FIGURE_11_TO_12: (
        GlobSpec("reproduce/logs/run_*/01_reproduce_correctness_main.log"),
        GlobSpec("reproduce/logs/single_*/correctness_main.log"),
        GlobSpec("reproduce/logs/single_*/correctness_main.meta.json"),
    ),
    FIGURE_13: (
        GlobSpec("reproduce/logs/run_*/02_reproduce_correctness_odd_degree_main.log"),
        GlobSpec("reproduce/logs/single_*/correctness_odd_degree.log"),
        GlobSpec("reproduce/logs/single_*/correctness_odd_degree.meta.json"),
    ),
    FIGURE_14_TO_16: (
        GlobSpec("reproduce/logs/run_*/03_reproduce_performance_main.log"),
        GlobSpec("reproduce/logs/single_*/performance_main.log"),
        GlobSpec("reproduce/logs/single_*/performance_main.meta.json"),
    ),
    TABLE_2: (
        GlobSpec("reproduce/logs/run_*/06_oeissequence_ganak_odd_degree.log"),
        GlobSpec("reproduce/logs/run_*/05_oeissequence_ganak_odd_degree.log"),
        GlobSpec("reproduce/logs/single_*/oeissequence_table2.log"),
        GlobSpec("reproduce/logs/single_*/oeissequence_table2.meta.json"),
    ),
}

ALL_TARGETS_EXTRA_LOG_GLOBS: tuple[GlobSpec, ...] = (
    GlobSpec("reproduce/logs/run_*/run_all.log"),
    GlobSpec("reproduce/logs/run_*/run_meta.json"),
)

ALIAS_TO_TARGETS: dict[str, tuple[str, ...]] = {
    # Preferred Figure/Table keys.
    "figure_5_to_8": (FIGURE_5_TO_8,),
    "figure-5-to-8": (FIGURE_5_TO_8,),
    "figure5to8": (FIGURE_5_TO_8,),
    "figure_9": (FIGURE_9,),
    "figure-9": (FIGURE_9,),
    "figure9": (FIGURE_9,),
    "figure_10": (FIGURE_10,),
    "figure-10": (FIGURE_10,),
    "figure10": (FIGURE_10,),
    "figure_11_to_12": (FIGURE_11_TO_12,),
    "figure-11-to-12": (FIGURE_11_TO_12,),
    "figure11to12": (FIGURE_11_TO_12,),
    "figure_13": (FIGURE_13,),
    "figure-13": (FIGURE_13,),
    "figure13": (FIGURE_13,),
    "figure_14_to_16": (FIGURE_14_TO_16,),
    "figure-14-to-16": (FIGURE_14_TO_16,),
    "figure14to16": (FIGURE_14_TO_16,),
    "table_2": (TABLE_2,),
    "table-2": (TABLE_2,),
    "table2": (TABLE_2,),
    "oeissequence": (TABLE_2,),
    "oesisequence": (TABLE_2,),
    # Backward-compatible aliases (deprecated in docs).
    "figure_5_to_9": (FIGURE_5_TO_8, FIGURE_9),
    "figure-5-to-9": (FIGURE_5_TO_8, FIGURE_9),
    "figure5to9": (FIGURE_5_TO_8, FIGURE_9),
    "correctness": (FIGURE_11_TO_12,),
    "correctness.main": (FIGURE_11_TO_12,),
    "reproduce.correctness.main": (FIGURE_11_TO_12,),
    "correctness-odd-degree": (FIGURE_13,),
    "correctness_odd_degree": (FIGURE_13,),
    "correctness.odd_degree": (FIGURE_13,),
    "correctness.odd_degree.main": (FIGURE_13,),
    "reproduce.correctness.odd_degree.main": (FIGURE_13,),
    "performance": (FIGURE_5_TO_8, FIGURE_9, FIGURE_14_TO_16),
    "performance.main": (FIGURE_5_TO_8, FIGURE_14_TO_16),
    "reproduce.performance.main": (FIGURE_5_TO_8, FIGURE_14_TO_16),
    "performance-rmodk": (FIGURE_9,),
    "performance_rmodk": (FIGURE_9,),
    "performance.rmodk": (FIGURE_9,),
    "performance.rmodk.main": (FIGURE_9,),
    "reproduce.performance.rmodk.main": (FIGURE_9,),
    "performance-odd-degree": (FIGURE_10,),
    "performance_odd_degree": (FIGURE_10,),
    "performance.odd_degree": (FIGURE_10,),
    "performance.odd_degree.main": (FIGURE_10,),
    "reproduce.performance.odd_degree.main": (FIGURE_10,),
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Delete generated outputs by paper artifact keys "
            "(Figure_X_to_Y / Table_2)."
        )
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help=(
            "Delete outputs for all supported Figure/Table targets and "
            "corresponding logs under reproduce/logs."
        ),
    )
    parser.add_argument(
        "--all-figure",
        action="store_true",
        help="Delete outputs for all supported Figure/Table targets.",
    )
    parser.add_argument(
        "--figure",
        action="append",
        default=[],
        metavar="FIGURE_OR_TABLE",
        help=(
            "Figure/Table key to delete. Repeatable. "
            "Examples: Figure_5_to_8, Figure_9, Figure_10, "
            "Figure_11_to_12, Figure_13, "
            "Figure_14_to_16, Table_2"
        ),
    )
    # Keep old CLI compatible, but hide it from help and docs.
    parser.add_argument(
        "--main",
        action="append",
        default=[],
        metavar="MAIN",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print matched paths but do not delete anything.",
    )
    parser.add_argument(
        "--logs",
        action="store_true",
        help=(
            "Also delete corresponding logs under reproduce/logs "
            "(step logs from run_*/ and single-run logs from single_*/). "
            "Note: --all already includes logs."
        ),
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Skip interactive confirmation.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List supported Figure/Table keys and exit.",
    )
    return parser


def _print_supported_targets() -> None:
    print("Supported Figure/Table keys:")
    print(f"  - {FIGURE_5_TO_8}: performance main (Section6)")
    print(f"  - {FIGURE_9}: performance rmodk")
    print(f"  - {FIGURE_10}: performance odd-degree")
    print(f"  - {FIGURE_11_TO_12}: correctness main")
    print(f"  - {FIGURE_13}: correctness odd-degree")
    print(f"  - {FIGURE_14_TO_16}: performance main (AppendixB.2)")
    print(f"  - {TABLE_2}: OEISsequence output")
    print()
    print("Module to artifact mapping:")
    print("  - reproduce.performance.main -> Figure_5_to_8 + Figure_14_to_16")
    print("  - reproduce.performance.rmodk.main -> Figure_9")
    print("  - reproduce.performance.odd_degree.main -> Figure_10")
    print("  - reproduce.correctness.main -> Figure_11_to_12")
    print("  - reproduce.correctness.odd_degree.main -> Figure_13")
    print("  - reproduce/OEISsequence/ganak_odd_degree.py -> Table_2")
    print()
    print("CLI shortcuts:")
    print("  - --all-figure: all Figure/Table outputs")
    print("  - --all: all Figure/Table outputs + corresponding logs")


def _resolve_selected_targets(args: argparse.Namespace) -> list[str]:
    if args.all or args.all_figure:
        return list(CANONICAL_TARGETS)

    selected: list[str] = []
    unknown: list[str] = []
    raw_keys = list(args.figure) + list(args.main)

    for item in raw_keys:
        key = item.strip().lower()
        targets = ALIAS_TO_TARGETS.get(key)
        if targets is None:
            unknown.append(item)
            continue
        for target in targets:
            if target not in selected:
                selected.append(target)

    if unknown:
        raise ValueError(
            "Unknown Figure/Table key(s): "
            + ", ".join(unknown)
            + ". Use --list to see supported keys."
        )

    if not selected:
        raise ValueError(
            "Please provide --all, --all-figure, or at least one --figure value."
        )

    return selected


def _is_excluded(relative_path: str, excludes: tuple[str, ...]) -> bool:
    for pattern in excludes:
        if fnmatch.fnmatch(relative_path, pattern):
            return True
    return False


def _add_glob_matches(repo_root: Path, targets: set[Path], glob_spec: GlobSpec) -> None:
    for path in repo_root.glob(glob_spec.pattern):
        relative_path = path.relative_to(repo_root).as_posix()
        if _is_excluded(relative_path, glob_spec.exclude):
            continue
        targets.add(path)


def _collect_targets(
    repo_root: Path,
    selected_targets: list[str],
    include_logs: bool,
) -> list[Path]:
    targets: set[Path] = set()

    for artifact_key in selected_targets:
        spec = TARGET_CLEANUP_MAP[artifact_key]
        for glob_spec in spec.globs:
            _add_glob_matches(repo_root=repo_root, targets=targets, glob_spec=glob_spec)

        if include_logs:
            for log_glob in TARGET_LOG_GLOBS.get(artifact_key, ()):
                _add_glob_matches(repo_root=repo_root, targets=targets, glob_spec=log_glob)

    if include_logs and set(selected_targets) == set(CANONICAL_TARGETS):
        for log_glob in ALL_TARGETS_EXTRA_LOG_GLOBS:
            _add_glob_matches(repo_root=repo_root, targets=targets, glob_spec=log_glob)

    # Delete deeper paths first to avoid parent-child ordering issues.
    return sorted(targets, key=lambda p: (len(p.parts), p.as_posix()), reverse=True)


def _print_plan(repo_root: Path, selected_targets: list[str], targets: list[Path]) -> None:
    print(f"Repository root: {repo_root}")
    print("Selected Figure/Table targets:")
    for key in selected_targets:
        print(f"  - {key}: {TARGET_CLEANUP_MAP[key].summary}")

    print(f"Matched paths: {len(targets)}")
    for path in targets:
        print(f"  - {path.relative_to(repo_root).as_posix()}")


def _confirm(total_paths: int) -> bool:
    answer = input(f"Delete {total_paths} matched path(s)? [y/N]: ").strip().lower()
    return answer in {"y", "yes"}


def _delete_path(path: Path, dry_run: bool) -> bool:
    if dry_run:
        return True

    try:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()
        return True
    except FileNotFoundError:
        return False


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list:
        _print_supported_targets()
        return 0

    try:
        selected_targets = _resolve_selected_targets(args)
    except ValueError as exc:
        parser.error(str(exc))

    repo_root = Path(__file__).resolve().parents[2]
    include_logs = args.logs or args.all
    targets = _collect_targets(
        repo_root=repo_root,
        selected_targets=selected_targets,
        include_logs=include_logs,
    )

    _print_plan(repo_root=repo_root, selected_targets=selected_targets, targets=targets)

    if not targets:
        print("No matching outputs found. Nothing to delete.")
        return 0

    if not args.dry_run and not args.yes and not _confirm(len(targets)):
        print("Cancelled. No files were deleted.")
        return 0

    removed = 0
    missing = 0

    for path in targets:
        if not path.exists() and not args.dry_run:
            missing += 1
            continue

        if _delete_path(path=path, dry_run=args.dry_run):
            removed += 1
        else:
            missing += 1

    if args.dry_run:
        print(f"Dry run complete. {removed} path(s) would be deleted.")
    else:
        print(f"Deletion complete. Removed {removed} path(s), missing {missing} path(s).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
