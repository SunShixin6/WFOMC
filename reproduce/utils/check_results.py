from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class CheckItem:
    name: str
    ok: bool
    detail: str


def _exists_nonempty_file(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, f"missing: {path}"
    if not path.is_file():
        return False, f"not a file: {path}"
    if path.stat().st_size <= 0:
        return False, f"empty file: {path}"
    return True, f"ok: {path}"


def _has_any_pdf(path: Path) -> tuple[bool, str]:
    if not path.exists() or not path.is_dir():
        return False, f"missing directory: {path}"
    count = sum(1 for _ in path.glob("*.pdf"))
    if count == 0:
        return False, f"no pdf files found under: {path}"
    return True, f"ok: {count} pdf(s) in {path}"


def _has_any_csv_recursive(path: Path) -> tuple[bool, str]:
    if not path.exists() or not path.is_dir():
        return False, f"missing directory: {path}"
    count = sum(1 for _ in path.rglob("*.csv"))
    if count == 0:
        return False, f"no csv files found under: {path}"
    return True, f"ok: {count} csv file(s) in {path}"


def _expected_strict_figure_files() -> dict[str, list[str]]:
    return {
        "Section6": [
            "Figure_5_a.pdf",
            "Figure_5_b.pdf",
            "Figure_5_c.pdf",
            "Figure_6_a.pdf",
            "Figure_6_b.pdf",
            "Figure_6_c.pdf",
            "Figure_6_d.pdf",
            "Figure_6_e.pdf",
            "Figure_6_f.pdf",
            "Figure_7_a.pdf",
            "Figure_7_b.pdf",
            "Figure_8_a.pdf",
            "Figure_8_b.pdf",
            "Figure_9_a.pdf",
            "Figure_9_b.pdf",
            "Figure_9_c.pdf",
            "Figure_10_a.pdf",
            "Figure_10_b.pdf",
            "Figure_10_c.pdf",
        ],
        "AppendixB.1": [
            "Figure_11_a.pdf",
            "Figure_11_b.pdf",
            "Figure_11_c.pdf",
            "Figure_12_a.pdf",
            "Figure_12_b.pdf",
            "Figure_12_c.pdf",
            "Figure_13_a.pdf",
            "Figure_13_b.pdf",
            "Figure_13_c.pdf",
        ],
        "AppendixB.2": [
            "Figure_14_a.pdf",
            "Figure_14_b.pdf",
            "Figure_14_c.pdf",
            "Figure_15_a.pdf",
            "Figure_15_b.pdf",
            "Figure_15_c.pdf",
            "Figure_16_a.pdf",
            "Figure_16_b.pdf",
        ],
    }


def run_checks(repo_root: Path, strict: bool) -> list[CheckItem]:
    checks: list[CheckItem] = []

    table2 = repo_root / "reproduce" / "results" / "Appendix.D" / "Table_2.csv"
    ok, detail = _exists_nonempty_file(table2)
    checks.append(CheckItem("table2_nonempty", ok, detail))

    section6 = repo_root / "reproduce" / "results" / "Section6"
    appendix_b1 = repo_root / "reproduce" / "results" / "AppendixB.1"
    appendix_b2 = repo_root / "reproduce" / "results" / "AppendixB.2"

    for name, path in [
        ("section6_has_pdf", section6),
        ("appendix_b1_has_pdf", appendix_b1),
        ("appendix_b2_has_pdf", appendix_b2),
    ]:
        ok, detail = _has_any_pdf(path)
        checks.append(CheckItem(name, ok, detail))

    if strict:
        strict_expected = _expected_strict_figure_files()
        for section, filenames in strict_expected.items():
            section_dir = repo_root / "reproduce" / "results" / section
            missing = [name for name in filenames if not (section_dir / name).is_file()]
            if missing:
                checks.append(
                    CheckItem(
                        f"strict_expected_{section}",
                        False,
                        f"missing {len(missing)} expected file(s) in {section_dir}: {', '.join(missing[:8])}" +
                        (" ..." if len(missing) > 8 else ""),
                    )
                )
            else:
                checks.append(
                    CheckItem(
                        f"strict_expected_{section}",
                        True,
                        f"ok: all {len(filenames)} expected files present in {section_dir}",
                    )
                )

        csv_roots = [
            repo_root / "reproduce" / "correctness" / "correctness_results" / "raw_data",
            repo_root / "reproduce" / "correctness" / "correctness_results" / "raw_data" / "odd_degree",
            repo_root / "reproduce" / "performance" / "performance_results" / "raw_data",
            repo_root / "reproduce" / "performance" / "performance_results" / "raw_data" / "odd_degree",
        ]
        for path in csv_roots:
            ok, detail = _has_any_csv_recursive(path)
            checks.append(
                CheckItem(
                    f"strict_csv_present::{path.relative_to(repo_root).as_posix()}",
                    ok,
                    detail,
                )
            )

    return checks


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check whether reproduction outputs are present and complete."
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Enable strict checks for expected figure filenames and CSV coverage.",
    )
    parser.add_argument(
        "--json",
        type=str,
        default=None,
        metavar="PATH",
        help="Write a machine-readable JSON report to PATH.",
    )
    parser.add_argument(
        "--repo-root",
        type=str,
        default=None,
        help="Repository root path (defaults to auto-detected root).",
    )
    return parser


def _print_report(checks: list[CheckItem], strict: bool) -> None:
    mode = "STRICT" if strict else "BASIC"
    print(f"Reproduction output check ({mode})")
    print("-" * 60)

    failed = 0
    for item in checks:
        status = "PASS" if item.ok else "FAIL"
        print(f"[{status}] {item.name}: {item.detail}")
        if not item.ok:
            failed += 1

    print("-" * 60)
    print(f"Summary: {len(checks) - failed}/{len(checks)} checks passed")
    if failed == 0:
        print("Overall: PASS")
    else:
        print("Overall: FAIL")


def _write_json_report(path: Path, repo_root: Path, strict: bool, checks: list[CheckItem]) -> None:
    payload: dict[str, Any] = {
        "repo_root": str(repo_root),
        "strict": strict,
        "overall_ok": all(item.ok for item in checks),
        "checks": [asdict(item) for item in checks],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.repo_root:
        repo_root = Path(args.repo_root).resolve()
    else:
        repo_root = Path(__file__).resolve().parents[2]

    checks = run_checks(repo_root=repo_root, strict=args.strict)
    _print_report(checks, strict=args.strict)

    if args.json:
        _write_json_report(Path(args.json).resolve(), repo_root=repo_root, strict=args.strict, checks=checks)

    return 0 if all(item.ok for item in checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
