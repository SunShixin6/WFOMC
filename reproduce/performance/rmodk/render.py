from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from reproduce.utils.figure_naming import resolve_publish_figure_filename

from .benchmarks import RuntimeConfig

PAPER_ALGORITHM_ORDER = ["incremental3", "ganak", "approxmc"]


def plot_comparison(
    csv_file_path: str,
    config: RuntimeConfig,
    metric: str = "time_sec",
    ylabel: str = "Runtime (s)",
    output_suffix: str = "time_comparison",
) -> None:
    """Draw Figure 9 runtime comparison charts from one rmodk CSV file."""
    csv_path = Path(csv_file_path)
    if (not csv_path.exists()) or csv_path.stat().st_size == 0:
        print(f"CSV is empty, skip plotting: {csv_path}")
        return

    try:
        data = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        print(f"CSV has no parsable content, skip plotting: {csv_path}")
        return

    valid_data = data[
        ~data["status"].astype(str).str.contains("timeout|error", case=False, na=False)
    ]
    if valid_data.empty:
        print(f"CSV has no successful rows, skip plotting: {csv_path}")
        return

    available_algorithms = [str(a).lower() for a in valid_data["algorithm"].dropna().unique()]
    algorithms = [a for a in PAPER_ALGORITHM_ORDER if a in available_algorithms]
    algorithms.extend(a for a in available_algorithms if a not in algorithms)
    domain_sizes = sorted(valid_data["domain_size"].unique())

    plt.figure(figsize=(12, 8))

    legend_labels = {
        "incremental3": "Ours",
        "ganak": "Ganak",
        "approxmc": "ApproxMC",
    }
    markers = {
        "incremental3": "o",
        "ganak": "s",
        "approxmc": "p",
    }
    colors = {
        "incremental3": "tab:orange",
        "ganak": "tab:green",
        "approxmc": "tab:blue",
    }

    for algorithm in algorithms:
        algo_data = valid_data[valid_data["algorithm"] == algorithm]
        values = [
            (
                algo_data[algo_data["domain_size"] == size][metric].iloc[0]
                if not algo_data[algo_data["domain_size"] == size].empty
                else float("nan")
            )
            for size in domain_sizes
        ]

        plt.plot(
            domain_sizes,
            values,
            marker=markers.get(algorithm, "o"),
            label=legend_labels.get(algorithm, algorithm),
            color=colors.get(algorithm, "black"),
            markersize=10,
            linewidth=2.5,
            markeredgecolor="black",
            markeredgewidth=1,
        )

    if metric == "time_sec":
        plt.yscale("log")

    ax = plt.gca()
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    plt.xlabel("Domain Size", fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.7)

    file_name = f"{csv_path.stem}_{output_suffix}.pdf"
    full_pdf_path = config.full_fig_results_root / file_name
    full_pdf_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(full_pdf_path, dpi=600, bbox_inches="tight")

    publish_file_name = resolve_publish_figure_filename("Section6", file_name)
    if publish_file_name is not None:
        publish_pdf_path = config.publish_fig_results_root / publish_file_name
        publish_pdf_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(publish_pdf_path, dpi=600, bbox_inches="tight")
        print(f"Chart saved to: {publish_pdf_path}")
    else:
        print(
            "No paper figure mapping for publish target, skip reproduce/results output: "
            f"Section6/{file_name}"
        )

    print(f"Chart saved to: {full_pdf_path}")
    plt.close()


def plot_metric(
    csv_file_path: str,
    config: RuntimeConfig,
    metric: str = "time_sec",
) -> None:
    if metric != "time_sec":
        raise ValueError("rmodk render only supports metric='time_sec'.")
    plot_comparison(
        csv_file_path,
        config=config,
        metric="time_sec",
        ylabel="Runtime (s)",
        output_suffix="time_comparison",
    )
