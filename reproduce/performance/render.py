from __future__ import annotations

from typing import Literal
import os

from reproduce.utils.figure_naming import resolve_publish_figure_filename
from reproduce.performance.benchmarks import RuntimeConfig
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd


def plot_comparison(
    csv_file_path: str,
    config: RuntimeConfig,
    metric: str = "time_sec",
    ylabel: str = "Runtime (s)",
    output_suffix: str = "time_comparison",
) -> None:
    """
    Based on the provided CSV file and the metrics (time or memory), create a graph comparing the performance of the algorithms.
    """
    if (not os.path.exists(csv_file_path)) or os.path.getsize(csv_file_path) == 0:
        print(f"CSV is empty, so skip the drawing.: {csv_file_path}")
        return

    try:
        data = pd.read_csv(csv_file_path)
    except pd.errors.EmptyDataError:
        print(f"CSV file has no content to be parsed, so the graph will be skipped: {csv_file_path}")
        return

    data = data[~data["status"].isin(["timeout", "error"])]

    algorithms = data["algorithm"].unique()
    domain_sizes = sorted(data["domain_size"].unique())

    plt.figure(figsize=(12, 8))

    legend_labels = {
        "incremental3": "Ours",
        "dr": "Ours",
        "fast": "Fast",
        "incremental": "Incremental",
        "recursive": "Recursive",
    }

    markers = {
        "incremental3": "o",
        "dr": "o",
        "fast": "s",
        "recursive": "p",
        "incremental": "^",
    }

    colors = {
        "incremental3": "tab:orange",
        "dr": "tab:orange",
        "fast": "tab:green",
        "recursive": "tab:blue",
        "incremental": "tab:brown",
    }

    for algo in algorithms:
        algo_data = data[data["algorithm"] == algo]
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
            marker=markers.get(algo, "o"),
            label=legend_labels.get(algo, algo),
            color=colors.get(algo, "black"),
            markersize=12,
            linewidth=3,
            markeredgecolor="black",
            markeredgewidth=1,
        )

    if metric == "time_sec":
        plt.yscale("log")

    ax = plt.gca()
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    plt.xlabel("Domain Size", fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.legend(fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(True, linestyle="--", alpha=0.7)

    file_name = os.path.basename(csv_file_path)
    base_name = os.path.splitext(file_name)[0]

    if metric == "time_sec":
        section = "Section6"
    else:
        section = "AppendixB.2"

    full_result_dir = config.full_fig_results_path / section
    full_result_dir.mkdir(parents=True, exist_ok=True)

    source_publish_name = f"{base_name}_{output_suffix}.pdf"
    full_pdf_path = full_result_dir / source_publish_name
    plt.savefig(full_pdf_path, dpi=600, bbox_inches="tight")

    publish_filename = resolve_publish_figure_filename(section, source_publish_name)
    if publish_filename is not None:
        publish_result_dir = config.publish_fig_results_path / section
        publish_result_dir.mkdir(parents=True, exist_ok=True)
        publish_pdf_path = publish_result_dir / publish_filename
        plt.savefig(publish_pdf_path, dpi=600, bbox_inches="tight")
        print(f"The chart has been saved to: {publish_pdf_path}")
    else:
        print(
            "No paper figure mapping for publish target, skip reproduce/results output: "
            f"{section}/{source_publish_name}"
        )

    # png_path = os.path.join(result_dir, f"{base_name}_{output_suffix}.png")
    # plt.savefig(png_path, dpi=600, bbox_inches="tight")

    print(f"The chart has been saved to: {full_pdf_path}")
    plt.close()


def plot_metric(
    csv_file_path: str,
    config: RuntimeConfig,
    metric: Literal["time_sec", "memory_mb"] = "time_sec",
) -> None:
    mapping = {
        "time_sec": ("Runtime (s)", "time_comparison"),
        "memory_mb": ("Memory Usage (MB)", "memory_comparison"),
    }
    ylabel, suffix = mapping[metric]
    plot_comparison(csv_file_path, config, metric, ylabel, suffix)
