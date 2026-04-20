from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from reproduce.utils.figure_naming import resolve_publish_figure_filename

from .benchmarks import RuntimeConfig


# Paper-aligned y-axis ranges for Figure13 subplots.
_PAPER_Y_LIMITS_BY_SWEEP: dict[str, tuple[float, float]] = {
    "vary_n": (1e-1, 1e2),
    "vary_k": (1e-1, 1e3),
    "vary_m": (0.0, 1e5),
}


def _paper_y_limits(sweep_name: str | None, variable_name: str) -> tuple[float, float] | None:
    if sweep_name is not None:
        y_limits = _PAPER_Y_LIMITS_BY_SWEEP.get(sweep_name)
        if y_limits is not None:
            return y_limits

    # Fallback when older CSVs do not contain `sweep`.
    fallback_by_variable = {
        "domain_size": (1e-1, 1e2),
        "k_value": (1e-1, 1e3),
        "m_value": (0.0, 1e5),
    }
    return fallback_by_variable.get(variable_name)


def _paper_x_label(sweep_name: str | None, variable_name: str) -> str:
    if sweep_name == "vary_n":
        return "Domain Size"
    if sweep_name == "vary_k":
        return "Number of Edges"
    if sweep_name == "vary_m":
        return "Number of odd-degree vertices"

    # Fallback for older CSVs without `sweep`.
    fallback_by_variable = {
        "domain_size": "Domain Size",
        "k_value": "Number of Edges",
        "m_value": "Number of odd-degree vertices",
    }
    return fallback_by_variable.get(variable_name, variable_name)


def _paper_pdf_basename(sweep_name: str | None, variable_name: str) -> str:
    if sweep_name == "vary_n":
        return "odd_degree_n"
    if sweep_name == "vary_k":
        return "odd_degree_k"
    if sweep_name == "vary_m":
        return "odd_degree_m"

    # Fallback for older CSVs without `sweep`.
    fallback_by_variable = {
        "domain_size": "odd_degree_n",
        "k_value": "odd_degree_k",
        "m_value": "odd_degree_m",
    }
    return fallback_by_variable.get(variable_name, csv_safe_name(variable_name))


def csv_safe_name(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in name)


def plot_model_count(csv_file_path: str, config: RuntimeConfig) -> None:
    """Plot one Figure13 CSV as model-count comparison across algorithms."""
    csv_path = Path(csv_file_path)
    if (not csv_path.exists()) or csv_path.stat().st_size == 0:
        print(f"CSV is empty, skip plotting: {csv_path}")
        return

    try:
        data = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        print(f"CSV has no parsable content, skip plotting: {csv_path}")
        return

    required_columns = {
        "algorithm",
        "variable_name",
        "variable_value",
        "model_count",
        "status",
    }
    if not required_columns.issubset(data.columns):
        missing = sorted(required_columns.difference(data.columns))
        print(f"CSV missing required columns {missing}, skip plotting: {csv_path}")
        return

    data = data[
        ~data["status"].astype(str).str.contains("timeout|error", case=False, na=False)
    ]
    if data.empty:
        print(f"CSV has no successful rows, skip plotting: {csv_path}")
        return

    data = data.copy()
    data["variable_value"] = pd.to_numeric(data["variable_value"], errors="coerce")
    data["model_count"] = pd.to_numeric(data["model_count"], errors="coerce")
    # Keep only numeric rows here; zero filtering is applied later for plotting.
    data = data.dropna(subset=["variable_value", "model_count"])
    if data.empty:
        print(f"CSV has no numeric rows, skip plotting: {csv_path}")
        return

    # Preserve full x-axis range, including domains that become empty after zero filtering.
    all_x_values = data["variable_value"].copy()

    # Do not render model_count == 0 points.
    plot_data = data[data["model_count"] > 0].copy()
    if plot_data.empty:
        print(f"All model counts are zero after filtering, skip plotting: {csv_path}")
        return

    variable_name = str(data["variable_name"].iloc[0])
    sweep_name = str(data["sweep"].iloc[0]) if "sweep" in data.columns else None
    plt.figure(figsize=(12, 8))

    ours = plot_data[
        plot_data["algorithm"].astype(str).str.lower().isin({"incremental3", "dr"})
    ].sort_values("variable_value")
    ganak = plot_data[plot_data["algorithm"].astype(str).str.lower() == "ganak"].sort_values("variable_value")
    approx = plot_data[plot_data["algorithm"].astype(str).str.lower() == "approxmc"].sort_values("variable_value")

    found_any = not (ours.empty and ganak.empty and approx.empty)
    if not ours.empty:
        plt.plot(
            ours["variable_value"],
            ours["model_count"],
            label="Ours",
            color="#ff7f0e",
            linewidth=3,
            zorder=1,
        )

    if not ganak.empty:
        plt.scatter(
            ganak["variable_value"],
            ganak["model_count"],
            color="#1f77b4",
            marker="s",
            s=200,
            linewidths=1,
            edgecolor="black",
            label="Ganak (Exact)",
            zorder=2,
        )

    if not approx.empty:
        approx_lower = approx["model_count"] * (1 - config.epsilon)
        approx_upper = approx["model_count"] * (1 + config.epsilon)
        plt.errorbar(
            approx["variable_value"],
            approx["model_count"],
            yerr=[approx["model_count"] - approx_lower, approx_upper - approx["model_count"]],
            fmt="o",
            color="#2ca02c",
            ecolor="red",
            elinewidth=3,
            capsize=5,
            markersize=8,
            markeredgewidth=1,
            label="ApproxMC (95% CI)",
            zorder=3,
        )

    if not found_any:
        plt.close()
        print(f"CSV has no known algorithms, skip plotting: {csv_path}")
        return

    y_limits = _paper_y_limits(sweep_name, variable_name)
    is_vary_m = sweep_name == "vary_m" or (sweep_name is None and variable_name == "m_value")

    if is_vary_m:
        # For vary_m, keep paper-style upper bound and force lower bound to 0.
        plt.yscale("symlog", linthresh=1)
        y_label = "Model Count (symlog scale)"
        if y_limits is not None:
            lower, upper = y_limits
            plt.ylim(lower, upper)
            if upper > 0:
                max_power = int(math.log10(upper))
                ticks = []
                if lower < 0:
                    ticks.extend([lower, 0.0])
                elif lower == 0:
                    ticks.append(0.0)
                ticks.extend(10**power for power in range(1, max_power + 1))
                plt.yticks(ticks)
    elif (plot_data["model_count"] > 0).all():
        plt.yscale("log")
        y_label = "Model Count (log scale)"
        if y_limits is not None:
            plt.ylim(*y_limits)
    else:
        plt.yscale("symlog", linthresh=1)
        y_label = "Model Count (symlog scale)"

    ax = plt.gca()
    min_x_value = int(math.floor(float(all_x_values.min())))
    max_x_value = int(math.ceil(float(all_x_values.max())))
    ax.set_xticks(list(range(min_x_value, max_x_value + 1)))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(round(x))}"))

    plt.xlabel(_paper_x_label(sweep_name, variable_name), fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.7)

    file_name = f"{_paper_pdf_basename(sweep_name, variable_name)}.pdf"
    full_path = config.full_fig_results_root / file_name
    full_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(full_path, dpi=600, bbox_inches="tight")

    publish_file_name = resolve_publish_figure_filename("AppendixB.1", file_name)
    if publish_file_name is not None:
        publish_path = config.publish_fig_results_root / publish_file_name
        publish_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(publish_path, dpi=600, bbox_inches="tight")
        print(f"Chart saved to: {publish_path}")
    else:
        print(
            "No paper figure mapping for publish target, skip reproduce/results output: "
            f"AppendixB.1/{file_name}"
        )

    print(f"Chart saved to: {full_path}")
    plt.close()
