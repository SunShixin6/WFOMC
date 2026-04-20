from __future__ import annotations

import math
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from reproduce.utils.figure_naming import resolve_publish_figure_filename

from .benchmarks import RuntimeConfig


def _to_numeric_frame(data: pd.DataFrame) -> pd.DataFrame:
    frame = data.copy()
    frame["domain_size"] = pd.to_numeric(frame["domain_size"], errors="coerce")
    frame["result"] = pd.to_numeric(frame["result"], errors="coerce")
    # Keep only numeric rows here; zero filtering is applied later for plotting.
    frame = frame.dropna(subset=["domain_size", "result"])
    return frame


def plot_model_count(
    csv_file_path: str,
    config: RuntimeConfig,
) -> None:
    """Plot per-model count comparison (Ours vs Ganak vs ApproxMC)."""
    if (not os.path.exists(csv_file_path)) or os.path.getsize(csv_file_path) == 0:
        print(f"CSV is empty, so skip plotting: {csv_file_path}")
        return

    try:
        data = pd.read_csv(csv_file_path)
    except pd.errors.EmptyDataError:
        print(f"CSV has no rows, so skip plotting: {csv_file_path}")
        return

    required_columns = {"formula", "algorithm", "domain_size", "result", "status"}
    if not required_columns.issubset(data.columns):
        missing = required_columns.difference(data.columns)
        print(f"CSV missing required columns {sorted(missing)}, skip plotting: {csv_file_path}")
        return

    data = data[~data["status"].astype(str).str.lower().isin(["timeout", "hard_timeout", "error"])]
    if data.empty:
        print(f"No completed rows after filtering, skip plotting: {csv_file_path}")
        return

    formula = str(data["formula"].iloc[0])
    ours_data = data[
        data["algorithm"].astype(str).str.contains("incremental3|dr", case=False, na=False)
    ]
    ganak_data = data[data["algorithm"].astype(str).str.contains("ganak", case=False, na=False)]
    approx_data = data[data["algorithm"].astype(str).str.contains("approxmc", case=False, na=False)]

    ours = _to_numeric_frame(ours_data)
    ganak = _to_numeric_frame(ganak_data)
    approx = _to_numeric_frame(approx_data)

    # Preserve full domain range for x-axis ticks, including domains with zero counts.
    all_domain_values = pd.concat(
        [ours["domain_size"], ganak["domain_size"], approx["domain_size"]],
        ignore_index=True,
    ).dropna()

    # Do not render model_count == 0 points.
    ours = ours[ours["result"] > 0].sort_values("domain_size")
    ganak = ganak[ganak["result"] > 0].sort_values("domain_size")
    approx = approx[approx["result"] > 0].sort_values("domain_size")

    if ours.empty and ganak.empty and approx.empty:
        print(f"All model counts are zero after filtering, skip plotting: {csv_file_path}")
        return

    approx_lower = approx_upper = None
    if not approx.empty:
        approx_lower = approx["result"] * (1 - config.epsilon)
        approx_upper = approx["result"] * (1 + config.epsilon)

    plt.close("all")
    plt.figure(figsize=(12, 8), dpi=300)

    plt.plot(
        ours["domain_size"],
        ours["result"],
        label="Ours",
        color="#ff7f0e",
        linewidth=3,
        zorder=1,
    )

    if not ganak.empty:
        plt.scatter(
            ganak["domain_size"],
            ganak["result"],
            color="#1f77b4",
            marker="s",
            s=200,
            linewidths=1,
            edgecolor="black",
            label="Ganak (Exact)",
            zorder=2,
        )

    if not approx.empty and approx_lower is not None and approx_upper is not None:
        plt.errorbar(
            approx["domain_size"],
            approx["result"],
            yerr=[approx["result"] - approx_lower, approx_upper - approx["result"]],
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

    all_result_values = pd.concat(
        [ours["result"], ganak["result"], approx["result"]],
        ignore_index=True,
    )
    has_non_positive = (all_result_values <= 0).any()

    if has_non_positive:
        plt.yscale("symlog", linthresh=1)
        y_label = "Model Count (symlog scale)"
    else:
        plt.yscale("log")
        y_label = "Model Count (log scale)"

    plt.xlabel("Domain Size", fontsize=20)
    plt.ylabel(y_label, fontsize=20)

    # Domain sizes are discrete integers; force integer ticks and labels on x-axis.
    ax = plt.gca()
    if not all_domain_values.empty:
        min_domain = int(math.floor(float(all_domain_values.min())))
        max_domain = int(math.ceil(float(all_domain_values.max())))
        ax.set_xticks(list(range(min_domain, max_domain + 1)))
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(round(x))}"))

    plt.legend(fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()

    max_ours = float(ours["result"].max()) if not ours.empty else 0.0
    max_ganak = float(ganak["result"].max()) if not ganak.empty else 0.0
    max_approx = float(approx["result"].max()) if not approx.empty else 0.0
    max_y_val = max(max_ours, max_ganak, max_approx)

    if max_y_val > 0:
        power = math.ceil(math.log10(max_y_val))
        if max_y_val == 10**power:
            power += 1
        upper_limit = 10**power
    else:
        upper_limit = 10

    if has_non_positive:
        plt.ylim(bottom=-0.5, top=upper_limit)
    else:
        positive_values = all_result_values[all_result_values > 0]
        if positive_values.empty:
            lower_limit = 0.1
        else:
            min_positive = float(positive_values.min())
            lower_limit = min(0.1, max(min_positive / 10, 1e-6))
        plt.ylim(bottom=lower_limit, top=upper_limit)

    results_dir = config.appendix_b1_results_path
    results_dir.mkdir(parents=True, exist_ok=True)

    safe_formula = formula.replace("/", "_")
    source_file_name = f"{safe_formula}_correctness.pdf"
    pdf_path = results_dir / source_file_name
    # png_path = results_dir / f"{safe_formula}_correctness.png"

    plt.savefig(pdf_path, format="pdf", dpi=600, bbox_inches="tight")

    publish_file_name = resolve_publish_figure_filename("AppendixB.1", source_file_name)
    if publish_file_name is not None:
        publish_results_dir = config.repo_root / "reproduce" / "results" / "AppendixB.1"
        publish_results_dir.mkdir(parents=True, exist_ok=True)
        publish_pdf_path = publish_results_dir / publish_file_name
        plt.savefig(publish_pdf_path, format="pdf", dpi=600, bbox_inches="tight")
        print(f"The chart has been saved to:\n{publish_pdf_path}\n")
    else:
        print(
            "No paper figure mapping for publish target, skip reproduce/results output: "
            f"AppendixB.1/{source_file_name}"
        )

    # plt.savefig(png_path, format="png", dpi=600, bbox_inches="tight")
    print(f"The chart has been saved to:\n{pdf_path}\n")
    plt.close()
