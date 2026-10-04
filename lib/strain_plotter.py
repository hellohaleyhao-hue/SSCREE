"""
Strain Plotter Module
---------------------

Automates batch steric strain calculations across protein variants using 
the SSCREE computational engine and renders comparative categorical 
beeswarm plots. Supports single/multi-chain inputs across PDB, MMCIF, 
and PyMOL session files.
"""

import os
import time
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
import lib.sscree as sc

CATEGORY_MAP = {
    'total': (2, 'Total Strain', 'Total Steric Overlap Score'),
    'protein': (3, 'Protein Strain', 'Protein Steric Overlap Score'),
    'nucleic': (4, 'Nucleic Strain', 'Nucleic Steric Overlap Score'),
    'organic': (5, 'Organic Strain', 'Organic Steric Overlap Score'),
    'inorganic': (6, 'Inorganic Strain', 'Inorganic Steric Overlap Score'),
    'solvent': (7, 'Solvent Strain', 'Solvent Steric Overlap Score'),
    'other': (8, 'Other Strain', 'Other Steric Overlap Score'),
}


def parse_variant_line(line: str, default_chain: str = "chain A") -> tuple[str, str]:
    """
    Parses a line into (variant, chain_selection).
    
    Supports formats:
      - "Thr5Ile     A"        -> ("Thr5Ile", "chain A")
      - "Thr5Ile, chain A"    -> ("Thr5Ile", "chain A")
      - "Thr5Ile, A"          -> ("Thr5Ile", "chain A")
      - "Thr5Ile"             -> ("Thr5Ile", default_chain)
    """
    line = line.strip()
    if not line:
        return "", ""

    if "," in line:
        parts = [p.strip() for p in line.split(",", 1)]
        variant, raw_chain = parts[0], parts[1]
    else:
        parts = line.split()
        variant = parts[0]
        raw_chain = parts[1] if len(parts) > 1 else default_chain

    # Normalize chain string for PyMOL selection syntax
    if not raw_chain.lower().startswith("chain"):
        chain = f"chain {raw_chain}"
    else:
        chain = raw_chain

    return variant, chain


def get_strain_data(variant_entries, col_idx, category_title, default_chain="chain A"):
    strain_data = {}
    print()

    for entry in variant_entries:
        variant, chain = parse_variant_line(entry, default_chain=default_chain)
        if not variant:
            continue

        # Execute calculation passing the specified chain
        results = sc.model.calculate_strain_value(variant, chain=chain)
        strain_data[variant] = [[row[0], row[1], row[col_idx]] for row in results]

    print(f"\n{category_title} data for all variants:")
    print("=" * 45)

    for variant, results in strain_data.items():
        print(f"\nVariant: {variant}")
        print(f"{'Rotamer':<8} | {'Frequency':<10} | {category_title:<15}")
        print("-" * 45)
        for row in results:
            print(f"{row[0]:<8} | {row[1]:<10} | {row[2]:<15.4f}")

    print("\n" + "=" * 45)
    print(f"All {category_title.lower()} data collected successfully.\n")

    return strain_data


def plot_strain_data(strain_data_1, strain_data_2, output_path, y_axis_label):
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 9), sharey=True)

    datasets = [
        (strain_data_1, axes[0], True),
        (strain_data_2, axes[1], False)
    ]

    for strain_data, ax, show_legend in datasets:
        rows = []

        for variant, rotamers in strain_data.items():
            for rotamer, frequency_str, strain_value in rotamers:
                rows.append([variant, rotamer, frequency_str, strain_value])

        df = pd.DataFrame(
            rows,
            columns=["Variant", "Rotamer", "Frequency_Str", "Selected_Strain"]
        )

        df["Is_NA"] = df["Frequency_Str"].astype(str).str.upper().str.contains("N/A")
        df["Frequency"] = (
            pd.to_numeric(df["Frequency_Str"].str.rstrip("%"), errors="coerce").fillna(0.0) / 100
        )

        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{x:g}'))
        ax.set_yscale("symlog", linthresh=0.001, linscale=0.3, subs=[2, 3, 4, 5, 6, 7, 8, 9])
        ax.set_yticks([0, 0.001, 0.01, 0.1, 1, 10, 100])
        ax.set_ylim(-0.0008, 150)

        sns.swarmplot(
            data=df,
            x="Variant",
            y="Selected_Strain",
            hue="Is_NA",
            palette={False: sns.color_palette()[0], True: "red"},
            size=6,
            edgecolor="black",
            linewidth=0.1,
            warn_thresh=0.01,
            ax=ax
        )

        if ax.get_legend() is not None:
            ax.get_legend().remove()

        min_size, max_size, na_fixed_size = 20, 500, 80
        df["MarkerSize"] = min_size + df["Frequency"] * (max_size - min_size)
        df.loc[df["Is_NA"], "MarkerSize"] = na_fixed_size

        variants = list(strain_data.keys())
        for collection, variant in zip(ax.collections, variants):
            variant_sizes = df.loc[df["Variant"] == variant, "MarkerSize"].to_numpy()
            collection.set_sizes(variant_sizes)

        # Force a canvas draw to freeze swarm point positions in place
        fig.canvas.draw()

        # Explicitly lock axis limits so adding lines won't trigger re-scaling
        curr_xlim = ax.get_xlim()
        curr_ylim = ax.get_ylim()
        ax.set_xlim(curr_xlim)
        ax.set_ylim(curr_ylim)

        for idx, variant in enumerate(variants):
            sub_df = df[(df["Variant"] == variant) & (~df["Is_NA"])]
            
            if not sub_df.empty and sub_df["Frequency"].sum() > 0:
                weighted_mean = (
                    sub_df["Selected_Strain"] * sub_df["Frequency"]
                ).sum() / sub_df["Frequency"].sum()

                ax.hlines(
                    y=weighted_mean,
                    xmin=idx - 0.2,
                    xmax=idx + 0.2,
                    colors="black",
                    linewidth=1.5,
                    alpha=1,
                    zorder=5
                )

        ax.set_xlabel(None)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
        ax.set_ylabel(f"{y_axis_label} ($\mathrm{{\AA}}^2$)")

        ax.grid(axis='y', linestyle='--', alpha=0.3, color='gray')
        ax.set_axisbelow(True)

        if show_legend:
            legend_frequencies = [0.05, 0.20, 0.40, 0.60]
            legend_handles = []

            for frequency in legend_frequencies:
                marker_size = min_size + frequency * (max_size - min_size)
                legend_handles.append(
                    Line2D(
                        [0], [0],
                        marker="o", linestyle="None",
                        markerfacecolor=sns.color_palette()[0],
                        markersize=marker_size ** 0.5,
                        label=f"{frequency:.0%}"
                    )
                )

            legend_handles.append(
                Line2D(
                    [0], [0],
                    marker="o", linestyle="None",
                    markerfacecolor="red", markeredgecolor="none",
                    markersize=na_fixed_size ** 0.5,
                    label="N/A"
                )
            )

            ax.legend(
                handles=legend_handles,
                title="Frequency",
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                borderaxespad=0,
                frameon=True,
                facecolor='white',
                framealpha=0.85
            )

    fig.subplots_adjust(left=0.08, right=0.88, top=0.93, bottom=0.08, hspace=0.25)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def run_pipeline(
    variant_file_path: str,
    structure_path: str,
    output_path: str,
    strain_category: str = 'protein',
    default_chain: str = 'chain A'
):
    """
    Executes the strain plotting pipeline.

    Args:
        strain_category: Category to plot for the specified structure. Allowed options:
            - `"protein"` (Default)
            - `"nucleic"`
            - `"organic"`
            - `"inorganic"`
            - `"solvent"`
            - `"other"`
            - `"total"`
    """
    sc.model.config(structure_path)

    category_key = strain_category.lower()
    if category_key not in CATEGORY_MAP:
        raise ValueError(f"Invalid category '{strain_category}'. Options: {list(CATEGORY_MAP.keys())}")

    col_idx, category_title, y_axis_label = CATEGORY_MAP[category_key]

    start_time = time.perf_counter()

    with open(variant_file_path, "r") as file:
        variant_list = [line.strip() for line in file if line.strip()]
        variants_1to15 = variant_list[:15]
        variants_16to30 = variant_list[15:30]

    data_1 = get_strain_data(variants_1to15, col_idx, category_title, default_chain=default_chain)
    data_2 = get_strain_data(variants_16to30, col_idx, category_title, default_chain=default_chain)

    end_time = time.perf_counter()
    print("\n" + "=" * 40)
    print(f"Execution complete. Total time taken: {end_time - start_time:.2f} seconds\n")

    all_zero = all(
        row[2] == 0 
        for dataset in (data_1, data_2) 
        for results in dataset.values() 
        for row in results
    )

    if all_zero:
        print("=" * 80 + "\n")
        print(f"NOTICE: ALL CALCULATED {category_title.upper()} VALUES ARE ZERO ACROSS ALL VARIANTS AND ROTAMERS.")
        print("Skipping plot generation.\n")
        print("=" * 80 + "\n")
    else:
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
        plot_strain_data(data_1, data_2, output_path, y_axis_label)
