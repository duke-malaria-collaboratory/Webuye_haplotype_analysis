#!/usr/bin/env python3

# Load the libraries
import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# Set figure fonts and font sizes
FIG_WIDTH = 7.2         
FIG_HEIGHT = 9.5         
POINT_SIZE = 15

AXIS_FONT = 12
TICK_FONT = 6
TAG_FONT = 11
LEGEND_FONT = 12
TIMEPOINT_FONT = 10

# Path to input and output paths
file_path = "/mnt/c/Users/enyambura/Downloads/ama_haplotype_monthly_frequency_for_human_mosquito_head.csv"
output_dir = "/mnt/c/Users/enyambura/Downloads/haplotype_plots/ama/"
os.makedirs(output_dir, exist_ok=True)

# Load the data
df = pd.read_csv(file_path)

df["month_year"] = pd.to_datetime(df["month_year"], errors="coerce")
df["haplotype_list"] = df["haplotype_list"].astype(str)

# Sort haplotypes numerically
df["hap_num"] = df["haplotype_list"].apply(lambda x: int(re.search(r"\d+", x).group()))
df = df.sort_values(["hap_num", "month_year"]).drop(columns="hap_num")

# Filter the August 2018 to December 2017 missing time gap
df = df[~((df["month_year"] >= "2018-08-01") & (df["month_year"] <= "2019-12-31"))]

pre_gap = pd.date_range("2017-06-01", "2018-07-01", freq="MS")
post_gap = pd.date_range("2020-01-01", "2021-11-01", freq="MS")
valid_months = pre_gap.union(post_gap)

df = df[df["month_year"].isin(valid_months)]

month_labels = [d.strftime("%b-%y") for d in valid_months]
df["month_str"] = pd.Categorical(
    df["month_year"].dt.strftime("%b-%y"),
    categories=month_labels,
    ordered=True
)

# Define the haplotype persistence groups/categories
hap_persist = (
    df.groupby("haplotype_list")["month_str"]
    .nunique()
    .reset_index(name="months_detected")
)
total_months = len(valid_months)
def categorize_persistence(m):
    if m >= int(0.8 * total_months):
        return "Persistent"
    elif m >= 2:
        return "Intermittent"
    else:
        return "Transient"
hap_persist["persistence_category"] = hap_persist["months_detected"].apply(categorize_persistence)
df = df.merge(hap_persist, on="haplotype_list", how="left")

# Define the colors for each persistence group
def persistence_color(cat):
    return {
        "Persistent": "#009E73",
        "Intermittent": "#E69F00",
        "Transient": "#0072B2"
    }.get(cat, "gray")
df["color"] = df["persistence_category"].apply(persistence_color)

# Split the haplotypes into four panels. 
unique_haps = df["haplotype_list"].unique().tolist()
haplotype_groups = [
    unique_haps[:90],
    unique_haps[90:180],
    unique_haps[180:270],
    unique_haps[270:360]
]
# Plot supplementary figure 3
def plot_haplotypes(data, hap_subset, panel_tag, filename):

    sub = data[data["haplotype_list"].isin(hap_subset)]
    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
    sns.scatterplot(
        data=sub,
        x="month_str",
        y="haplotype_list",
        s=POINT_SIZE,
        hue="color",
        palette={
            "#009E73": "#009E73",
            "#E69F00": "#E69F00",
            "#0072B2": "#0072B2"
        },
        alpha=0.8,
        edgecolor=None,
        legend=False,
        ax=ax
    )
    # Label the axis titles
    ax.set_xlabel("Month–Year", fontsize=AXIS_FONT)
    ax.set_ylabel("ama1 haplotype", fontsize=AXIS_FONT) 
    ax.set_xticks(np.arange(len(month_labels)))
    ax.set_xticklabels(month_labels, rotation=90, fontsize=TICK_FONT)
    ax.tick_params(axis="y", labelsize=TICK_FONT)
    ax.tick_params(axis="both", which="major", length=4, width=0.4)
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    ax.axvline(month_labels.index("Jul-18") + 0.5, color="gray", linestyle="--", linewidth=1)

    # Add the time period labels on top of each panel
    ax.text(month_labels.index("Oct-17"), -1.5, "2017–2018",
            fontsize=TIMEPOINT_FONT, ha="center")
    ax.text(month_labels.index("Jun-20"), -1.5, "2020–2021",
            fontsize=TIMEPOINT_FONT, ha="center")

    # Label the panels starting with A
    ax.text(
        0.01, 0.99,
        panel_tag,
        transform=ax.transAxes,
        fontsize=TAG_FONT,
        fontweight="bold",
        ha="left",
        va="top"
    )
    # Legend 
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#009E73", markersize=6,
                   label="Persistent (≥30 months)"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#E69F00", markersize=6,
                   label="Intermittent"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#0072B2", markersize=6,
                   label="Transient (1 month)")
    ]
    ax.legend(
        handles=handles,
        title="Persistence category",
        fontsize=LEGEND_FONT,
        title_fontsize=LEGEND_FONT,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
        frameon=False
    )
    # Save the supplementary figure as tif
    temp_png = os.path.join(output_dir, "temp.png")
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(temp_png, dpi=600, bbox_inches="tight")
    plt.close()
    img = Image.open(temp_png)
    tif_path = os.path.join(output_dir, filename)
    img.save(tif_path, format="TIFF", compression="tiff_lzw")
    os.remove(temp_png)

# Generate supplementary figure 3
plot_haplotypes(df, haplotype_groups[0], "A", "FigS3A.tif")
plot_haplotypes(df, haplotype_groups[1], "B", "FigS3B.tif")
plot_haplotypes(df, haplotype_groups[2], "C", "FigS3C.tif")
plot_haplotypes(df, haplotype_groups[3], "D", "FigS3D.tif")

