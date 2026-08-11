#!/usr/bin/env python3

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import os
import matplotlib as mpl

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# Load data
df = pd.read_csv("ama_human_mosquito_head_haplotype_occurrences_total_months_active.csv")

# Sort the haplotypes 
haplotype_order = sorted(df['haplotype'], key=lambda x: int(x[1:]))  # Sort H1, H2, H3...
df['haplotype'] = pd.Categorical(df['haplotype'], categories=haplotype_order, ordered=True)
df = df.sort_values('haplotype')

# Normalize haplotype occurrences by months active
df['norm_occurrences_2017_2018'] = df['number_of_occurrences_2017_2018'] / df['total_months_active_2017_2018']
df['norm_occurrences_2020_2021'] = df['number_of_occurrences_2020_2021'] / df['total_months_active_2020_2021']

# Plot the diverging plot
def create_diverging_bar_plot(df_subset, title, ax, xlim=40):
    if df_subset.empty:
        return
    df_plot = df_subset.copy()
    df_plot['norm_occurrences_2017_2018'] = -df_plot['norm_occurrences_2017_2018']  # negative for left bars
    df_plot = df_plot[::-1]  

    y = np.arange(len(df_plot))
    bar_height = 0.20  
    neg_bars = ax.barh(y, df_plot['norm_occurrences_2017_2018'], height=bar_height, color='#1f77b4', label='2017–2018')
    pos_bars = ax.barh(y, df_plot['norm_occurrences_2020_2021'], height=bar_height, color='#ff7f0e', label='2020–2021')

    # Axis labels and title
    ax.set_yticks(y)
    ax.set_yticklabels(df_plot['haplotype'], fontsize=9)
    ax.set_ylabel('ama1 haplotype ID', fontsize=12)
    ax.set_xlabel('Time normalised frequencies', fontsize=12)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_ylim(-0.5, len(df_plot) - 0.5)
    ax.margins(x=0)
    ax.set_xlim(-xlim, xlim)
    xticks = np.arange(-xlim, xlim + 1, 5)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(abs(int(x))) if x != 0 else '0' for x in xticks], fontsize=8)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    handles = [Patch(color='#1f77b4', label='2017_2018'),
               Patch(color='#ff7f0e', label='2020_2021')]
    legend = ax.legend(handles=handles, title='Time period', loc='upper right', frameon=False, fontsize=8)
    plt.setp(legend.get_title(), fontsize=10, fontweight='bold')

# Define the output folder
output_folder = "/figures/"
os.makedirs(output_folder, exist_ok=True)

from PIL import Image

declining_haplotypes = df[
    df['norm_occurrences_2020_2021'] < df['norm_occurrences_2017_2018']
]

fig, ax = plt.subplots(figsize=(7.2, 6.5))  
create_diverging_bar_plot(
    declining_haplotypes,
    "Declining Haplotypes",
    ax,
    xlim=50
)

# Panel tag 
ax.text(
    -0.08, 1.02,          
    "A",                  
    transform=ax.transAxes,
    fontsize=14,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False
)
plt.tight_layout()
temp_png = os.path.join(output_folder, "temp.png")
fig.savefig(temp_png, dpi=600)
plt.close(fig)

img = Image.open(temp_png)
img.save(
    os.path.join(output_folder, "Fig4A.tif"),
    format="TIFF",
    compression="tiff_lzw"
)
os.remove(temp_png)


increasing_haplotypes = df[
    df['norm_occurrences_2020_2021'] > df['norm_occurrences_2017_2018']
]
fig, ax = plt.subplots(figsize=(7.2, 6.5))
create_diverging_bar_plot(
    increasing_haplotypes,
    "Increasing Haplotypes",
    ax,
    xlim=50
)
ax.text(
    -0.08, 1.02,          
    "B",                  
    transform=ax.transAxes,
    fontsize=14,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False
)
plt.tight_layout()
temp_png = os.path.join(output_folder, "temp.png")
fig.savefig(temp_png, dpi=600)
plt.close(fig)
img = Image.open(temp_png)
img.save(
    os.path.join(output_folder, "Fig4B.tif"),
    format="TIFF",
    compression="tiff_lzw"
)
os.remove(temp_png)