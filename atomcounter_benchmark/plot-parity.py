"""
Plot a parity line for the data from two CSV files.
Here, they'll be the outputs of the atomcounter vs the atomistic model

The two input files must share the same column headers
(e.g., Perimeter, Interface, Surface, Total).
"""

from sys import exit
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_args():
    parser = ArgumentParser()
    parser.add_argument("csv_a", help="first CSV file. Atomistic output")
    parser.add_argument("csv_b", help="second CSV file. AtomCounter output")
    parser.add_argument(
        "csv_c", help="third CSV file. Input to both AtomCounter and Atomistic"
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show the plots on screen (in addition to saving)",
    )
    parser.add_argument(
        "--output_dir", type=str, help="directory where parity plots should be saved"
    )

    return parser.parse_args()


def main():
    args = read_args()
    df_a = pd.read_csv(args.csv_a)  # atomistic output
    df_b = pd.read_csv(args.csv_b)  # atomcounter output
    df_c = pd.read_csv(args.csv_c)  # input to both

    df_a.columns = df_a.columns.str.strip()
    df_b.columns = df_b.columns.str.strip()
    df_c.columns = df_c.columns.str.strip()

    common_cols = df_a.columns.intersection(df_b.columns)
    if common_cols.empty:
        exit("Error: the two files share no common column headers.")

    out_dir = args.output_dir
    figsize = (6, 6)
    plt.rcParams.update({"font.size": 15,
        "font.family": "arial"}
        )
    vmax, vmin = 100, 0

    for col in common_cols:

        x = df_a[col]
        y = df_b[col]

#        ss_res = ((x - y) ** 2).sum()
#        ss_tot = ((x - x.mean()) ** 2).sum()
#        r2 = 1 - ss_res / ss_tot if ss_tot != 0 else float("nan")
#        print(col, r2)

        diff = np.abs(100 * (x - y) / y)
        diff[~np.isfinite(diff)] = np.nan

        theta = df_c["Theta"].values
        curv_radii = df_c["R (A)"].values

        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        pad = 0.05 * (hi - lo)
        lo, hi = lo - pad, hi + pad

        ##create parity plots
        plt.figure(figsize=figsize)
        plt.scatter(x, y, c="tab:blue")
        plt.plot([lo, hi], [lo, hi], "k--", lw=1)
        plt.xlabel(f"{col} atoms (atomistic)")
        plt.ylabel(f"{col} atoms (atomcounter)")
        plt.title(f"Parity plot of results: {col}")
        plt.xlim(lo, hi)
        plt.ylim(lo, hi)
        if col == "Total":
            ax = plt.gca()
            for tick in ax.get_xticklabels():
                tick.set_rotation(25)
                tick.set_ha("right")

        plt.tight_layout()
        out_path_1 = out_dir + f"/parity_{col.lower()}.png"
        plt.savefig(out_path_1, dpi=150)
        if args.show:
            plt.show()
        print(f"Saved {out_path_1}")

        ##create heatmap of differences
        plt.close()
        plt.figure(figsize=figsize)
        sc = plt.scatter(theta, curv_radii, c=diff, cmap="coolwarm", s=60, vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(sc)
        cbar.set_label("Absolute Percent Difference (%)")
        plt.ylabel(r"$R$ (Å)")
        plt.xlabel("θ (°)")
        plt.title(f"Heatmap of differences in results: {col}")
        plt.tight_layout()
        out_path_2 = out_dir + f"/heatmap_{col.lower()}.png"
        plt.savefig(out_path_2, dpi=150)
        if args.show:
            plt.show()
        plt.close()
        print(f"Saved {out_path_2}")


if __name__ == "__main__":
    main()
