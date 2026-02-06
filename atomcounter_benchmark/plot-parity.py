"""
Plot a parity line for the data from two CSV files.
Here, they'll be the outputs of the atomcounter vs the atomistic model

The two input files must share the same column headers
(e.g., Perimeter, Interface, Surface, Total).
"""

import sys
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def read_args():
    parser = ArgumentParser()
    parser.add_argument("csv_a", help="first CSV file. Atomistic output")
    parser.add_argument("csv_b", help="second CSV file. AtomCounter output")
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
    df_a = pd.read_csv(args.csv_a)
    df_b = pd.read_csv(args.csv_b)

    df_a.columns = df_a.columns.str.strip()
    df_b.columns = df_b.columns.str.strip()

    common_cols = df_a.columns.intersection(df_b.columns)
    if common_cols.empty:
        sys.exit("Error: the two files share no common column headers.")

    out_dir = args.output_dir

    ##create parity plots
    for col in common_cols:
        x = df_a[col]
        y = df_b[col]

        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        pad = 0.05 * (hi - lo)
        lo, hi = lo - pad, hi + pad

        plt.figure(figsize=(4, 4))
        plt.scatter(x, y, c="tab:blue")
        plt.plot([lo, hi], [lo, hi], "k--", lw=1)
        plt.xlabel(f"{col} – file A")
        plt.ylabel(f"{col} – file B")
        plt.title(f"Parity plot of results: {col}")
        plt.xlim(lo, hi)
        plt.ylim(lo, hi)
        plt.tight_layout()

        out_path = out_dir + f"/parity_{col.lower()}.png"
        plt.savefig(out_path, dpi=120)
        if args.show:
            plt.show()
        else:
            plt.close()

        print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
