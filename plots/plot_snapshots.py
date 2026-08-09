import matplotlib
matplotlib.use("Agg")

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def plot_snapshot(csv_path):
    df = pd.read_csv(csv_path)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(df["x"], df["y"], df["z"], s=15)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(csv_path.stem)

    fig.savefig(csv_path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_snapshots.py demo_step0000.csv [more files...]")
        sys.exit(1)

    for arg in sys.argv[1:]:
        plot_snapshot(Path(arg))


if __name__ == "__main__":
    main()
