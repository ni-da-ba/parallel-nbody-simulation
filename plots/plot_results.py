import matplotlib
matplotlib.use("Agg")

import re
import sys
import matplotlib.pyplot as plt


def parse_summary_file(path):
    rows = []
    pattern = re.compile(r"SUMMARY\s+(.*)")

    with open(path, "r") as f:
        for line in f:
            m = pattern.search(line)
            if not m:
                continue

            row = {}
            for token in m.group(1).split():
                key, value = token.split("=", 1)
                row[key] = value
            rows.append(row)

    if not rows:
        raise RuntimeError(f"No SUMMARY lines found in {path}")

    for row in rows:
        row["ranks"] = int(row["ranks"])
        row["N"] = int(row["N"])
        row["steps"] = int(row["steps"])
        row["dt"] = float(row["dt"])
        row["theta"] = float(row["theta"])
        row["elapsed"] = float(row["elapsed"])
        row["avg_step"] = float(row["avg_step"])
        row["mass_change"] = float(row["mass_change"])
        row["px_change"] = float(row["px_change"])
        row["py_change"] = float(row["py_change"])
        row["pz_change"] = float(row["pz_change"])
        row["p_mag_change"] = float(row["p_mag_change"])

    return rows


def filter_rows(rows, **kwargs):
    out = []
    for row in rows:
        ok = True
        for key, value in kwargs.items():
            if row[key] != value:
                ok = False
                break
        if ok:
            out.append(row)
    return out


def sort_by(rows, key):
    return sorted(rows, key=lambda r: r[key])


def plot_size_sweep(rows, out_prefix="size_sweep"):
    plt.figure()

    configs = [
        ("allpairs", 0.5, "All-Pairs"),
        ("barneshut", 0.5, "Barnes-Hut (θ=0.5)"),
        ("barneshut", 0.2, "Barnes-Hut (θ=0.2)")
    ]

    for mode, theta, label in configs:
        sub = filter_rows(rows, mode=mode, theta=theta)
        sub = sort_by(sub, "N")
        x = [r["N"] for r in sub]
        y = [r["elapsed"] for r in sub]
        plt.plot(x, y, marker="o", label=label)

    plt.xlabel("Number of Particles (N)")
    plt.ylabel("Elapsed Time (s)")
    plt.title("Runtime vs Particle Count")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{out_prefix}_runtime.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure()

    for mode, theta, label in configs:
        sub = filter_rows(rows, mode=mode, theta=theta)
        sub = sort_by(sub, "N")
        x = [r["N"] for r in sub]
        y = [r["p_mag_change"] for r in sub]
        plt.plot(x, y, marker="o", label=label)

    plt.xlabel("Number of Particles (N)")
    plt.ylabel("Momentum Drift Magnitude")
    plt.title("Accuracy vs Particle Count")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{out_prefix}_drift.png", dpi=200, bbox_inches="tight")
    plt.close()


def plot_scaling(rows, out_prefix="scaling"):
    plt.figure()

    configs = [
        (5000, "allpairs", 0.5, "All-Pairs, N=5000"),
        (5000, "barneshut", 0.5, "Barnes-Hut, N=5000"),
        (10000, "allpairs", 0.5, "All-Pairs, N=10000"),
        (10000, "barneshut", 0.5, "Barnes-Hut, N=10000")
    ]

    for N, mode, theta, label in configs:
        sub = filter_rows(rows, N=N, mode=mode, theta=theta)
        sub = sort_by(sub, "ranks")
        x = [r["ranks"] for r in sub]
        y = [r["elapsed"] for r in sub]
        plt.plot(x, y, marker="o", label=label)

    plt.xlabel("MPI Ranks")
    plt.ylabel("Elapsed Time (s)")
    plt.title("Runtime vs MPI Ranks")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{out_prefix}_runtime.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure()

    ranks_present = sorted(set(r["ranks"] for r in rows))
    plt.plot(ranks_present, ranks_present, linestyle="--", marker="o", label="Ideal")

    for N, mode, theta, label in configs:
        sub = filter_rows(rows, N=N, mode=mode, theta=theta)
        sub = sort_by(sub, "ranks")

        t1 = None
        for r in sub:
            if r["ranks"] == 1:
                t1 = r["elapsed"]
                break

        if t1 is None:
            continue

        x = [r["ranks"] for r in sub]
        y = [t1 / r["elapsed"] for r in sub]
        plt.plot(x, y, marker="o", label=label)

    plt.xlabel("MPI Ranks")
    plt.ylabel("Speedup")
    plt.title("Speedup vs MPI Ranks")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{out_prefix}_speedup.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure()

    for N, mode, theta, label in configs:
        sub = filter_rows(rows, N=N, mode=mode, theta=theta)
        sub = sort_by(sub, "ranks")

        t1 = None
        for r in sub:
            if r["ranks"] == 1:
                t1 = r["elapsed"]
                break

        if t1 is None:
            continue

        x = [r["ranks"] for r in sub]
        speedup = [t1 / r["elapsed"] for r in sub]
        efficiency = [s / p for s, p in zip(speedup, x)]
        plt.plot(x, efficiency, marker="o", label=label)

    plt.xlabel("MPI Ranks")
    plt.ylabel("Parallel Efficiency")
    plt.title("Parallel Efficiency vs MPI Ranks")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"{out_prefix}_efficiency.png", dpi=200, bbox_inches="tight")
    plt.close()


def main():
    if len(sys.argv) != 3:
        print("Usage:")
        print("  python3 plot_results.py size size_results.txt")
        print("  python3 plot_results.py scale scale_results.txt")
        sys.exit(1)

    mode = sys.argv[1]
    path = sys.argv[2]

    rows = parse_summary_file(path)

    if mode == "size":
        plot_size_sweep(rows)
    elif mode == "scale":
        plot_scaling(rows)
    else:
        raise RuntimeError("First argument must be 'size' or 'scale'.")


if __name__ == "__main__":
    main()
