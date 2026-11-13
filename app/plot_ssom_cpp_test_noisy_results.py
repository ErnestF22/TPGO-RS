import os
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Note: relative paths are relative to the som project root, not this file's location

LEN_TIMESTAMP = 17  # length of the timestamp suffix in each result folder name
NUM_TESTS_PER_INSTANCE = 10  # must match what was used in the C++ tests
MINDEG_TO_PLOT = 3

# Each sigma appears once (the original list had 0.01, 0.02, 0.05 twice)
XPOINTS = np.array([0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.5, 1, 2])

SAVE_DIR = "/home/rimlab/workspace/matlab_ws/som/plot_cpp/"
SAVE_FORMATS = ("png", "jpg", "pdf", "eps")


def save_figure(fig, basename):
    """Saves the figure as png, pdf and eps in SAVE_DIR."""
    for ext in SAVE_FORMATS:
        fig.savefig(os.path.join(SAVE_DIR, f"{basename}.{ext}"), dpi=150)




def read_every_other_line(path):
    """Values are stored on every 2nd line (the odd lines are labels)."""
    with open(path, "r") as fh:
        return np.array(
            [float(line) for count, line in enumerate(fh, start=1) if count % 2 == 0],
            dtype=np.float64,
        )


def parse_testid(name):
    """Returns (n, mindeg, sigma, file_prefix) from a result folder name."""
    testid = name[: len(name) - LEN_TIMESTAMP]
    n_s, mindeg_s, sigma_s = testid.split("_")[:3]
    n = int(n_s[1:])
    mindeg = int(mindeg_s[6:])
    sigma = float(sigma_s[5:]) / 100  # !! /100
    prefix = testid if testid.endswith("_") else testid + "_"  # temporary fix
    return n, mindeg, sigma, prefix


def parse_bool(token):
    token = token.strip().lower()
    if token in ("true", "t", "yes"):
        return 1
    if token in ("false", "f", "no"):
        return 0
    return int(float(token) != 0)


def load_bool_results(results_path, suffix):
    """
    Loads a per-test boolean file (one value per test, every 2nd line) from every folder.
    Returns a list of (n, mindeg, sigma, array_of_0_1).
    """
    out = []
    for name in sorted(os.listdir(results_path)):
        n, mindeg, sigma, prefix = parse_testid(name)
        path = os.path.join(results_path, name, prefix + suffix)
        with open(path, "r") as fh:
            vals = [parse_bool(line) for count, line in enumerate(fh, start=1) if count % 2 == 0]
        out.append((n, mindeg, sigma, np.array(vals, dtype=int)))
    return out


def load_results(results_path, load_lambdas=False):
    """
    Returns (tuples, tuples_lambdas_acceptable).
    tuples: (n, mindeg, sigma, rot, transl, scale_mean, scale_max, exec_time)
    tuples_lambdas_acceptable: (n, mindeg, sigma, lambdas_acceptable) -- only if load_lambdas
    """
    tuples = []
    tuples_lambdas = []

    for name in sorted(os.listdir(results_path)):
        folder = os.path.join(results_path, name)
        n, mindeg, sigma, prefix = parse_testid(name)

        def mean_of(filename):
            return read_every_other_line(os.path.join(folder, prefix + filename)).mean()

        tuples.append(
            (
                n,
                mindeg,
                sigma,
                mean_of("rot_errors_mean.txt"),
                mean_of("transl_errors_mean.txt"),
                mean_of("lambda_mean_errors_mean.txt"),
                mean_of("lambda_max_errors_mean.txt"),
                mean_of("exec_times.txt"),
            )
        )

        if load_lambdas:
            lam = read_every_other_line(os.path.join(folder, prefix + "lambdas_acceptable.txt"))
            tuples_lambdas.append((n, mindeg, sigma, lam))

    return tuples, tuples_lambdas


def to_series(tuples, mindeg=None):
    """Map result tuples (filtered by mindeg) onto XPOINTS. Missing sigmas stay NaN."""
    if mindeg is None:
        mindeg = MINDEG_TO_PLOT
    keys = ["rot", "transl", "scale_mean", "scale_max", "exec"]
    series = {k: np.full_like(XPOINTS, np.nan, dtype=np.float64) for k in keys}

    for t in tuples:
        if t[1] != mindeg:
            continue
        idx = np.where(np.isclose(XPOINTS, t[2]))[0]  # tolerant float matching
        if idx.size == 0:
            print(f"Warning: sigma={t[2]} not in XPOINTS, skipping")
            continue
        for k, value in zip(keys, t[3:8]):
            series[k][idx] = value
    return series


# ############################ LOAD ############################

tuples_rs, tuples_lambdas_acceptable = load_results("results_lsom/", load_lambdas=True)
tuples_icp, _ = load_results("results_ssom_icp/")
tuples_procrustes, _ = load_results("results_ssom_procrustes/")
tuples_mm, _ = load_results("results_ssom_mm/")
tuples_rs_useful = load_bool_results("results_ssom_mm/", "rs_actually_useful.txt")

# # ############################ ACCEPTABLE LAMBDAS (RS) ############################

# rows = [(t[0], idx) for idx, t in enumerate(tuples_lambdas_acceptable) if t[1] == MINDEG_TO_PLOT]
# rows.sort(key=lambda x: x[0])

# if not rows:
#     print("No acceptable-lambda rows to plot for mindeg", MINDEG_TO_PLOT)
# else:
#     mat = np.zeros((len(rows), NUM_TESTS_PER_INSTANCE), dtype=int)
#     ylabels = []
#     for i, (n_val, orig_idx) in enumerate(rows):
#         t = tuples_lambdas_acceptable[orig_idx]
#         arr = np.asarray(t[3], dtype=int)
#         k = min(arr.size, NUM_TESTS_PER_INSTANCE)  # pad/truncate to expected length
#         mat[i, :k] = arr[:k]
#         ylabels.append(f"n={n_val}, σ={t[2]}")

#     fig_table, ax_table = plt.subplots(
#         figsize=(max(6, NUM_TESTS_PER_INSTANCE), max(3, len(rows) * 0.3))
#     )
#     ax_table.imshow(
#         mat, cmap=ListedColormap(["red", "green"]), aspect="auto",
#         interpolation="nearest", vmin=0, vmax=1,
#     )
#     ax_table.set_yticks(np.arange(len(rows)))
#     ax_table.set_yticklabels(ylabels, fontsize=8)
#     ax_table.set_xticks(np.arange(NUM_TESTS_PER_INSTANCE))
#     ax_table.set_xticklabels([f"run {i + 1}" for i in range(NUM_TESTS_PER_INSTANCE)])
#     ax_table.set_title("Acceptable lambdas (red=0, green=1)")
#     fig_table.tight_layout()

# ############################ RS USEFULNESS (SSOM MM) ############################

MAX_SIGMA_USEFULNESS = 0.2  # exclude the "too noisy" results
tuples_rs_useful = [t for t in tuples_rs_useful if t[2] <= MAX_SIGMA_USEFULNESS + 1e-9]

# Overall percentage of useful tests
all_flags = np.concatenate([t[3] for t in tuples_rs_useful]) if tuples_rs_useful else np.array([])
if all_flags.size:
    print(f"RS usefulness, all mindeg: {100 * all_flags.mean():.1f}% ({all_flags.sum()}/{all_flags.size})")

# Use all N / sigma values so every N appears on the axis (empty if no data for given mindeg)
ns = sorted({t[0] for t in tuples_rs_useful})
sigmas = sorted({t[2] for t in tuples_rs_useful})

for useful_mindeg in (2, 3):
    useful_rows = [t for t in tuples_rs_useful if t[1] == useful_mindeg]
    if not useful_rows:
        print("No RS usefulness rows to plot for mindeg", useful_mindeg)
        continue

    flags = np.concatenate([t[3] for t in useful_rows])
    pct = 100 * flags.mean()
    print(f"RS usefulness, mindeg {useful_mindeg}: {pct:.1f}% ({flags.sum()}/{flags.size})")

    cell = {(t[0], t[2]): t[3] for t in useful_rows}

    fig_useful, ax_u = plt.subplots(
        figsize=(max(6, len(ns) * 1.2), max(4, len(sigmas) * 0.8))
    )
    pad = 0.1  # gap between neighbouring cells
    for xi, n_val in enumerate(ns):
        for yi, s_val in enumerate(sigmas):
            arr = cell.get((n_val, s_val))
            if arr is None:
                continue  # no data for this (N, sigma): leave empty
            k = max(NUM_TESTS_PER_INSTANCE, arr.size)
            bar_w = (1 - 2 * pad) / k
            for j in range(k):
                if j >= arr.size:
                    continue  # missing test result: leave empty
                color = "green" if arr[j] else "red"
                ax_u.add_patch(
                    plt.Rectangle(
                        (xi - 0.5 + pad + j * bar_w, yi - 0.5 + pad),
                        bar_w * 0.85, 1 - 2 * pad,
                        facecolor=color, edgecolor="none",
                    )
                )

    ax_u.set_xlim(-0.5, len(ns) - 0.5)
    ax_u.set_ylim(-0.5, len(sigmas) - 0.5)
    ax_u.set_xticks(np.arange(len(ns)))
    ax_u.set_xticklabels([str(n_val) for n_val in ns])
    ax_u.set_yticks(np.arange(len(sigmas)))
    ax_u.set_yticklabels([str(s_val) for s_val in sigmas])
    ax_u.set_xlabel("N")
    ax_u.set_ylabel("σ (noise)")
    ax_u.set_title(
        f"SSOM MM: RS usefulness, mindeg {useful_mindeg} "
        f"(overall {pct:.1f}%; green=useful, red=not)"
    )
    ax_u.set_xticks(np.arange(len(ns) + 1) - 0.5, minor=True)
    ax_u.set_yticks(np.arange(len(sigmas) + 1) - 0.5, minor=True)
    ax_u.grid(which="minor", color="lightgray", linewidth=0.5)
    ax_u.tick_params(which="minor", length=0)
    fig_useful.tight_layout()
    save_figure(fig_useful, f"rs_usefulness_mindeg{useful_mindeg}")
    plt.close(fig_useful)

# ############################ MAIN PLOT ############################

methods = {
    "TPGO-RS": (to_series(tuples_rs), "green"),
    "TPGO-ICP": (to_series(tuples_icp), "red"),
    "TPGO-PROCR": (to_series(tuples_procrustes), "blue"),
    "SSOM MM": (to_series(tuples_mm), "orange"),
}

panels = [
    ("rot", "Rotation errors", None),
    ("transl", "Translation errors", None),
    ("scale_mean", "Scale errors avg (mean)", 10),
    ("scale_max", "Scale errors avg (max)", 10),
    ("exec", "Execution times [ms]", None),
]

fig = plt.figure()
gs = fig.add_gridspec(len(panels), hspace=1)
axs = gs.subplots(sharex=True)
fig.suptitle("mindeg " + str(MINDEG_TO_PLOT))



# ############################ PER-N ERROR PLOTS (ROT, TRANSL, LAMBDAS) ############################

METHODS_PER_N = {
    "LSOM": tuples_rs,
    "SSOM ICP": tuples_icp,
    "SSOM Procrustes": tuples_procrustes,
    "SSOM MM": tuples_mm,
}
METRICS_PER_N = {
    "rot": "Rotation error",
    "transl": "Translation error",
    "scale_mean": "Lambda mean error",
    "scale_max": "Lambda max error",
}

all_ns = sorted({t[0] for tuples in METHODS_PER_N.values() for t in tuples})
for mindeg_val in (2, 3):
    for n_val in all_ns:
        series_per_method = {
            name: to_series([t for t in tuples if t[0] == n_val], mindeg=mindeg_val)
            for name, tuples in METHODS_PER_N.items()
        }
        for metric, metric_label in METRICS_PER_N.items():
            fig_n, ax_n = plt.subplots(figsize=(7, 4.5))
            for name, series in series_per_method.items():
                ax_n.plot(XPOINTS, series[metric], marker="o", linestyle="none", label=name)
            ax_n.set_xlabel("sigma")
            ax_n.set_ylabel(metric_label)
            ax_n.set_title(f"{metric_label} vs sigma (n={n_val}, mindeg={mindeg_val})")
            ax_n.grid(True, color="0.85")  # solid gray instead of alpha
            ax_n.legend(framealpha=1.0)  # opaque legend background
            fig_n.tight_layout()

            save_figure(fig_n, f"n{n_val}_mindeg{mindeg_val}_{metric}")
            plt.close(fig_n)

for ax, (key, title, clip) in zip(axs, panels):
    for label, (series, color) in methods.items():
        y = series[key]
        if clip is not None:
            y = np.minimum(y, clip)  # NaN stays NaN
        ax.plot(XPOINTS, y, "o", label=label, color=color)
    if clip is not None:
        ax.set_ylim(top=clip)
    ax.set_title(title)

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
    ax.legend(framealpha=1.0)  # opaque legend background

# save_figure(fig, f"main_mindeg{MINDEG_TO_PLOT}")

def close_on_key(event):
    plt.close("all")


for num in plt.get_fignums():
    plt.figure(num).canvas.mpl_connect("key_press_event", close_on_key)

plt.show()

# per_key = Counter()
# for n, mindeg, sigma, arr in tuples_rs_useful:
#     per_key[(mindeg, n, sigma)] += arr.size
#     print(f"mindeg={mindeg} n={n} sigma={sigma} tests={arr.size}")

# for mindeg in (2, 3):
#     rows = [(k, v) for k, v in per_key.items() if k[0] == mindeg]
#     print(f"mindeg {mindeg}: {len(rows)} (n, sigma) keys, {sum(v for _, v in rows)} tests")

# dupes = [k for k, v in Counter((t[1], t[0], t[2]) for t in tuples_rs_useful).items() if v > 1]
# print("Duplicate (mindeg, n, sigma) folders:", dupes)
