#!/usr/bin/env python3
"""
Cluster every hyper‑graph inside a folder in one go.

Usage
-----
Default behaviour (identical to the previous script):
    python cluster_hypergraphs_folder.py

Batch‑cluster every graph in a folder that follows the naming scheme
```
    rep<i>_he.hgr      # i = 1 … 25 (hyper‑graph file)
    rep<i>_assign.txt  #            (ground‑truth partition)
```
    python cluster_hypergraphs_folder.py -f /path/to/DCHSBM/scenA1

All result files and CSVs are written to
    experiments/DCHSBM/scenA1/<theta>/ …

Only the segment that changes in the input path (e.g. DCHSBM/scenA1)
is reproduced below the global ``experiments`` directory.
"""

import argparse
import csv
import glob
import os
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    normalized_mutual_info_score,
    v_measure_score,
)

# ------------- CONSTANT DEFAULT CONFIGURATION -------------
# These remain the fall‑back values when *no* folder is given -----------------
DEFAULT_GRAPHS = [
    "primary_school.hgr",
    "high-school.hgr",
    "house-bills.hgr",
    "house-committees.hgr",
    "senate-bills.hgr",
    "trivago-clicks.hgr",
]

DEFAULT_GROUND_TRUTH_FILES = {
    "primary_school.hgr": "data/node-labels-contact-primary-school.txt",
    "high-school.hgr": "data/node-labels-contact-high-school.txt",
    "house-bills.hgr": "data/node-labels-house-bills.txt",
    "house-committees.hgr": "data/node-labels-house-committees.txt",
    "senate-bills.hgr": "data/node-labels-senate-bills.txt",
    "trivago-clicks.hgr": "data/node-labels-trivago-clicks.txt",
}

THETA_VALUES = [round(x * 0.1, 1) for x in range(3, 11)]  # 0.3 … 0.7
#SEEDS = [902394, 781321, 529182, 342198, 439820, 123981, 298342, 618234, 849302, 719283]
SEEDS = [42]
EXECUTABLE = "./build/mt-kahypar/application/MtKaHyPar"
NUM_THREADS = 1
VERBOSE = 0

# ------------- COMMAND‑LINE PARSING / DYNAMIC CONFIGURATION -------------

def discover_graphs(folder: Path):
    """Return a list of graph file names and the corresponding ground‑truth map."""
    pattern = str(folder / "rep*_he.hgr")
    graph_paths = sorted(Path(p) for p in glob.glob(pattern))
    if not graph_paths:
        raise FileNotFoundError(f"No *_he.hgr files found in {folder}")

    graphs = [p.name for p in graph_paths]
    gt_map = {}
    missing_gt = []
    for gp in graph_paths:
        gt_file = gp.with_name(gp.stem.replace("_he", "_assign") + ".txt")
        if not gt_file.exists():
            missing_gt.append(gt_file.name)
        gt_map[gp.name] = str(gt_file)

    if missing_gt:
        print("[WARN] Missing ground‑truth files:")
        for f in missing_gt:
            print(f"   → {f}")
    return graphs, gt_map


def parse_cli() -> tuple[list[str], dict[str, str], str]:
    """Handle --folder/-f argument and return (GRAPHS, GROUND_TRUTH_FILES, DATA_FOLDER)."""
    parser = argparse.ArgumentParser(
        description="Cluster every hyper‑graph in a folder that follows the rep<i> naming scheme.",
    )
    parser.add_argument(
        "-f",
        "--folder",
        metavar="DIR",
        type=str,
        help="Directory that contains rep*_he.hgr and rep*_assign.txt files.",
    )
    args = parser.parse_args()

    if args.folder:
        folder = Path(args.folder).expanduser().resolve()
        graphs, gt_map = discover_graphs(folder)
        data_folder = str(folder)  # hypergraphs live here
    else:
        # fall‑back to the hard‑coded demo graphs
        graphs = list(DEFAULT_GRAPHS)
        gt_map = dict(DEFAULT_GROUND_TRUTH_FILES)
        data_folder = "data"

    return graphs, gt_map, data_folder, args.folder


GRAPHS, GROUND_TRUTH_FILES, DATA_FOLDER, _cli_folder = parse_cli()

# Where to write results ------------------------------------------------------
if _cli_folder:
    # Keep only the two trailing components (e.g. DCHSBM/scenA1)
    cli_path = Path(_cli_folder).resolve()
    tail = Path(*cli_path.parts[-2:]) if len(cli_path.parts) >= 2 else cli_path.name
    OUTPUT_BASE_FOLDER = Path("experiments_strLoyaltyWeighted") / tail
else:
    OUTPUT_BASE_FOLDER = Path("experiments_pimodalt")

OUTPUT_BASE_FOLDER.mkdir(parents=True, exist_ok=True)

# ============= METRIC HELPERS =============

def read_clusters(file_path):
    with open(file_path, "r") as f:
        return [int(line.strip()) for line in f]


def compute_purity(predicted, truth):
    total = len(truth)
    contingency = Counter(zip(predicted, truth))
    return sum(max(contingency.get((c, g), 0) for g in set(truth)) for c in set(predicted)) / total


def count_pairs(clusters):
    cnt = Counter(clusters)
    return sum(v * (v - 1) // 2 for v in cnt.values())


def compute_f1(predicted, truth):
    n = len(predicted)
    TP = sum(
        1
        for i in range(n)
        for j in range(i + 1, n)
        if predicted[i] == predicted[j] and truth[i] == truth[j]
    )
    P = count_pairs(predicted)
    T = count_pairs(truth)
    precision = TP / P if P else 0
    recall = TP / T if T else 0
    return 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0


def compute_metrics(predicted, truth):
    return {
        "AMI": adjusted_mutual_info_score(truth, predicted),
        "ARI": adjusted_rand_score(truth, predicted),
        "NMI": normalized_mutual_info_score(truth, predicted),
        "Homogeneity": homogeneity_score(truth, predicted),
        "Completeness": completeness_score(truth, predicted),
        "V-Measure": v_measure_score(truth, predicted),
        "Purity": compute_purity(predicted, truth),
        "F1": compute_f1(predicted, truth),
    }


def parse_output_line(line):
    parts = line.strip().split(",")
    if len(parts) < 20:
        raise ValueError(f"Unexpected MT‑KaHyPar CSV line: {line}")
    return {
        "graph": parts[2],
        "theta": float(parts[8]),
        "seed": int(parts[4]),
        "objective": float(parts[9]),
        "modularity": float(parts[10]),
        "time": float(parts[14]),
    }

# CSV column order ------------------------------------------------------------
summary_columns = [
    "objective",
    "modularity",
    "time",
    "AMI",
    "ARI",
    "NMI",
    "Homogeneity",
    "Completeness",
    "V-Measure",
    "Purity",
    "F1",
]

# -------------------- MAIN LOOP --------------------

total_runs = len(GRAPHS) * len(THETA_VALUES) * len(SEEDS)
run_counter = 1

global_summary_rows = []

for theta in THETA_VALUES:
    theta_str = f"{theta:.1f}"
    theta_float_str = f"{theta:.6f}"

    theta_folder = OUTPUT_BASE_FOLDER / theta_str
    theta_folder.mkdir(parents=True, exist_ok=True)

    csv_path = theta_folder / f"results_theta_{theta_str}.csv"

    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["graph", "theta", "seed"] + summary_columns)
        writer.writeheader()

        for graph in GRAPHS:
            if graph not in GROUND_TRUTH_FILES:
                print(f"[WARN] No ground‑truth for {graph}, skipping.")
                continue

            graph_path = Path(DATA_FOLDER) / graph
            gt_path = Path(GROUND_TRUTH_FILES[graph])

            if not graph_path.exists():
                print(f"[WARN] Graph file does not exist: {graph_path}")
                continue
            if not gt_path.exists():
                print(f"[WARN] Ground‑truth file does not exist: {gt_path}")
                continue

            ground_truth = read_clusters(gt_path)

            for seed in SEEDS:
                print(f"[{run_counter}/{total_runs}] ▶️ {graph} | θ={theta_str} | seed={seed}")
                cmd = [
                    EXECUTABLE,
                    "-h",
                    str(graph_path),
                    "--preset-type=cluster",
                    "-o",
                    "pimod",
                    "--write-partition-file=true",
                    "--r-fm-type=do_nothing",
                    "-t",
                    str(NUM_THREADS),
                    "--theta",
                    str(theta),
                    "--v",
                    str(VERBOSE),
                    "--csv",
                    "1",
                    "--partition-output-folder",
                    str(theta_folder),
                    "--seed",
                    str(seed),
                ]

                try:
                    result = subprocess.run(
                        cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                    )
                    output_lines = result.stdout.strip().splitlines()
                    for line in output_lines:
                        if line.startswith("MT-KaHyPar"):
                            record = parse_output_line(line)
                            # Mt‑KaHyPar writes the partition file next to --partition-output-folder
                            out_file = theta_folder / f"{graph}.theta{theta_float_str}.seed{seed}.KaHyPar"
                            if not out_file.exists():
                                print(f"[WARN] Missing cluster file: {out_file}")
                                break

                            predicted = read_clusters(out_file)
                            if len(predicted) != len(ground_truth):
                                print(
                                    f"[WARN] Cluster size mismatch ({len(predicted)} vs {len(ground_truth)}) in {out_file}"
                                )
                                break

                            metrics = compute_metrics(predicted, ground_truth)
                            record.update(metrics)
                            writer.writerow(record)
                            print(f"[{run_counter}/{total_runs}] ✅ done")
                            break
                    else:
                        print(f"[WARN] No MT‑KaHyPar CSV line for {graph}, seed {seed}")
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] Mt‑KaHyPar failed for {graph}, seed {seed}: {e}")
                    print(e.stderr)

                run_counter += 1

    # -------- summary for this θ --------
    try:
        df = pd.read_csv(csv_path)
        summary_row = {"theta": theta}
        print(f"\n📊 Arithmetic means for θ = {theta_str}")
        for col in summary_columns:
            mean_val = df[col].dropna().mean() if col in df.columns else np.nan
            summary_row[col] = mean_val
            print(f"  {col:>15}: {mean_val:.6f}" if not np.isnan(mean_val) else f"  {col:>15}: [nan]")
        print("")
        global_summary_rows.append(summary_row)
    except Exception as e:
        print(f"[ERROR] Could not compute summary stats for {csv_path}: {e}")

# -------------------- GLOBAL SUMMARY CSVs --------------------
summary_path = OUTPUT_BASE_FOLDER / "global_summary.csv"
try:
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["theta"] + summary_columns)
        writer.writeheader()
        writer.writerows(global_summary_rows)
    print(f"📄 Global summary written to {summary_path}")
except Exception as e:
    print(f"[ERROR] Could not write {summary_path}: {e}")

# Per‑graph summary -----------------------------------------------------------
summary_by_graph_path = OUTPUT_BASE_FOLDER / "global_summary_by_graph.csv"
per_graph_rows = []
try:
    for theta in THETA_VALUES:
        theta_str = f"{theta:.1f}"
        csv_path = OUTPUT_BASE_FOLDER / theta_str / f"results_theta_{theta_str}.csv"
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        if df.empty:
            continue
        for graph_name in df["graph"].unique():
            df_g = df[df["graph"] == graph_name]
            row = {"graph": graph_name, "theta": theta}
            for col in summary_columns:
                row[col] = df_g[col].dropna().mean() if col in df_g.columns else np.nan
            per_graph_rows.append(row)
    with open(summary_by_graph_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["graph", "theta"] + summary_columns)
        writer.writeheader()
        writer.writerows(per_graph_rows)
    print(f"📄 Global per‑graph summary written to {summary_by_graph_path}")
except Exception as e:
    print(f"[ERROR] Could not write per‑graph summary: {e}")

print("🎉 All experiments complete.")
