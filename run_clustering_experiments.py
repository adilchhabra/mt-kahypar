import subprocess
import os
import csv
import pandas as pd
import numpy as np
from collections import Counter
from sklearn.metrics import (
    adjusted_mutual_info_score, adjusted_rand_score,
    normalized_mutual_info_score, homogeneity_score,
    completeness_score, v_measure_score
)

# === CONFIGURATION ===
GRAPHS = [
    "primary_school.hgr",
    "high-school.hgr",
    "house-bills.hgr",
    "house-committees.hgr",
    "senate-bills.hgr",
    "trivago-clicks.hgr"#,
    #"walmart-trips.hgr"
]

GROUND_TRUTH_FILES = {
    "primary_school.hgr": "data/node-labels-contact-primary-school.txt",
    "high-school.hgr": "data/node-labels-contact-high-school.txt",
    "house-bills.hgr": "data/node-labels-house-bills.txt",
    "house-committees.hgr": "data/node-labels-house-committees.txt",
    "senate-bills.hgr": "data/node-labels-senate-bills.txt",
    "trivago-clicks.hgr": "data/node-labels-trivago-clicks.txt"#,
    #"walmart-trips.hgr": "data/node-labels-walmart-trips.txt"
}
THETA_VALUES = [round(x * 0.1, 1) for x in range(3, 8)]
SEEDS = [902394, 781321, 529182, 342198, 439820, 123981, 298342, 618234, 849302, 719283]
DATA_FOLDER = "data"
OUTPUT_BASE_FOLDER = "experiments_pimodalt"
EXECUTABLE = "./build/mt-kahypar/application/MtKaHyPar"
NUM_THREADS = 1
VERBOSE = 0

# === METRIC COLUMNS ===
summary_columns = [
    "objective", "modularity", "time",
    "AMI", "ARI", "NMI", "Homogeneity", "Completeness", "V-Measure", "Purity", "F1"
]

# === CLUSTER METRIC HELPERS ===
def read_clusters(file_path):
    with open(file_path, "r") as f:
        return [int(line.strip()) for line in f]

def compute_purity(predicted, truth):
    total = len(truth)
    contingency = Counter(zip(predicted, truth))
    return sum(max(contingency.get((c, g), 0) for g in set(truth)) for c in set(predicted)) / total

def compute_f1(predicted, truth):
    #return 1.0
    n = len(predicted)
    def count_pairs(clusters):
        cnt = Counter(clusters)
        return sum(v * (v - 1) // 2 for v in cnt.values())
    TP = sum(1 for i in range(n) for j in range(i + 1, n)
             if predicted[i] == predicted[j] and truth[i] == truth[j])
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
    parts = line.strip().split(',')
    if len(parts) < 20:
        raise ValueError(f"Unexpected output format: {line}")
    return {
        "graph": parts[2],
        "theta": float(parts[8]),
        "seed": int(parts[4]),
        "objective": float(parts[9]),
        "modularity": float(parts[10]),
        "time": float(parts[14])
    }

# === MAIN LOOP ===
os.makedirs(OUTPUT_BASE_FOLDER, exist_ok=True)
total_runs = len(GRAPHS) * len(THETA_VALUES) * len(SEEDS)
run_counter = 1
global_summary_rows = []

for theta in THETA_VALUES:
    theta_str = f"{theta:.1f}"
    theta_float_str = f"{theta:.6f}"
    folder = os.path.join(OUTPUT_BASE_FOLDER, theta_str)
    os.makedirs(folder, exist_ok=True)
    csv_path = os.path.join(folder, f"results_theta_{theta_str}.csv")

    with open(csv_path, "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["graph", "theta", "seed"] + summary_columns)
        writer.writeheader()

        for graph in GRAPHS:
            if graph not in GROUND_TRUTH_FILES:
                print(f"[WARN] No ground truth for {graph}, skipping.")
                continue

            graph_path = os.path.join(DATA_FOLDER, graph)
            gt_path = GROUND_TRUTH_FILES[graph]
            ground_truth = read_clusters(gt_path)

            for seed in SEEDS:
                print(f"[{run_counter}/{total_runs}] ▶️ Running {graph} θ={theta} seed={seed}")
                cmd = [
                    EXECUTABLE, "-h", graph_path,
                    "--preset-type=cluster", "-o", "pimod",
                    "--write-partition-file=true",
                    "--r-fm-type=do_nothing",
                    "-t", str(NUM_THREADS),
                    "--theta", str(theta),
                    "--v", str(VERBOSE),
                    "--csv", "1",
                    "--partition-output-folder", folder,
                    "--seed", str(seed)
                ]

                try:
                    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    output_lines = result.stdout.strip().splitlines()
                    for line in output_lines:
                        if line.startswith("MT-KaHyPar"):
                            record = parse_output_line(line)
                            output_file = os.path.join(folder, f"{graph}.theta{theta_float_str}.seed{seed}.KaHyPar")
                            if not os.path.exists(output_file):
                                print(f"[WARN] Missing cluster file: {output_file}")
                                break

                            predicted = read_clusters(output_file)
                            if len(predicted) != len(ground_truth):
                                print(f"[WARN] Cluster size mismatch: predicted={len(predicted)}, truth={len(ground_truth)}")
                                break

                            metrics = compute_metrics(predicted, ground_truth)
                            record.update(metrics)
                            writer.writerow(record)
                            print(f"[{run_counter}/{total_runs}] ✅ Done")
                            break
                    else:
                        print(f"[WARN] No valid MT-KaHyPar output for {graph}, seed={seed}")

                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] MtKaHyPar failed: {e}")
                    print(e.stderr)

                run_counter += 1

    # === Summary stats for this theta
    try:
        df = pd.read_csv(csv_path)
        summary_row = {"theta": theta}
        print(f"\n📊 Arithmetic Means for θ = {theta_str}")
        for col in summary_columns:
            if col in df.columns:
                mean_val = df[col].dropna().mean()
                summary_row[col] = mean_val
                print(f"  {col:>15}: {mean_val:.6f}")
            else:
                summary_row[col] = float('nan')
                print(f"  {col:>15}: [missing]")
        print("")
        global_summary_rows.append(summary_row)
    except Exception as e:
        print(f"[ERROR] Failed to compute summary stats for {csv_path}: {e}")

# === GLOBAL SUMMARY CSV ===
summary_path = os.path.join(OUTPUT_BASE_FOLDER, "global_summary.csv")
try:
    with open(summary_path, "w", newline="") as summary_file:
        writer = csv.DictWriter(summary_file, fieldnames=["theta"] + summary_columns)
        writer.writeheader()
        for row in global_summary_rows:
            writer.writerow(row)
    print(f"📄 Global summary written to: {summary_path}")
except Exception as e:
    print(f"[ERROR] Could not write global summary CSV: {e}")

# === GLOBAL SUMMARY BY GRAPH ===
summary_by_graph_path = os.path.join(OUTPUT_BASE_FOLDER, "global_summary_by_graph.csv")
grouped_rows = []

try:
    for theta in THETA_VALUES:
        theta_str = f"{theta:.1f}"
        csv_path = os.path.join(OUTPUT_BASE_FOLDER, theta_str, f"results_theta_{theta_str}.csv")
        df = pd.read_csv(csv_path)
        if df.empty:
            continue

        for graph_name in df["graph"].unique():
            df_graph = df[df["graph"] == graph_name]
            row = {"graph": graph_name, "theta": theta}
            for col in summary_columns:
                if col in df_graph.columns:
                    row[col] = df_graph[col].dropna().mean()
                else:
                    row[col] = float('nan')
            grouped_rows.append(row)

    with open(summary_by_graph_path, "w", newline="") as summary_file:
        writer = csv.DictWriter(summary_file, fieldnames=["graph", "theta"] + summary_columns)
        writer.writeheader()
        for row in grouped_rows:
            writer.writerow(row)

    print(f"📄 Global per-graph summary written to: {summary_by_graph_path}")
except Exception as e:
    print(f"[ERROR] Could not write per-graph global summary CSV: {e}")

print("🎉 All experiments complete.")
