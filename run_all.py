

import subprocess
import time
import os
import re
import csv

reference_path = "data/test_input/salmonella_sub.fa"
reads_path = "data/test_input/reads_sal_sub.fq"
naive_map_path = "data/test_output/sal_sub_reads.naive.map"
output_csv = "data/results/results.csv"

os.makedirs("data", exist_ok=True)
os.makedirs("data/output", exist_ok=True)
os.makedirs("data/results", exist_ok=True)

header = ["w\\k"] + [str(k) for k in range(3, 11)]
results = []

for w in range(4, 12):
    row = [str(w)]
    for k in range(3, w):
        build_output = f"data/output/build_k{k}_w{w}.bin"
        query_output = f"data/output/query_k{k}_w{w}.map"

        # Run build command and time it
        build_cmd = [
            "cargo", "run", "--bin", "buildsa", reference_path,
            str(k), str(w), build_output
        ]
        start_build = time.time()
        build_proc = subprocess.run(build_cmd, capture_output=True, text=True)
        build_time = time.time() - start_build

        # Parse stdout to get reduction info
        stdout = build_proc.stdout
        orig_len = int(re.search(r"Original sequence length: (\d+)", stdout).group(1))
        mini_len = int(re.search(r"Minimizer sequence length \(indices stored\): (\d+)", stdout).group(1))
        reduction = mini_len / orig_len

        # Get index file size
        index_size = os.path.getsize(build_output)

        # Run query command and time it
        query_cmd = [
            "cargo", "run", "--bin", "querysa", build_output,
            reads_path, query_output
        ]
        start_query = time.time()
        subprocess.run(query_cmd, check=True)
        query_time = time.time() - start_query

        # Run parse_false_neg.py
        neg_proc = subprocess.run(
            ["python3", "parse_false_neg.py", naive_map_path, query_output],
            capture_output=True, text=True, check=True
        )
        fn = float(re.search(r"Average # of false negatives per query: ([0-9.]+)", neg_proc.stdout).group(1))

        # Run parse_false_pos.py
        pos_proc = subprocess.run(
            ["python3", "parse_false_pos.py", naive_map_path, query_output],
            capture_output=True, text=True, check=True
        )
        fp = float(re.search(r"Average # of false positives per query: ([0-9.]+)", pos_proc.stdout).group(1))
        fpr = float(re.search(r"Average # of queries with false positive: ([0-9.]+)", pos_proc.stdout).group(1))

        # Fill result cell
        cell = f"R={reduction:.3f}, BT={build_time:.2f}s, QT={query_time:.2f}s, FS={index_size//1024}KB, FN={fn:.3f}, FP={fp:.3f}, FPR={fpr:.3f}"
        row.append(cell)
    results.append(row)

# Write CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(results)


# Define a helper function to extract specific metrics from a cell string
def extract_metric(cell, label):
    match = re.search(fr"{label}=([0-9.]+)", cell)
    return match.group(1) if match else ""

# Metric labels to extract and output filenames
metrics = {
    "R": "data/results/reduction.csv",
    "BT": "data/results/build_time.csv",
    "QT": "data/results/query_time.csv",
    "FS": "data/results/file_size.csv",
    "FN": "data/results/false_negatives.csv",
    "FP": "data/results/false_positives.csv",
    "FPR": "data/results/false_positive_rate.csv"
}

for label, filename in metrics.items():
    metric_rows = []
    for row in results:
        metric_row = [row[0]]  # first cell is w
        for cell in row[1:]:
            metric_row.append(extract_metric(cell, label))
        metric_rows.append(metric_row)

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(metric_rows)

print("Separate metric CSVs written.")