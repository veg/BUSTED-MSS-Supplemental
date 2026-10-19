import argparse
import json
import csv
import traceback
from pathlib import Path

def extract_fields(json_path):
    try:
        with open(json_path) as f:
            data = json.load(f)

        fit = data["fits"]["Unconstrained model"]
        aic = fit["AIC-c"]
        lrt = data["test results"]["LRT"]
        pval = data["test results"]["p-value"]
        log_likelihood = fit["Log Likelihood"]

        # dS distributions: 3 expected classes (0, 1, 2)
        ds_raw = fit["Rate Distributions"].get("Synonymous site-to-site rates", [])
        ds = {}
        for i in range(3):
            class_data = ds_raw[str(i)] if i < len(ds_raw) else {}
            ds[f"dS_class_{i}_proportion"] = class_data.get("proportion", None)
            ds[f"dS_class_{i}_rate"] = class_data.get("rate", None)

        # omega distributions: ignore class 0 (error sink), map class 1–3 to dS class 0–2
        omega_raw = fit["Rate Distributions"].get("Test", [])
        omega = {}
        for i in range(1, 4):
            class_data = omega_raw[str(i)] if i < len(omega_raw) else {}
            omega[f"omega_class_{i-1}_omega"] = class_data.get("omega", None)
            omega[f"omega_class_{i-1}_proportion"] = class_data.get("proportion", None)

        return {
            "AIC": aic,
            "LogLikelihood": log_likelihood,
            "LRT": lrt,
            "pval": pval,
            **ds,
            **omega
        }

    except Exception as e:
        print(f"Error parsing {json_path}: {e}")
        traceback.print_exc()
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--busted", required=True, help="Path to busted folder")
    parser.add_argument("--busted-mss", required=True, help="Path to busted-mss folder")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    args = parser.parse_args()

    busted_dir = Path(args.busted)
    busted_mss_dir = Path(args.busted_mss)

    # Use only .json files that match the pattern [GENE].[FOLDER].json
    busted_files = {
        f.stem.rsplit(".", 1)[0]: f for f in busted_dir.glob("*.json") if "." in f.stem
    }
    busted_mss_files = {
        f.stem.rsplit(".", 1)[0]: f for f in busted_mss_dir.glob("*.json") if "." in f.stem
    }

    all_genes = sorted(set(busted_files.keys()) & set(busted_mss_files.keys()))

    results = []

    for gene in all_genes:
        busted_data = extract_fields(busted_files[gene])
        busted_mss_data = extract_fields(busted_mss_files[gene])

        if busted_data is None or busted_mss_data is None:
            continue

        row = {"gene": gene}
        for k, v in busted_data.items():
            row[f"busted_{k}"] = v
        for k, v in busted_mss_data.items():
            row[f"busted_mss_{k}"] = v

        results.append(row)

    if results:
        fieldnames = ["gene"] + sorted(k for k in results[0] if k != "gene")
        with open(args.output, "w", newline="") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
    else:
        print("No valid data found to write.")

if __name__ == "__main__":
    main()

