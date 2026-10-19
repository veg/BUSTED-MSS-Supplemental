#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Optional, Tuple

# -------- Helpers --------

def get_gene_id(path: Path) -> str:
    """Gene id = filename up to first '.' (e.g., 'aaeB_13_species')."""
    return path.stem.split('.', 1)[0]

def find_log_likelihood(obj) -> Optional[float]:
    """
    Extract the fitted model log-likelihood from a HyPhy BUSTED JSON.
    Primary path: ["fits"]["Unconstrained model"]["Log Likelihood"]
    Falls back to searching common variants if needed.
    """
    try:
        return float(obj["fits"]["Unconstrained model"]["Log Likelihood"])
    except Exception:
        pass

    # Flexible fallbacks (just in case)
    # 1) Any 'fits' child that has a 'Log Likelihood' number
    try:
        print("Couldn't find log likelihood value for a gene. Searching the file for it (bad idea).")
        fits = obj.get("fits", {})
        for _, fit in fits.items():
            if isinstance(fit, dict):
                for k, v in fit.items():
                    if isinstance(k, str) and "lik" in k.lower():
                        try:
                            return float(v)
                        except Exception:
                            continue
    except Exception:
        pass

    # 2) Last-resort recursive search for a single 'Log Likelihood' float
    def _walk(x):
        print("Really couldn't find log likelihood value for a gene. Searching the WHOLE file for it (really bad idea).")
        if isinstance(x, dict):
            for k, v in x.items():
                if isinstance(k, str) and "lik" in k.lower():
                    try:
                        return float(v)
                    except Exception:
                        pass
                found = _walk(v)
                if found is not None:
                    return found
        elif isinstance(x, list):
            for it in x:
                found = _walk(it)
                if found is not None:
                    return found
        return None

    return _walk(obj)

def collect_lnl_from_dir(d: Path) -> Dict[str, float]:
    """Read *.json files from a directory and map gene_id -> log-likelihood."""
    results = {}
    for p in sorted(d.glob("*.json")):
        try:
            with p.open() as f:
                data = json.load(f)
            lnl = find_log_likelihood(data)
            if lnl is not None and math.isfinite(lnl):
                results[get_gene_id(p)] = lnl
        except Exception:
            # Skip unreadable/bad JSON files gracefully
            continue
    return results

# -------- Main --------

def main():
    ap = argparse.ArgumentParser(
        description="Extract BUSTED log-likelihoods per gene from four model folders and compute differences."
    )
    ap.add_argument("--busted-nosrv-e", required=True, type=Path,
                    help="Folder with BUSTED-noSRV-E JSONs")
    ap.add_argument("--busted-srv-e", required=True, type=Path,
                    help="Folder with BUSTED-SRV-E JSONs")
    ap.add_argument("--busted-mss-2classga-nosrv", required=True, type=Path,
                    help="Folder with BUSTED-MSS-2classGA-noSRV JSONs")
    ap.add_argument("--busted-mss-2classga-srv", required=True, type=Path,
                    help="Folder with BUSTED-MSS-2classGA-SRV JSONs")
    ap.add_argument("-o", "--output", required=True, type=Path,
                    help="Output CSV path")

    args = ap.parse_args()

    # Load each directory
    lnl_nosrv_e   = collect_lnl_from_dir(args.busted_nosrv_e)
    lnl_srv_e     = collect_lnl_from_dir(args.busted_srv_e)
    lnl_mss_nsrv  = collect_lnl_from_dir(args.busted_mss_2classga_nosrv)
    lnl_mss_srv   = collect_lnl_from_dir(args.busted_mss_2classga_srv)

    # Union of gene ids from any model (don’t assume all completed)
    all_genes = sorted(set().union(lnl_nosrv_e, lnl_srv_e, lnl_mss_nsrv, lnl_mss_srv))

    # Prepare CSV
    header = [
        "gene",
        "BUSTED-noSRV-E_lnL",
        "BUSTED-SRV-E_lnL",
        "BUSTED-MSS-2classGA-noSRV_lnL",
        "BUSTED-MSS-2classGA-SRV_lnL",
        # Requested differences:
        "diff__BUSTED-SRV-E_minus_BUSTED-noSRV-E",
        "diff__BUSTED-MSS-2classGA-noSRV_minus_BUSTED-noSRV-E",
        "diff__BUSTED-MSS-2classGA-SRV_minus_BUSTED-noSRV-E",
    ]

    with args.output.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for g in all_genes:
            nosrv = lnl_nosrv_e.get(g)
            srv   = lnl_srv_e.get(g)
            mns   = lnl_mss_nsrv.get(g)
            msr   = lnl_mss_srv.get(g)

            # Compute differences only when both terms exist
            def _diff(a: Optional[float], b: Optional[float]) -> Optional[float]:
                return None if (a is None or b is None) else (a - b)

            row = [
                g,
                nosrv,
                srv,
                mns,
                msr,
                _diff(srv, nosrv),
                _diff(mns, nosrv),
                _diff(msr, nosrv),
            ]
            w.writerow(row)

    print(f"Wrote: {args.output}")

if __name__ == "__main__":
    main()
