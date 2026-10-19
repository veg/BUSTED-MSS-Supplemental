import csv
import os
import statistics

# Configuration
data_dir = "data"
files = {
    "Caenorhabditis": "20250819_BUSTED_vs_BUSTED-MSS_codonGA_caenor.csv",
    "Drosophila": "20250415_BUSTED_vs_BUSTED-MSS_codonGA_droso.csv",
    "Enterobacteria": "20250819_BUSTED_vs_BUSTED-MSS_codonGA_entero.csv",
    "Primate": "20250820_BUSTED_vs_BUSTED-MSS_codonGA_primate.csv",
    "Saccharomyces": "20250819_BUSTED_vs_BUSTED-MSS_codonGA_sacch.csv"
}

def get_stats(vals):
    if not vals:
        return 0, 0, 0, 0
    return len(vals), statistics.mean(vals), min(vals), max(vals)

print(f"{'Dataset':<15} | {'Group':<12} | {'N':<5} | {'w3':<15} | {'p3':<15} | {'dAICc':<15}")
print("-" * 90)

for clade, filename in files.items():
    filepath = os.path.join(data_dir, filename)
    if not os.path.exists(filepath):
        continue
    
    # Storage
    busted_only_w3 = []
    busted_only_p3 = []
    busted_only_aic = []
    
    mss_only_w3 = []
    mss_only_p3 = []
    mss_only_aic = []
    
    both_w3 = []
    both_p3 = []
    both_aic = []
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # P-values
                p_busted = float(row.get('busted_pval', 1.0))
                p_mss = float(row.get('busted_mss_pval', 1.0))
                
                # Parameters BUSTED
                w3_b = float(row.get('busted_omega_class_2_omega', 0))
                p3_b = float(row.get('busted_omega_class_2_proportion', 0))
                aic_b = float(row.get('busted_AIC', 0))
                
                # Parameters MSS
                w3_m = float(row.get('busted_mss_omega_class_2_omega', 0))
                p3_m = float(row.get('busted_mss_omega_class_2_proportion', 0))
                aic_m = float(row.get('busted_mss_AIC', 0))
                
                delta_aic = aic_b - aic_m
                
                sig_b = p_busted <= 0.05
                sig_m = p_mss <= 0.05
                
                # Classification
                if sig_b and not sig_m:
                    if w3_b <= 100: # Filter outliers based on the source model
                        busted_only_w3.append(w3_b)
                        busted_only_p3.append(p3_b)
                        busted_only_aic.append(delta_aic)
                        
                elif not sig_b and sig_m:
                    if w3_m <= 100:
                        mss_only_w3.append(w3_m)
                        mss_only_p3.append(p3_m)
                        mss_only_aic.append(delta_aic)
                        
                elif sig_b and sig_m:
                    if w3_m <= 100: # Use MSS params for Both
                        both_w3.append(w3_m)
                        both_p3.append(p3_m)
                        both_aic.append(delta_aic)
                    
            except (ValueError, TypeError):
                continue

    # Process and Print
    groups = [
        ("BUSTED_Only", busted_only_w3, busted_only_p3, busted_only_aic),
        ("MSS_Only", mss_only_w3, mss_only_p3, mss_only_aic),
        ("Both", both_w3, both_p3, both_aic)
    ]

    # Calculate proportions as percentages
    # Need total counts. Hardcoding them based on Table 2 for consistency.
    totals = {
        "Caenorhabditis": 7696,
        "Drosophila": 12062,
        "Enterobacteria": 1594,
        "Primate": 14371,
        "Saccharomyces": 4916
    }
    
    for name, w_list, p_list, aic_list in groups:
        n, mw, minw, maxw = get_stats(w_list)
        _, mp, minp, maxp = get_stats(p_list)
        _, maic, minaic, maxaic = get_stats(aic_list)
        
        pct = (n / totals.get(clade, 1)) * 100
        
        n_str = f"{n} ({pct:.1f}\\%)" # Latex formatted percentage
        w_str = f"{mw:.2f} ({minw:.0f}-{maxw:.0f})"
        p_str = f"{mp:.3f} ({minp:.2g}-{maxp:.2g})"
        aic_str = f"{maic:.1f} ({minaic:.0f}-{maxaic:.0f})"
        
        print(f"{clade:<15} | {name:<12} | {n_str:<12} | {w_str:<15} | {p_str:<15} | {aic_str:<15}")
    print("-" * 90)