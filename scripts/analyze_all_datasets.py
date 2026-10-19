
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os

datasets = {
    'Drosophila': 'data/20250415_BUSTED_vs_BUSTED-MSS_codonGA_droso.csv',
    'Caenorhabditis': 'data/20250819_BUSTED_vs_BUSTED-MSS_codonGA_caenor.csv',
    'Enterobacteria': 'data/20250819_BUSTED_vs_BUSTED-MSS_codonGA_entero.csv',
    'Saccharomyces': 'data/20250819_BUSTED_vs_BUSTED-MSS_codonGA_sacch.csv',
    'Primate': 'data/20250820_BUSTED_vs_BUSTED-MSS_codonGA_primate.csv'
}

print(f"{'Dataset':<15} | {'N':<5} | {'rho (Omega)':<12} | {'p (Omega)':<10} | {'rho (Sig)':<12} | {'p (Sig)':<10}")
print("-" * 85)

for name, filepath in datasets.items():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue
        
    df = pd.read_csv(filepath)
    
    # Calculate Delta AIC (MSS improvement)
    df['Delta_AIC'] = df['busted_AIC'] - df['busted_mss_AIC']
    
    # Filter for genes where BUSTED found evidence of selection (p <= 0.05)
    significant_busted = df[df['busted_pval'] <= 0.05].copy()
    
    if len(significant_busted) < 10:
        print(f"{name:<15} | {len(significant_busted):<5} | {'N/A':<12} | {'N/A':<10} | {'N/A':<12} | {'N/A':<10}")
        continue

    # Calculate shift in omega_3 (assuming class 2 is the selection class)
    # Note: Column names are 'busted_omega_class_2_omega' etc.
    # omega_3 is usually class 2 (0, 1, 2)
    significant_busted['Omega_Diff'] = significant_busted['busted_omega_class_2_omega'] - significant_busted['busted_mss_omega_class_2_omega']

    # Also calculate shift in p-value (log scale)
    significant_busted['Log_P_Diff'] = -np.log10(significant_busted['busted_pval'] + 1e-10) - (-np.log10(significant_busted['busted_mss_pval'] + 1e-10))

    # Correlation
    corr_omega, p_omega = spearmanr(significant_busted['Delta_AIC'], significant_busted['Omega_Diff'])
    corr_p, p_p = spearmanr(significant_busted['Delta_AIC'], significant_busted['Log_P_Diff'])
    
    print(f"{name:<15} | {len(significant_busted):<5} | {corr_omega:<12.3f} | {p_omega:<10.2e} | {corr_p:<12.3f} | {p_p:<10.2e}")
