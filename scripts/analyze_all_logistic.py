
import pandas as pd
import numpy as np
import statsmodels.api as sm
import os

datasets = {
    'Drosophila': 'data/20250415_BUSTED_vs_BUSTED-MSS_codonGA_droso.csv',
    'Caenorhabditis': 'data/20250819_BUSTED_vs_BUSTED-MSS_codonGA_caenor.csv',
    'Enterobacteria': 'data/20250819_BUSTED_vs_BUSTED-MSS_codonGA_entero.csv',
    'Saccharomyces': 'data/20250819_BUSTED_vs_BUSTED-MSS_codonGA_sacch.csv',
    'Primate': 'data/20250820_BUSTED_vs_BUSTED-MSS_codonGA_primate.csv'
}

print(f"{'Dataset':<15} | {'N':<5} | {'Lost':<5} | {'DeltaAIC_Coef':<12} | {'DeltaAIC_P':<10} | {'LRT_Coef':<10} | {'LRT_P':<10}")
print("-" * 90)

for name, filepath in datasets.items():
    if not os.path.exists(filepath):
        continue
        
    df = pd.read_csv(filepath)
    
    # Calculate Delta AIC
    df['Delta_AIC'] = df['busted_AIC'] - df['busted_mss_AIC']
    
    # Filter for significant BUSTED genes
    initial_sig = df[df['busted_pval'] <= 0.05].copy()
    
    if len(initial_sig) < 20: # Need sufficient samples for logistic regression
        print(f"{name:<15} | {len(initial_sig):<5} | {'N/A':<5} | {'N/A':<12} | {'N/A':<10} | {'N/A':<10} | {'N/A':<10}")
        continue

    # Define Target: Lost Significance (1 if MSS p > 0.05)
    initial_sig['Lost_Signal'] = (initial_sig['busted_mss_pval'] > 0.05).astype(int)
    
    lost_count = initial_sig['Lost_Signal'].sum()
    if lost_count < 5 or (len(initial_sig) - lost_count) < 5:
         print(f"{name:<15} | {len(initial_sig):<5} | {lost_count:<5} | {'Skewed':<12} | {'N/A':<10} | {'N/A':<10} | {'N/A':<10}")
         continue

    # Predictors
    initial_sig['Log_LRT'] = np.log(initial_sig['busted_LRT'] + 1e-10)
    initial_sig['Omega3_BUSTED'] = initial_sig['busted_omega_class_2_omega']
    initial_sig['p3_BUSTED'] = initial_sig['busted_omega_class_2_proportion']
    
    # Features
    features = ['Delta_AIC', 'Log_LRT', 'Omega3_BUSTED', 'p3_BUSTED']
    X = initial_sig[features]
    X = sm.add_constant(X)
    y = initial_sig['Lost_Signal']
    
    try:
        model = sm.Logit(y, X)
        result = model.fit(disp=0)
        
        aic_coef = result.params['Delta_AIC']
        aic_p = result.pvalues['Delta_AIC']
        lrt_coef = result.params['Log_LRT']
        lrt_p = result.pvalues['Log_LRT']
        
        print(f"{name:<15} | {len(initial_sig):<5} | {lost_count:<5} | {aic_coef:<12.4f} | {aic_p:<10.3e} | {lrt_coef:<10.4f} | {lrt_p:<10.3e}")
        
    except Exception as e:
        print(f"{name:<15} | Error: {e}")
