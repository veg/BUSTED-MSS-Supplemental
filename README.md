# BUSTED+S+MSS Supplemental Materials and Data

This repository accompanies the manuscript:

> **Accounting for Global Synonymous Substitution-Rate Heterogeneity Improves the Calibration of Episodic Positive Selection Inference**  
> Hannah E. Verdonk, Arthur M. Pivirotto, Jody Hey, and Sergei L. Kosakovsky Pond

It provides the clade-specific two-class codon-pair synonymous rate partitions (`--mss-file` TSVs and raw genetic algorithm JSON outputs), training gene lists, per-gene model fit summary tables, simulation benchmarks, and analysis/workflow scripts used in the manuscript.

---

## Repository Structure

### 1. `models/` — Clade-Specific 2-Class MSS Codon-Pair Partitions
Pre-computed 2-class (`NEUTRAL` vs. `SELECTED`) codon-pair synonymous rate partitions inferred via the HyPhy genetic algorithm (`MSS-selector-2.bf` / `MSS-selector-codon.bf`) and model-averaged (`processor-codon.bf`):
- `caenor-codonGA-2class-subsample1.tsv` / `.json` (*Caenorhabditis*)
- `droso-codonGA-2class-subsample1.tsv` / `.json` (*Drosophila*)
- `enterobacteria-codonGA-2class-subsample1.tsv` / `.json` (*Enterobacteria*)
- `primate-codonGA-2class-subsample1.tsv` / `.json` (*Primate*)
- `saccharomyces-codonGA-2class-subsample1.tsv` / `.json` (*Saccharomyces*)

### 2. `training_genes/` — Training Alignment Lists
Lists of randomly sampled alignments used to train the clade-specific GA partitions (`*-codonGA-2class-subsample1.txt`, 100 alignments per clade) and the 500-gene joint-fit reference sets (`*_500RandomAlign_0.txt`, 500 alignments per clade):
- `caenor-codonGA-2class-subsample1.txt` & `caenorhabditis_500RandomAlign_0.txt`
- `droso-codonGA-2class-subsample1.txt` & `drosophila_500RandomAlign_0.txt`
- `enterobacteria-codonGA-2class-subsample1.txt` & `enterobacteria_500RandomAlign_0.txt`
- `primate-codonGA-2class-subsample1.txt` & `primate_500RandomAlign_0.txt`
- `saccharomyces-codonGA-2class-subsample1.txt` & `saccharomyces_500RandomAlign_0.txt`

### 3. `results/` — Per-Gene Empirical Fits, Interaction Decomposition, and Simulations
- **Empirical BUSTED+S vs. BUSTED+S+MSS per-gene comparisons across all five clades:**
  - `20250415_BUSTED_vs_BUSTED-MSS_codonGA_droso.csv` (*Drosophila*, $N = 12,061$)
  - `20250819_BUSTED_vs_BUSTED-MSS_codonGA_caenor.csv` (*Caenorhabditis*, $N = 7,695$)
  - `20250819_BUSTED_vs_BUSTED-MSS_codonGA_entero.csv` (*Enterobacteria*, $N = 1,594$)
  - `20250819_BUSTED_vs_BUSTED-MSS_codonGA_sacch.csv` (*Saccharomyces*, $N = 4,916$)
  - `20250820_BUSTED_vs_BUSTED-MSS_codonGA_primate.csv` (*Primate*, $N = 14,371$)
- **2x2 Factorial Model Comparison (BUSTED, BUSTED+S, BUSTED+MSS, BUSTED+S+MSS) in *Enterobacteria*:**
  - `20251013_entero_busted_MSS_SRV_interactions.csv`
- **Neutral Codon-Selection False Positive Rate Simulation Grid ($dS_s \in [0.05, 1.0]$):**
  - `20250917_droso_2classCodonGA1aln_neutral_grid-dSs.csv`

### 4. `alignments/` — Empirical Codon Alignments and Phylogenetic Trees
Compressed archives containing the filtered coding sequence alignments (`.msa`) and RAxML-NG maximum-likelihood phylogenetic trees (`.msa.raxml.bestTree`) for each clade:
- `caenorhabditis-data.tar.gz`
- `drosophila-data.tar.gz`
- `enterobacteria-data.tar.gz`
- `primate-data.tar.gz`
- `saccharomyces-data.tar.gz`

### 5. `scripts/` — Pipeline Workflows and Statistical Analysis Scripts
- `snakefile-GA`, `snakefile-BUSTED`, `snakefile-BUSTED-MSS`: Snakemake workflows for running the GA partition search, BUSTED+S, and BUSTED+S+MSS.
- `processor-codon.bf`: HyPhy batch script to compute Akaike/BIC-weighted model-averaged codon-pair partitions (`.tsv`) from raw GA `.json` outputs.
- `compare-BUSTED-vs-BUSTEDMSS.py`: Extracts per-gene log-likelihoods, $\text{AIC}_c$, $\omega$ distributions, branch lengths, and LRT $p$-values from BUSTED+S and BUSTED+S+MSS JSON outputs.
- `MSS_SRV_interaction.py` & `compute_interaction_stats.py`: Computes the log-likelihood interaction decomposition $I$ across the 4-model factorial design.
- `calculate_table3_stats.py`, `analyze_all_datasets.py`, `analyze_all_logistic.py`: Statistical summary and regression scripts used to generate the tables and figures in the manuscript.

---

## Running BUSTED+S+MSS in HyPhy

To run **BUSTED+S+MSS** on a codon alignment and tree using one of the pre-computed clade partitions in `models/`:

```bash
hyphy busted \
  --alignment <path/to/alignment.msa> \
  --tree <path/to/tree.nwk> \
  --srv Yes \
  --error-sink Yes \
  --mss Yes \
  --mss-type Codon-file \
  --mss-file models/droso-codonGA-2class-subsample1.tsv \
  --mss-neutral NEUTRAL \
  --starting-points 5 \
  --output <path/to/output.BUSTED-MSS.json>
```
