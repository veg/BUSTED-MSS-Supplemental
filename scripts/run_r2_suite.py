#!/usr/bin/env python3
"""
Generate all 800 simulated codon alignments (8 scenarios x 100 empirical 14-species Enterobacteria trees)
and provide a worker function to run paired BUSTED+S and BUSTED+S+MSS fits using hyphy3-busted.
"""

import os
import re
import sys
import json
import subprocess
import numpy as np
from scipy.linalg import expm

NUC = ['A', 'C', 'G', 'T']
CODON_TABLE = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'TGT':'C', 'TGC':'C', 'TGA':'*',
    'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R',
    'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
    'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G',
    'GGG':'G'
}
SENSE_CODONS = [c1+c2+c3 for c1 in NUC for c2 in NUC for c3 in NUC if CODON_TABLE[c1+c2+c3] != '*']
AMINO_ACIDS = sorted(list(set(CODON_TABLE[c] for c in SENSE_CODONS)))
AA_TO_IDX = {a: i for i, a in enumerate(AMINO_ACIDS)}
CODON_AA_IDX = np.array([AA_TO_IDX[CODON_TABLE[c]] for c in SENSE_CODONS])

OPTIMAL_CODONS = {
    'AAA', 'AAC', 'ACC', 'ACG', 'AGC', 'ATC', 'CAG', 'CCG', 'CGC', 'CGT',
    'CTG', 'GAA', 'GAC', 'GCC', 'GCG', 'GGT', 'GGC', 'GTT', 'GTA', 'TCT', 'TCC', 'TTC', 'TAC'
}

AA_HYDRO = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}
AA_VOL = {
    'A': 88.6, 'C': 108.5, 'D': 111.1, 'E': 138.4, 'F': 189.9,
    'G': 60.1, 'H': 153.2, 'I': 166.7, 'K': 168.6, 'L': 166.7,
    'M': 162.9, 'N': 114.1, 'P': 112.7, 'Q': 143.8, 'R': 173.4,
    'S': 89.0, 'T': 116.1, 'V': 140.0, 'W': 227.8, 'Y': 193.6
}
H_VEC = np.array([AA_HYDRO[a] for a in AMINO_ACIDS])
H_VEC = (H_VEC - H_VEC.mean()) / H_VEC.std()
V_VEC = np.array([AA_VOL[a] for a in AMINO_ACIDS])
V_VEC = (V_VEC - V_VEC.mean()) / V_VEC.std()


def parse_newick(newick_str):
    s = newick_str.strip().rstrip(';')
    tokens = re.findall(r'\(|\)|,|[^(),:;]+|:[^(),:;]+', s)
    nodes = [{'name': 'root', 'parent': None, 'length': 0.0, 'children': []}]
    curr = 0
    for tok in tokens:
        if tok == '(':
            new_idx = len(nodes)
            nodes.append({'name': f'internal_{new_idx}', 'parent': curr, 'length': 0.0, 'children': []})
            nodes[curr]['children'].append(new_idx)
            curr = new_idx
        elif tok == ',':
            curr = nodes[curr]['parent']
            new_idx = len(nodes)
            nodes.append({'name': f'internal_{new_idx}', 'parent': curr, 'length': 0.0, 'children': []})
            nodes[curr]['children'].append(new_idx)
            curr = new_idx
        elif tok == ')':
            curr = nodes[curr]['parent']
        elif tok.startswith(':'):
            nodes[curr]['length'] = max(float(tok[1:]), 1e-7)
        else:
            nodes[curr]['name'] = tok.strip("'")
    return nodes


def load_mss_partition(tsv_path):
    selected_pairs = set()
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[3] in ('SELECTED', 'NEUTRAL'):
                aa, c1, c2, cls = parts[0], parts[1], parts[2], parts[3]
                if cls == 'SELECTED':
                    selected_pairs.add((c1, c2))
                    selected_pairs.add((c2, c1))
    return selected_pairs


def build_nucleotide_mutation_matrix():
    pi_nuc = np.array([0.24, 0.26, 0.27, 0.23])
    theta = {
        ('A','C'): 0.35, ('A','G'): 1.00, ('A','T'): 0.30,
        ('C','G'): 0.40, ('C','T'): 1.20, ('G','T'): 0.35
    }
    M = np.zeros((4, 4))
    for i, n1 in enumerate(NUC):
        for j, n2 in enumerate(NUC):
            if i != j:
                r = theta.get((n1, n2), theta.get((n2, n1)))
                M[i, j] = r * pi_nuc[j]
    M -= np.diag(np.sum(M, axis=1))
    M /= (-np.sum(pi_nuc * np.diag(M)))
    pi_codon_mut = np.zeros(61)
    for idx, c in enumerate(SENSE_CODONS):
        pi_codon_mut[idx] = pi_nuc[NUC.index(c[0])] * pi_nuc[NUC.index(c[1])] * pi_nuc[NUC.index(c[2])]
    pi_codon_mut /= np.sum(pi_codon_mut)
    return M, pi_nuc, pi_codon_mut


def fix_prob(s):
    if abs(s) < 1e-6:
        return 1.0 + 0.5 * s + (s * s) / 12.0
    return s / (1.0 - np.exp(-s))


def simulate_mutsel_alignment(nodes, L, selected_pairs, S_syn, kappa_syn_sel, rng):
    M_nuc, pi_nuc, pi_codon_mut = build_nucleotide_mutation_matrix()
    F_syn = np.array([S_syn if c in OPTIMAL_CODONS else 0.0 for c in SENSE_CODONS])

    n_archetypes = 25
    srv_rates = [0.25, 0.85, 2.15]
    srv_probs = [0.25, 0.55, 0.20]

    Q_archetypes = []
    pi_archetypes = []
    for _ in range(n_archetypes):
        w_h = rng.normal(0.0, 2.2)
        w_v = rng.normal(0.0, 2.2)
        pref_aa = rng.integers(0, 20)
        F_aa = w_h * H_VEC + w_v * V_VEC + rng.normal(0.0, 1.2, size=20)
        F_aa[pref_aa] += rng.uniform(3.0, 6.0)
        F_codon = F_aa[CODON_AA_IDX] + F_syn

        log_pi = np.log(pi_codon_mut) + F_codon
        log_pi -= np.max(log_pi)
        pi_c = np.exp(log_pi)
        pi_c /= np.sum(pi_c)

        Q = np.zeros((61, 61))
        for i, c1 in enumerate(SENSE_CODONS):
            for j, c2 in enumerate(SENSE_CODONS):
                if i == j:
                    continue
                diffs = [pos for pos in range(3) if c1[pos] != c2[pos]]
                if len(diffs) == 1:
                    pos = diffs[0]
                    mu_ij = M_nuc[NUC.index(c1[pos]), NUC.index(c2[pos])]
                    s_ij = F_codon[j] - F_codon[i]
                    rate = mu_ij * fix_prob(s_ij)
                    if CODON_TABLE[c1] == CODON_TABLE[c2] and (c1, c2) in selected_pairs:
                        rate *= kappa_syn_sel
                    Q[i, j] = rate
        Q -= np.diag(np.sum(Q, axis=1))
        Q *= 3.0
        Q_archetypes.append(Q)
        pi_archetypes.append(pi_c)

    site_arch = rng.integers(0, n_archetypes, size=L)
    site_srv = rng.choice(3, size=L, p=srv_probs)

    node_states = np.zeros((len(nodes), L), dtype=np.int32)
    for s in range(L):
        node_states[0, s] = rng.choice(61, p=pi_archetypes[site_arch[s]])

    groups = {}
    for s in range(L):
        groups.setdefault((site_arch[s], site_srv[s]), []).append(s)

    for (arch_idx, srv_idx), site_indices in groups.items():
        Q_scaled = Q_archetypes[arch_idx] * srv_rates[srv_idx]
        for n_idx in range(1, len(nodes)):
            parent_idx = nodes[n_idx]['parent']
            bl = nodes[n_idx]['length']
            P = expm(Q_scaled * bl)
            P = np.maximum(P, 0.0)
            P /= np.sum(P, axis=1, keepdims=True)
            for s in site_indices:
                node_states[n_idx, s] = rng.choice(61, p=P[node_states[parent_idx, s]])

    return {
        n['name']: ''.join(SENSE_CODONS[node_states[n_idx, s]] for s in range(L))
        for n_idx, n in enumerate(nodes) if len(n['children']) == 0
    }


def build_mg94_mss_matrix(M_nuc, selected_pairs, dS_selected, omega):
    Q = np.zeros((61, 61))
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            if i == j:
                continue
            diffs = [pos for pos in range(3) if c1[pos] != c2[pos]]
            if len(diffs) == 1:
                pos = diffs[0]
                rate = M_nuc[NUC.index(c1[pos]), NUC.index(c2[pos])]
                if CODON_TABLE[c1] == CODON_TABLE[c2]:
                    if (c1, c2) in selected_pairs:
                        rate *= dS_selected
                else:
                    rate *= omega
                Q[i, j] = rate
    Q -= np.diag(np.sum(Q, axis=1))
    return Q


def simulate_mg94_busted_alignment(nodes, L, selected_pairs, dS_selected, omegas, omega_probs, rng):
    M_nuc, pi_nuc, pi_codon = build_nucleotide_mutation_matrix()
    srv_rates = [0.25, 0.85, 2.15]
    srv_probs = [0.25, 0.55, 0.20]

    Q_ref = build_mg94_mss_matrix(M_nuc, selected_pairs, 1.0, 1.0)
    ref_flux = -np.sum(pi_codon * np.diag(Q_ref))
    scale_factor = 3.0 / ref_flux

    Q_by_omega = [build_mg94_mss_matrix(M_nuc, selected_pairs, dS_selected, w) * scale_factor for w in omegas]

    site_srv = rng.choice(3, size=L, p=srv_probs)
    node_states = np.zeros((len(nodes), L), dtype=np.int32)
    node_states[0, :] = rng.choice(61, size=L, p=pi_codon)

    for n_idx in range(1, len(nodes)):
        parent_idx = nodes[n_idx]['parent']
        bl = nodes[n_idx]['length']
        branch_omega_cls = rng.choice(3, size=L, p=omega_probs)
        for srv_idx in range(3):
            for w_idx in range(3):
                mask = np.where((site_srv == srv_idx) & (branch_omega_cls == w_idx))[0]
                if len(mask) == 0:
                    continue
                P = expm(Q_by_omega[w_idx] * (srv_rates[srv_idx] * bl))
                P = np.maximum(P, 0.0)
                P /= np.sum(P, axis=1, keepdims=True)
                for s in mask:
                    node_states[n_idx, s] = rng.choice(61, p=P[node_states[parent_idx, s]])

    return {
        n['name']: ''.join(SENSE_CODONS[node_states[n_idx, s]] for s in range(L))
        for n_idx, n in enumerate(nodes) if len(n['children']) == 0
    }


def write_fasta(seqs, path):
    with open(path, 'w') as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


SCENARIOS = [
    ("mutsel_aa_only",           "mutsel", {"S_syn": 0.0, "kappa_syn_sel": 1.00}),
    ("mutsel_aa_plus_syn",       "mutsel", {"S_syn": 1.2, "kappa_syn_sel": 0.45}),
    ("mg94_null_dS100",          "mg94",   {"dS": 1.00, "omegas": [0.02, 0.25, 1.00], "probs": [0.70, 0.25, 0.05]}),
    ("mg94_null_dS050",          "mg94",   {"dS": 0.50, "omegas": [0.02, 0.25, 1.00], "probs": [0.70, 0.25, 0.05]}),
    ("mg94_pos_w25_dS100",       "mg94",   {"dS": 1.00, "omegas": [0.02, 0.25, 2.50], "probs": [0.70, 0.26, 0.04]}),
    ("mg94_pos_w25_dS050",       "mg94",   {"dS": 0.50, "omegas": [0.02, 0.25, 2.50], "probs": [0.70, 0.26, 0.04]}),
    ("mg94_pos_w40_dS100",       "mg94",   {"dS": 1.00, "omegas": [0.02, 0.25, 4.00], "probs": [0.70, 0.25, 0.05]}),
    ("mg94_pos_w40_dS050",       "mg94",   {"dS": 0.50, "omegas": [0.02, 0.25, 4.00], "probs": [0.70, 0.25, 0.05]}),
]


def run_task(task_idx, base_dir, hyphy3_bin, tsv_path):
    """Run a single (scenario_idx, rep_idx) task in [0..799]."""
    scen_idx = task_idx // 100
    rep_idx = task_idx % 100
    scen_name, scen_type, params = SCENARIOS[scen_idx]

    with open(os.path.join(base_dir, "selected_100_trees.json")) as f:
        trees = json.load(f)
    t_info = trees[rep_idx]
    gene = t_info["gene"]
    L = t_info["sites"]
    nwk = t_info["newick"]

    scen_dir = os.path.join(base_dir, "runs", scen_name)
    os.makedirs(scen_dir, exist_ok=True)
    fna_path = os.path.join(scen_dir, f"rep_{rep_idx:03d}_{gene}.fna")
    nwk_path = os.path.join(scen_dir, f"rep_{rep_idx:03d}_{gene}.nwk")
    json_s = os.path.join(scen_dir, f"rep_{rep_idx:03d}_{gene}.busted_s.json")
    json_sm = os.path.join(scen_dir, f"rep_{rep_idx:03d}_{gene}.busted_s_mss.json")

    if os.path.exists(json_s) and os.path.exists(json_sm) and os.path.getsize(json_s) > 500 and os.path.getsize(json_sm) > 500:
        return

    nodes = parse_newick(nwk)
    selected_pairs = load_mss_partition(tsv_path)
    # Use paired seed across conditions so rep_idx has identical background seed!
    rng = np.random.default_rng(20261019 + rep_idx * 1000)

    if scen_type == "mutsel":
        seqs = simulate_mutsel_alignment(nodes, L, selected_pairs, params["S_syn"], params["kappa_syn_sel"], rng)
    else:
        seqs = simulate_mg94_busted_alignment(nodes, L, selected_pairs, params["dS"], params["omegas"], params["probs"], rng)

    write_fasta(seqs, fna_path)
    with open(nwk_path, "w") as f:
        f.write(nwk + "\n")

    # Fit BUSTED+S
    if not (os.path.exists(json_s) and os.path.getsize(json_s) > 500):
        subprocess.run([
            hyphy3_bin,
            "--alignment", fna_path,
            "--tree", nwk_path,
            "--srv",
            "--threads", "2",
            "--no-progress",
            "--output", json_s
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

    # Fit BUSTED+S+MSS
    if not (os.path.exists(json_sm) and os.path.getsize(json_sm) > 500):
        subprocess.run([
            hyphy3_bin,
            "--alignment", fna_path,
            "--tree", nwk_path,
            "--srv",
            "--mss", tsv_path,
            "--mss-neutral", "NEUTRAL",
            "--mss-norm", "neutral_ref",
            "--threads", "2",
            "--no-progress",
            "--output", json_sm
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)


if __name__ == "__main__":
    # Usage: run_r2_suite.py <base_dir> <hyphy3_bin> <tsv_path> <start_idx> <end_idx> <num_workers>
    from concurrent.futures import ProcessPoolExecutor
    base_dir = sys.argv[1]
    hyphy3_bin = sys.argv[2]
    tsv_path = sys.argv[3]
    start_idx = int(sys.argv[4])
    end_idx = int(sys.argv[5])
    num_workers = int(sys.argv[6])

    tasks = list(range(start_idx, end_idx))
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = [ex.submit(run_task, t, base_dir, hyphy3_bin, tsv_path) for t in tasks]
        for i, fut in enumerate(futures, 1):
            fut.result()
            if i % 10 == 0 or i == len(tasks):
                print(f"[{start_idx}..{end_idx}] Completed {i}/{len(tasks)}", flush=True)
