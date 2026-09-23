"""
phylograph.py
=============
Single source of truth for stellar phylogenetic tree construction and
spectral comparison. Every run (MW-GSE, MW-Sgr, GSE-Sgr, Nissen) imports
the SAME functions from here, so the tree-building and distance code can
never drift out of sync again. Only per-run differences (which frames,
cleaning, sizes, intra mechanism) live in the notebook via RunConfig.

Cleaning (FE_H cuts, KDE reweighting, M54 removal, etc.) is intentionally
NOT in this module: it legitimately differs per comparison and stays in the
notebook. This module receives ready-to-sample DataFrames.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from collections import deque

import numpy as np
import pandas as pd
import networkx as nx
from ete3 import Tree
from biotite.sequence.phylo import neighbor_joining


# ----------------------------------------------------------------------
# Core: distance matrix, NJ tree, graph, spectra  (IDENTICAL FOR ALL RUNS)
# ----------------------------------------------------------------------

def chemical_distance(X: np.ndarray, ord: int = 1) -> np.ndarray:
    """Pairwise distance matrix in abundance space (ord=1 -> Manhattan)."""
    diff = np.abs(X[:, None, :] - X[None, :, :]) ** ord
    return (diff.sum(axis=2)) ** (1.0 / ord)


def nj_newick(X: np.ndarray, ord: int = 1) -> str:
    """Neighbor-Joining tree (Newick) from an abundance matrix."""
    return neighbor_joining(chemical_distance(X, ord=ord)).to_newick(include_distance=True)


def assign_names_to_nodes(tree: Tree) -> Tree:
    c = 1
    queue = deque([tree])
    while queue:
        node = queue.popleft()
        if not node.name:
            node.name = "root" if node.is_root() else f"node{c}"
            c += 1
        queue.extend(node.children)
    return tree
def branch_length_distribution(newick_file: str, max_trees: int = 50) -> np.ndarray:
    """Collect all branch lengths from the first max_trees trees in a newick file."""
    pairs = _read_pairs(newick_file)[:max_trees]
    lengths = []
    for t1_str, t2_str in pairs:
        for newick in (t1_str, t2_str):
            t = Tree(newick.strip())
            for node in t.traverse():
                if not node.is_root():
                    lengths.append(node.dist)
    return np.array(lengths)
# ----------------------------------------------------------------------
# Edge-weight kernels: choose how branch length d maps to graph weight
# ----------------------------------------------------------------------

def affinity_kernel(sigma: float = 1.0):
    """w = exp(-d/sigma).  The e^{-d} convention adopted in the paper."""
    return lambda d: np.exp(-d / sigma)

def distance_kernel():
    """w = d.  The alternative A_ij = d_ij convention (divergence)."""
    return lambda d: d


def tree_to_graph(tree: Tree, G: nx.Graph,
                  weight_fn=None,
                  clip_negative: bool = True) -> nx.Graph:
    """Weighted graph from an ete3 tree.

    weight_fn : callable d -> w mapping branch length to edge weight.
                Defaults to affinity_kernel(1.0) (the paper's e^{-d}).
                Pass distance_kernel() for the A_ij = d_ij convention.
    clip_negative : if True, NJ branch lengths < 0 are set to 0 before
                    weighting (keeps affinity weights in (0,1]); if False,
                    a warning is issued and the raw value is used.
    """
    if weight_fn is None:
        weight_fn = affinity_kernel(1.0)
    if tree.is_leaf():
        return G
    for child in tree.children:
        d = tree.get_distance(child)
        if d < 0:
            if clip_negative:
                d = 0.0
            else:
                warnings.warn(f"negative branch length {d:.4g}; weight may exceed 1")
        G.add_edge(tree.name, child.name, weight=weight_fn(d))
        tree_to_graph(child, G, weight_fn=weight_fn, clip_negative=clip_negative)
    return G


def _graph_from_newick(newick, weight_fn=None, clip_negative=True):
    t = assign_names_to_nodes(Tree(newick.strip()))
    return tree_to_graph(t, nx.Graph(), weight_fn=weight_fn, clip_negative=clip_negative)


def _spectral_distance(eigs1, eigs2):
    k = min(len(eigs1), len(eigs2))
    return float(np.linalg.norm(eigs1[:k] - eigs2[:k]))


def adjacency_spectral_distance(G1, G2):
    # adjacency: sort DESCENDING (largest eigenvalues first)
    e1 = np.sort(np.linalg.eigvalsh(nx.to_numpy_array(G1)))[::-1]
    e2 = np.sort(np.linalg.eigvalsh(nx.to_numpy_array(G2)))[::-1]
    return _spectral_distance(e1, e2)


def laplacian_spectral_distance(G1, G2):
    # combinatorial Laplacian: sort ASCENDING
    e1 = np.sort(np.linalg.eigvalsh(nx.laplacian_matrix(G1).toarray()))
    e2 = np.sort(np.linalg.eigvalsh(nx.laplacian_matrix(G2).toarray()))
    return _spectral_distance(e1, e2)


def normalized_laplacian_distance(G1, G2):
    e1 = np.sort(np.linalg.eigvalsh(nx.normalized_laplacian_matrix(G1).toarray()))
    e2 = np.sort(np.linalg.eigvalsh(nx.normalized_laplacian_matrix(G2).toarray()))
    return _spectral_distance(e1, e2)


def pair_distances(newick1, newick2, weight_fn=None, clip_negative=True):
    G1 = _graph_from_newick(newick1, weight_fn, clip_negative)
    G2 = _graph_from_newick(newick2, weight_fn, clip_negative)
    return (adjacency_spectral_distance(G1, G2),
            laplacian_spectral_distance(G1, G2),
            normalized_laplacian_distance(G1, G2))



# ----------------------------------------------------------------------
# Run configuration  (the ONLY thing that changes between comparisons)
# ----------------------------------------------------------------------

@dataclass
class RunConfig:
    name: str
    pop_a: tuple
    pop_b: tuple
    elements: list
    out_dir: str = "data/history_trees"
    n_nodes: int = 120
    n_iterations: int = 200
    dist_ord: int = 1
    weight_fn: callable = field(default_factory=lambda: affinity_kernel(1.0))
    clip_negative: bool = True
    perturb_sigma: float = 0.0
    shared_base_intra: bool = False
    seed: int = 0

    def frames(self):
        (la, da), (lb, db) = self.pop_a, self.pop_b
        da = da[self.elements].dropna().reset_index(drop=True)
        db = db[self.elements].dropna().reset_index(drop=True)
        print(f"[{self.name}] {la}: {len(da)} stars | {lb}: {len(db)} stars")
        return (la, da), (lb, db)


# ----------------------------------------------------------------------
# Tree generation  (driven entirely by RunConfig)
# ----------------------------------------------------------------------

def _sample_matrix(df, n, rng):
    idx = rng.choice(len(df), size=n, replace=False)
    return df.values[idx]

def _two_matrices(dfa, dfb, intra, cfg, rng):
    """Sample (and optionally perturb) the two abundance matrices for one
    iteration. Extracted from _two_trees so the abundance baseline can be
    computed from the SAME matrices that build the trees."""
    n = cfg.n_nodes
    Xa = _sample_matrix(dfa, n, rng)
    if intra and cfg.shared_base_intra:
        Xb = Xa.copy()
    else:
        Xb = _sample_matrix(dfb, n, rng)
    if cfg.perturb_sigma > 0:
        Xa = Xa + rng.normal(0, cfg.perturb_sigma, Xa.shape)
        Xb = Xb + rng.normal(0, cfg.perturb_sigma, Xb.shape)
    return Xa, Xb


def _two_trees(dfa, dfb, intra, cfg, rng):
    Xa, Xb = _two_matrices(dfa, dfb, intra, cfg, rng)
    return (nj_newick(Xa, ord=cfg.dist_ord),
            nj_newick(Xb, ord=cfg.dist_ord))


def generate_trees(cfg: RunConfig) -> dict:
    """Write {A-A, B-B, A-B}.newick files. Returns {pair_name: path}."""
    os.makedirs(cfg.out_dir, exist_ok=True)
    (la, dfa), (lb, dfb) = cfg.frames()
    if len(dfa) < cfg.n_nodes or len(dfb) < cfg.n_nodes:
        raise ValueError(f"need >= {cfg.n_nodes} stars; have {len(dfa)} / {len(dfb)}")

    pairs = [
        (f"{la}-{la}", dfa, dfa, True),
        (f"{lb}-{lb}", dfb, dfb, True),
        (f"{la}-{lb}", dfa, dfb, False),
    ]
    rng = np.random.default_rng(cfg.seed)
    paths = {}
    for pair_name, d1, d2, intra in pairs:
        path = os.path.join(cfg.out_dir, f"{pair_name}.newick")
        with open(path, "w") as f:
            for _ in range(cfg.n_iterations):
                t1, t2 = _two_trees(d1, d2, intra, cfg, rng)
                f.write(f"{t1}\n{t2}\n")
        paths[pair_name] = path
    return paths


# ----------------------------------------------------------------------
# Distance computation over saved trees
# ----------------------------------------------------------------------

def _read_pairs(path):
    with open(path) as f:
        lines = f.readlines()
    return [(lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]


def compute_distances(cfg: RunConfig, paths: dict | None = None,
                      tree_dir: str | None = None) -> dict:
    """Return {pair_name: DataFrame[ASD, LSD, NLSD]} and save CSVs.

    tree_dir : directory to READ .newick files from. Defaults to cfg.out_dir.
               Set it to reuse one set of trees while writing CSVs elsewhere
               (e.g. the same trees under a different weight kernel).
    """
    read_dir = tree_dir if tree_dir is not None else cfg.out_dir
    if paths is None:
        (la, _), (lb, _) = cfg.pop_a, cfg.pop_b
        names = [f"{la}-{la}", f"{lb}-{lb}", f"{la}-{lb}"]
        paths = {n: os.path.join(read_dir, f"{n}.newick") for n in names}

    os.makedirs(cfg.out_dir, exist_ok=True)
    results = {}
    for pair_name, path in paths.items():
        rows = [pair_distances(t1, t2,
                               weight_fn=cfg.weight_fn,
                               clip_negative=cfg.clip_negative)
                for t1, t2 in _read_pairs(path)]
        df = pd.DataFrame(rows, columns=["ASD", "LSD", "NLSD"])
        df.to_csv(os.path.join(cfg.out_dir, f"{pair_name}_distances.csv"), index=False)
        results[pair_name] = df
    return results

# ----------------------------------------------------------------------
# AUC separation metric  (rank-based; invariant to per-metric rescaling)
# ----------------------------------------------------------------------
def auc_table(panels, measures=("ASD", "LSD", "NLSD"), suffix="distances"):
    """panels: {panel_label: (dir, [intra_csv_names], inter_csv_name)}.

    AUC = P(D_inter > D_intra), the normalized Mann-Whitney U
    (Hanley & McNeil 1982). 1 = clean separation, 0.5 = none, <0.5 = inversion.
    """
    from sklearn.metrics import roc_auc_score

    def load(d, n):
        return pd.read_csv(os.path.join(d, f"{n}_{suffix}.csv"))

    rows = []
    for panel, (d, intra_names, inter_name) in panels.items():
        intra = pd.concat([load(d, n) for n in intra_names], ignore_index=True)
        inter = load(d, inter_name)
        for m in measures:
            y = np.r_[np.zeros(len(intra)), np.ones(len(inter))]
            s = np.r_[intra[m].values, inter[m].values]
            rows.append({"panel": panel, "measure": m, "AUC": roc_auc_score(y, s)})
    return (pd.DataFrame(rows)
            .pivot(index="measure", columns="panel", values="AUC")
            .reindex(index=list(measures)))

def auc_table_latex(panels: dict, measures=("ASD", "LSD", "NLSD")) -> str:
    """Return a copy-pasteable A&A LaTeX table from the AUC results."""
    df = auc_table(panels, measures=measures).round(2)

    # panel display names (unicode dashes -> LaTeX)
    col_labels = {c: c.replace("-", "--") for c in df.columns}

    n_cols = len(df.columns)
    col_fmt = "l" + "c" * n_cols

    lines = []
    lines.append(r"\begin{table}")
    lines.append(r"\caption{AUC values for the separation between intra- and")
    lines.append(r"inter-population distance distributions.}")
    lines.append(r"\label{tab:auc}")
    lines.append(r"\centering")
    lines.append(r"\begin{tabular}{" + col_fmt + r"}")
    lines.append(r"\hline\hline")

    # header row
    header = "Measure & " + " & ".join(col_labels.get(c, c) for c in df.columns) + r" \\"
    lines.append(header)
    lines.append(r"\hline")

    # data rows
    for measure, row in df.iterrows():
        vals = " & ".join(f"{v:.2f}" for v in row.values)
        lines.append(f"{measure} & {vals}" + r" \\")



    return "\n".join(lines)