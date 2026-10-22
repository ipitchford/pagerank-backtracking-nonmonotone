"""Prevalence of non-monotone backtracking-weighted PageRank (NUMERICAL).

SPDX-License-Identifier: MIT
Copyright (c) 2026 Anonymous
Author: Anonymous

Model. G is a connected simple undirected graph with minimum degree >= 2
(we use the largest connected component of the 2-core). States are directed
edges (i,j). The column-stochastic kernel Q_mu moves (i,j) -> (j,k), k ~ j,
with weight (1 + (mu-1)[k == i]) / (d_j - 1 + mu). Restart
u_(i,j) = v_i / d_i with v uniform. Edge PageRank q solves
q = (1-eps) u + eps Q_mu q, and the vertex score is the tail projection
pi_mu(i) = sum_{j~i} q(i,j). Endpoints: pi_1 is ordinary PageRank
(1-eps)(I - eps A D^{-1})^{-1} v and pi_inf = (v + eps P v)/(1 + eps).

Usage: python3 -B -O prevalence.py OUTDIR   (writes OUTDIR/prevalence.json)
All results are numerical (float64, scipy.sparse.linalg.spsolve).
"""

import argparse
import json
import math
import os
from fractions import Fraction

import networkx as nx
import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.stats import kendalltau

EPS_LIST = (0.5, 0.85)
TOL_REL = 1e-12
MU_GRID = np.unique(np.concatenate([np.linspace(0.0, 1.0, 101),
                                    np.geomspace(1.0, 1000.0, 200)]))


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def r6(x):
    """Round to 6 significant figures; None for non-finite / None."""
    if x is None:
        return None
    x = float(x)
    if not math.isfinite(x):
        return None
    return float(f"{x:.6g}")


# ---------------------------------------------------------------- model
class EdgeModel:
    def __init__(self, G):
        self.G = G
        self.n = G.number_of_nodes()
        nodes = list(G.nodes())
        check(sorted(nodes) == list(range(self.n)), "nodes must be 0..n-1")
        self.deg = np.array([G.degree(i) for i in range(self.n)], dtype=float)
        check(self.deg.min() >= 2, "min degree must be >= 2")
        states = []
        for a, b in sorted(G.edges()):
            a, b = min(a, b), max(a, b)
            states.append((a, b))
            states.append((b, a))
        states.sort()
        self.states = states
        idx = {s: t for t, s in enumerate(states)}
        S = len(states)
        self.S = S
        rN, cN, rR, cR = [], [], [], []
        for t, (i, j) in enumerate(states):
            for k in sorted(G.neighbors(j)):
                tgt = idx[(j, k)]
                if k == i:
                    rR.append(tgt)
                    cR.append(t)
                else:
                    rN.append(tgt)
                    cN.append(t)
        self.N = sp.csc_matrix((np.ones(len(rN)), (rN, cN)), shape=(S, S))
        self.R = sp.csc_matrix((np.ones(len(rR)), (rR, cR)), shape=(S, S))
        tails = np.array([s[0] for s in states])
        heads = np.array([s[1] for s in states])
        self.tail_deg = self.deg[tails]
        self.head_deg = self.deg[heads]
        self.Tproj = sp.csr_matrix((np.ones(S), (tails, np.arange(S))),
                                   shape=(self.n, S))
        self.v = np.full(self.n, 1.0 / self.n)
        self.u = self.v[tails] / self.tail_deg
        self.I = sp.identity(S, format="csc")
        A = nx.to_scipy_sparse_array(G, nodelist=range(self.n), format="csr",
                                     dtype=float)
        self.A = sp.csr_matrix(A)
        self.P = sp.csr_matrix(self.A @ sp.diags(1.0 / self.deg))

    def Q(self, mu):
        denom = self.head_deg - 1.0 + mu
        return (self.N + mu * self.R) @ sp.diags(1.0 / denom)

    def pi(self, mu, eps):
        M = (self.I - eps * self.Q(mu)).tocsc()
        q = spla.spsolve(M, (1.0 - eps) * self.u)
        p = self.Tproj @ q
        check(abs(p.sum() - 1.0) < 1e-10, f"sum pi != 1 at mu={mu}")
        check(p.min() > -1e-14, "negative score")
        return p

    def pi_pagerank(self, eps):
        M = (sp.identity(self.n, format="csc") - eps * self.P).tocsc()
        return (1.0 - eps) * spla.spsolve(M, self.v)

    def pi_inf(self, eps):
        return (self.v + eps * (self.P @ self.v)) / (1.0 + eps)


def validate_house():
    G = nx.Graph([(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3)])
    m = EdgeModel(G)
    eps = 0.5
    p0 = m.pi(0.0, eps)
    p1 = m.pi(1.0, eps)
    pbig = m.pi(1e8, eps)
    pinf = m.pi_inf(eps)
    t0 = float(Fraction(638, 3615))
    t1 = float(Fraction(22, 125))
    ti = float(Fraction(8, 45))
    check(abs(p0[4] - t0) < 1e-12, f"house pi_0(4)={p0[4]} != 638/3615")
    check(abs(p1[4] - t1) < 1e-12, f"house pi_1(4)={p1[4]} != 22/125")
    check(abs(pinf[4] - ti) < 1e-14, "house pi_inf formula(4) != 8/45")
    check(abs(pbig[4] - ti) < 1e-6, f"house pi_1e8(4)={pbig[4]} !-> 8/45")
    check(np.max(np.abs(pbig - pinf)) < 1e-6, "house mu->inf limit mismatch")
    check(np.max(np.abs(p1 - m.pi_pagerank(eps))) < 1e-13,
          "house pi_1 != PageRank formula")
    return {
        "graph": "house edges 01,03,04,12,14,23",
        "eps": 0.5,
        "pi_0_vertex4": r6(p0[4]), "pi_0_vertex4_target": "638/3615",
        "pi_1_vertex4": r6(p1[4]), "pi_1_vertex4_target": "22/125",
        "pi_mu1e8_vertex4": r6(pbig[4]), "pi_inf_vertex4_target": "8/45",
        "abs_err_pi0": r6(abs(p0[4] - t0)),
        "abs_err_pi1": r6(abs(p1[4] - t1)),
        "abs_err_mu1e8_vs_inf": r6(abs(pbig[4] - ti)),
        "passed": True,
    }


# ---------------------------------------------------------------- graphs
def core_lcc(G):
    H = nx.Graph(G)
    H.remove_edges_from(list(nx.selfloop_edges(H)))
    C = nx.k_core(H, 2)
    if C.number_of_nodes() == 0:
        return None
    comps = sorted(nx.connected_components(C),
                   key=lambda c: (-len(c), sorted(c)[0]))
    L = nx.Graph(C.subgraph(comps[0]))
    if L.number_of_nodes() < 3:
        return None
    mapping = {u: k for k, u in enumerate(sorted(L.nodes()))}
    L2 = nx.Graph()
    L2.add_nodes_from(range(len(mapping)))
    L2.add_edges_from((mapping[a], mapping[b]) for a, b in L.edges())
    return L2


def relabel_sorted(G):
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    H.add_edges_from(G.edges())
    return nx.convert_node_labels_to_integers(H, ordering="sorted")


def powerlaw_config(n=300, gamma=2.5, dmin=2, seed=1):
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    deg = np.floor(dmin * (1.0 - u) ** (-1.0 / (gamma - 1.0))).astype(int)
    deg = np.minimum(deg, n - 1)
    if deg.sum() % 2 == 1:
        deg[int(np.argmax(deg))] -= 1 if deg.max() > dmin else -1
    M = nx.configuration_model([int(d) for d in deg], seed=seed)
    G = nx.Graph(M)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    return G


def named_graphs():
    return [
        ("karate", relabel_sorted(nx.karate_club_graph()),
         "nx.karate_club_graph, unweighted"),
        ("les_miserables", relabel_sorted(nx.les_miserables_graph()),
         "nx.les_miserables_graph, unweighted, labels sorted"),
        ("florentine", relabel_sorted(nx.florentine_families_graph()),
         "nx.florentine_families_graph, labels sorted"),
        ("er_200", nx.gnp_random_graph(200, 4 / 200, seed=1),
         "nx.gnp_random_graph(200, 4/200, seed=1)"),
        ("ba_300", nx.barabasi_albert_graph(300, 2, seed=1),
         "nx.barabasi_albert_graph(300, 2, seed=1)"),
        ("holme_kim_300", nx.powerlaw_cluster_graph(300, 2, 0.3, seed=1),
         "nx.powerlaw_cluster_graph(300, 2, 0.3, seed=1)"),
        ("regular3_100", nx.random_regular_graph(3, 100, seed=1),
         "nx.random_regular_graph(3, 100, seed=1) (control)"),
        ("config_powerlaw_300", powerlaw_config(),
         "degrees floor(2*(1-U)^(-1/1.5)), U~numpy default_rng(1).random(300),"
         " capped at n-1, parity fixed on max degree; "
         "nx.configuration_model(seed=1), simplified (multi-edges merged,"
         " self-loops removed)"),
    ]


# ---------------------------------------------------------------- analysis
def scores_on_grid(m, eps):
    rows = [m.pi(mu, eps) for mu in MU_GRID]
    pinf = m.pi_inf(eps)
    S = np.vstack(rows + [pinf])  # (len(grid)+1, n), last row = mu->inf
    i1 = int(np.where(MU_GRID == 1.0)[0][0])
    check(np.max(np.abs(S[i1] - m.pi_pagerank(eps))) < 1e-12,
          "edge solver pi_1 != PageRank formula")
    return S, i1


def nonmonotone_mask(S, tol):
    d = np.diff(S, axis=0)
    return np.any(d > tol, axis=0) & np.any(d < -tol, axis=0)


def sgn(x, tol):
    return np.where(x > tol, 1, np.where(x < -tol, -1, 0))


def analyse(m, eps, full=True):
    S, i1 = scores_on_grid(m, eps)
    n = m.n
    p0, p1, pinf = S[0], S[i1], S[-1]
    tol = TOL_REL * float(S.max())
    nm = nonmonotone_mask(S, tol)
    n01 = int(np.sum(MU_GRID <= 1.0))
    nm01 = nonmonotone_mask(S[:n01], tol)
    out = {
        "n": n, "m": m.G.number_of_edges(),
        "frac_nonmonotone": r6(nm.mean()),
        "count_nonmonotone": int(nm.sum()),
        "frac_nonmonotone_mu_0_1": r6(nm01.mean()),
        "count_nonmonotone_mu_0_1": int(nm01.sum()),
        "tol_abs": r6(tol),
    }
    if not full:
        return out
    hi = S.max(axis=0)
    lo = S.min(axis=0)
    exc = (hi - np.maximum(p0, pinf)) + (np.minimum(p0, pinf) - lo)
    rel = exc / p1
    if nm.any():
        out["rel_excursion_median"] = r6(np.median(rel[nm]))
        out["rel_excursion_max"] = r6(np.max(rel[nm]))
        out["median_degree_nonmonotone"] = r6(np.median(m.deg[nm]))
    else:
        out["rel_excursion_median"] = None
        out["rel_excursion_max"] = None
        out["median_degree_nonmonotone"] = None
    out["median_degree_all"] = r6(np.median(m.deg))
    s0 = sgn(p0 - p1, tol)
    si = sgn(pinf - p1, tol)
    elig = (s0 != 0) & (si != 0)
    wrong = elig & (s0 != -si)
    out["shortcut_wrong_fraction"] = r6(wrong.mean())
    out["shortcut_wrong_count"] = int(wrong.sum())
    out["shortcut_eligible_count"] = int(elig.sum())
    # Rank statistics are undefined (pure float noise) when a score vector
    # is constant within tol, e.g. on a regular graph: report null.
    flat = [float(np.ptp(x)) <= tol for x in (p0, p1, pinf)]
    t01 = None if (flat[0] or flat[1]) else kendalltau(p0, p1).statistic
    t1i = None if (flat[1] or flat[2]) else kendalltau(p1, pinf).statistic
    out["kendall_tau_pi0_pi1"] = r6(t01)
    out["kendall_tau_pi1_piinf"] = r6(t1i)
    if flat[0] or flat[1]:
        out["top10_overlap_pi0_pi1"] = None
    else:
        top0 = set(np.argsort(-p0, kind="stable")[:10].tolist())
        top1 = set(np.argsort(-p1, kind="stable")[:10].tolist())
        out["top10_overlap_pi0_pi1"] = len(top0 & top1)
    # pairwise quantities, all pairs i<j
    iu, ju = np.triu_indices(n, k=1)
    npairs = len(iu)
    rev = 0
    cert = np.zeros(npairs, dtype=bool)
    r = m.P @ (m.P @ m.v) - m.v
    bound = eps ** 2 / (1.0 - eps ** 2) * float(np.abs(r).sum())
    cert_violate = 0
    chunk = 20000
    for a in range(0, npairs, chunk):
        I = iu[a:a + chunk]
        J = ju[a:a + chunk]
        D = S[:, I] - S[:, J]
        sg = sgn(D, tol)
        same = (sg[0] == sg[-1]) & (sg[0] != 0)
        opp = np.any(sg == -sg[0][None, :], axis=0)
        rev += int(np.sum(same & opp))
        c = np.abs(D[-1]) > bound
        cert[a:a + chunk] = c
        flip = np.any(sg != sg[-1][None, :], axis=0)
        cert_violate += int(np.sum(c & flip))
    out["transient_reversal_pairs"] = rev
    out["transient_reversal_fraction"] = r6(rev / npairs)
    out["pairs_total"] = npairs
    out["pairs_restriction"] = "none (all vertex pairs)"
    out["residual_l1_bound"] = r6(bound)
    out["residual_certified_pair_fraction"] = r6(cert.mean())
    out["residual_certified_pairs_with_grid_order_change"] = cert_violate
    return out


def er_ensemble():
    res = {}
    for n in (20, 100):
        per_eps = {}
        graphs = []
        skipped = 0
        for seed in range(50):
            L = core_lcc(nx.gnp_random_graph(n, 4 / n, seed=seed))
            if L is None:
                skipped += 1
                continue
            graphs.append(EdgeModel(L))
        for eps in EPS_LIST:
            fr, fr01, sizes, anyfrac = [], [], [], []
            for m in graphs:
                a = analyse(m, eps, full=False)
                fr.append(a["count_nonmonotone"] / a["n"])
                fr01.append(a["count_nonmonotone_mu_0_1"] / a["n"])
                sizes.append(a["n"])
                anyfrac.append(a["count_nonmonotone"] > 0)
            fr = np.array(fr)
            fr01 = np.array(fr01)
            per_eps[f"eps={eps}"] = {
                "graphs_used": len(fr),
                "frac_nonmonotone_mean": r6(fr.mean()),
                "frac_nonmonotone_median": r6(np.median(fr)),
                "frac_nonmonotone_max": r6(fr.max()),
                "frac_nonmonotone_mu_0_1_mean": r6(fr01.mean()),
                "frac_nonmonotone_mu_0_1_median": r6(np.median(fr01)),
                "fraction_graphs_with_any_nonmonotone": r6(np.mean(anyfrac)),
                "core_lcc_n_mean": r6(np.mean(sizes)),
                "core_lcc_n_median": r6(np.median(sizes)),
            }
        res[f"n={n}"] = {
            "generator": f"nx.gnp_random_graph({n}, 4/{n}, seed=s), s=0..49",
            "skipped_empty_or_tiny_core": skipped,
            "results": per_eps,
        }
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir", help="output directory (must exist or be "
                    "creatable); writes prevalence.json there")
    args = ap.parse_args()
    outdir = os.path.abspath(args.outdir)
    here = os.path.dirname(os.path.abspath(__file__))
    check(outdir != here, "refusing to write beside the script")
    os.makedirs(outdir, exist_ok=True)

    result = {
        "status": "numerical",
        "description": "Prevalence of non-monotone vertex scores pi_mu in "
                       "backtracking-weighted PageRank (tail projection), "
                       "uniform restart v_i=1/n, graphs = LCC of 2-core.",
        "definitions": {
            "kernel": "(i,j)->(j,k), k~j, weight (1+(mu-1)[k==i])/(d_j-1+mu)",
            "restart": "u_(i,j)=v_i/d_i, v uniform",
            "score": "pi_mu(i)=sum_{j~i} q(i,j), q=(1-eps)u+eps Q_mu q",
            "pi_inf": "(v+eps P v)/(1+eps), P=A D^-1",
            "mu_grid": "np.unique(concat(linspace(0,1,101), "
                       "geomspace(1,1000,200))) then mu->inf (pi_inf)",
            "mu_grid_size_finite": int(len(MU_GRID)),
            "tolerance": "tol = 1e-12 * max over grid and vertices of pi",
            "nonmonotone": "consecutive differences along grid (incl. pi_inf) "
                           "contain one > tol and one < -tol",
            "nonmonotone_mu_0_1": "same, restricted to grid points mu<=1",
            "excursion": "(max_grid pi - max(pi_0,pi_inf)) + "
                         "(min(pi_0,pi_inf) - min_grid pi); relative = "
                         "excursion/pi_1(i); stats over nonmonotone vertices",
            "shortcut_wrong": "vertices with sign(pi_0-pi_1) != "
                              "-sign(pi_inf-pi_1), both nonzero beyond tol; "
                              "fraction over all n vertices",
            "transient_reversal": "pairs with sign(pi_0 diff)==sign(pi_inf "
                                  "diff)!=0 and some grid point with the "
                                  "opposite sign (margins > tol)",
            "residual_bound": "r=P^2 v - v; bound=eps^2/(1-eps^2)*||r||_1; "
                              "certified pair: |pi_inf(i)-pi_inf(j)|>bound",
            "residual_certified_pairs_with_grid_order_change":
                "certified pairs whose sign (beyond tol) at some grid point "
                "differs from the pi_inf sign (consistency check)",
            "rounding": "floats rounded to 6 significant figures",
        },
        "eps_values": list(EPS_LIST),
        "software": {"numpy": np.__version__, "scipy": scipy.__version__,
                     "networkx": nx.__version__},
        "validation_house": validate_house(),
        "graphs": {},
    }
    for name, G, desc in named_graphs():
        L = core_lcc(G)
        check(L is not None, f"{name}: empty 2-core")
        m = EdgeModel(L)
        entry = {"source": desc,
                 "original_n": G.number_of_nodes(),
                 "original_m": G.number_of_edges(),
                 "by_eps": {}}
        for eps in EPS_LIST:
            entry["by_eps"][f"eps={eps}"] = analyse(m, eps)
        result["graphs"][name] = entry
    reg = result["graphs"]["regular3_100"]["by_eps"]
    for k in reg:
        check(reg[k]["count_nonmonotone"] == 0,
              "3-regular control is nonmonotone")
    result["er_ensemble"] = er_ensemble()
    path = os.path.join(outdir, "prevalence.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, sort_keys=True, indent=1)
        fh.write("\n")


if __name__ == "__main__":
    main()
