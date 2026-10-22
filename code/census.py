# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Certified census of coordinatewise monotonicity on small graphs.

For every connected simple graph on 2..7 vertices (networkx graph atlas,
which lists every graph on at most seven vertices up to isomorphism) and
every damping value in DAMPING, decide exactly, for each vertex, whether
mu -> pi_mu(i) is monotone on [0, inf) and on [0, 1].

Certification: pi_mu(i) = N_i(mu)/Den(mu) with integer polynomials obtained
by exact interpolation (degree bound 2m proved; two extra exact checks).
Den > 0 on [0, inf).  pi is non-monotone on an open interval iff
W = N' Den - N Den' has a root of odd multiplicity there; odd-multiplicity
roots are counted exactly on each open interval (squarefree decomposition,
exact removal of roots at 0 and 1, Sturm sequence).
"""
from fractions import Fraction as F
import networkx as nx
import bt_exact as B

DAMPING = [F(1, 10), F(1, 4), F(1, 3), F(1, 2), F(17, 20)]
MAX_ORDER = 7


def atlas_graphs():
    out = []
    for idx, G in enumerate(nx.graph_atlas_g()):
        n = G.number_of_nodes()
        if n < 2 or n > MAX_ORDER or not nx.is_connected(G):
            continue
        out.append((idx, n, sorted(tuple(sorted(e)) for e in G.edges())))
    return out


def analyse_graph(item):
    idx, n, edges = item
    degs = [0] * n
    for a, b in edges:
        degs[a] += 1
        degs[b] += 1
    rec = {"atlas": idx, "n": n, "m": len(edges), "edges": [list(e) for e in edges],
           "min_degree": min(degs), "per_eps": {}}
    for eps in DAMPING:
        N, D = B.rational_response(n, edges, eps)
        bad, bad01 = [], []
        for i in range(n):
            r = B.sign_change_analysis(N[i], D)
            if r["roots_pos"] > 0:
                bad.append(i)
            if r["roots_01"] > 0:
                bad01.append(i)
        rec["per_eps"][str(eps)] = {"nonmonotone": bad, "nonmonotone_on_01": bad01}
    return rec


def run(jobs=1):
    items = atlas_graphs()
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(jobs) as pool:
            recs = pool.map(analyse_graph, items, chunksize=4)
    else:
        recs = [analyse_graph(it) for it in items]
    recs.sort(key=lambda r: r["atlas"])
    summary = {}
    for n in range(2, MAX_ORDER + 1):
        for cls in ("all", "mindeg2"):
            sel = [r for r in recs if r["n"] == n and (cls == "all" or r["min_degree"] >= 2)]
            row = {"graphs": len(sel)}
            for eps in DAMPING:
                k = str(eps)
                row["nonmonotone_" + k] = sum(1 for r in sel if r["per_eps"][k]["nonmonotone"])
                row["nonmonotone_on_01_" + k] = sum(1 for r in sel if r["per_eps"][k]["nonmonotone_on_01"])
            row["nonmonotone_some_eps"] = sum(
                1 for r in sel if any(r["per_eps"][str(e)]["nonmonotone"] for e in DAMPING))
            row["nonmonotone_all_eps"] = sum(
                1 for r in sel if all(r["per_eps"][str(e)]["nonmonotone"] for e in DAMPING))
            vtx = sum(r["n"] for r in sel)
            row["vertices"] = vtx
            for eps in DAMPING:
                k = str(eps)
                row["nonmonotone_vertices_" + k] = sum(len(r["per_eps"][k]["nonmonotone"]) for r in sel)
            summary["n%d_%s" % (n, cls)] = row
    return {"damping": [str(e) for e in DAMPING], "graphs": recs, "summary": summary,
            "method": "exact: interpolation of N_i, Den (degree <= 2m), squarefree "
                      "decomposition of W = N'Den - N Den', odd part, roots at 0 and 1 "
                      "divided out exactly, Sturm counts on (0, inf) and (0, 1)"}
