# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Exact checks for the endpoint-derivative formula, the endpoint tests,
the small-damping two-step criterion, the spread clauses, a transient rank
reversal, a finite exact rank-reversal search, the tree and unique-extrema
counterexamples, root-counter regressions, and the degree-tail constructions."""
from fractions import Fraction as F
import flint
import bt_exact as B

check = B.check


# ------------------------------------------------------------ endpoint formula
def endpoint_derivative_formula(n, edges, eps):
    """d pi / d mu at mu = 1 (min degree >= 2):
       eps^2 (I - eps P)^{-1} [ h o D^{-1} pi_1 - P D^{-1} P pi_1 ]."""
    eps = F(eps)
    nb = B.adjacency(n, edges)
    d = [len(x) for x in nb]
    p1 = B.pi_one(n, edges, eps)
    h = [sum(F(1, d[j]) for j in nb[k]) for k in range(n)]
    Pp = [sum(p1[j] / d[j] for j in nb[k]) for k in range(n)]
    DPp = [Pp[k] / d[k] for k in range(n)]
    PDPp = [sum(DPp[j] / d[j] for j in nb[k]) for k in range(n)]
    g = [h[k] * p1[k] / d[k] - PDPp[k] for k in range(n)]
    # solve (I - eps P) y = g
    A = [[F(int(r == c)) - (eps * F(1, d[c]) if c in nb[r] else 0) for c in range(n)] for r in range(n)]
    y, _ = B.solve_det(A, g)
    return [eps * eps * yk for yk in y], g


def deriv_at(N, D, x):
    Np, Dp = N.derivative(), D.derivative()
    Dx = B.peval(D, x)
    return (B.peval(Np, x) * Dx - B.peval(N, x) * B.peval(Dp, x)) / (Dx * Dx)


def two_step_profile(n, edges, k):
    """f_k(mu) = n * M_k(mu) as (numerator, denominator) fmpq data, uniform restart."""
    nb = B.adjacency(n, edges)
    d = [len(x) for x in nb]
    h = [sum(F(1, d[j]) for j in nb[i]) for i in range(n)]
    terms = []
    for j in nb[k]:
        # (h_j - 1/d_k + mu/d_k) / (d_j - 1 + mu)
        terms.append(((h[j] - F(1, d[k]), F(1, d[k])), (F(d[j] - 1), F(1))))
    sigma = [F(d[j], d[k]) - h[j] for j in nb[k]]
    return terms, sigma


def profile_polys(terms):
    num = flint.fmpq_poly([0])
    den = flint.fmpq_poly([1])
    for (a0, a1), (b0, b1) in terms:
        q = lambda x: flint.fmpq(x.numerator, x.denominator)
        tn = flint.fmpq_poly([q(a0), q(a1)])
        td = flint.fmpq_poly([q(b0), q(b1)])
        num = num * td + tn * den
        den = den * td
    return num, den


def fmpq_to_fmpz(p):
    c = [F(int(x.p), int(x.q)) for x in p.coeffs()]
    return B.to_fmpz(c)


def profile_nonmonotone(terms):
    num, den = profile_polys(terms)
    r = B.sign_change_analysis(fmpq_to_fmpz(num), fmpq_to_fmpz(den))
    return r["roots_pos"] > 0


def run_endpoint_and_twostep(census_result):
    out = {}
    recs = [r for r in census_result["graphs"] if r["min_degree"] >= 2]
    eps = F(1, 2)
    nverified = 0
    detect = {"total_nonmonotone": 0, "test_same_side": 0, "test_derivative_vs_infinity": 0,
              "test_derivative_vs_zero": 0, "any_test": 0}
    for r in recs:
        n, edges = r["n"], [tuple(e) for e in r["edges"]]
        formula, _ = endpoint_derivative_formula(n, edges, eps)
        N, D = B.rational_response(n, edges, eps)
        p0 = B.solve_point(n, edges, eps, 0)[0]
        p1 = B.solve_point(n, edges, eps, 1)[0]
        pinf = B.pi_inf(n, edges, eps)
        bad = set(r["per_eps"][str(eps)]["nonmonotone"])
        for i in range(n):
            dv = deriv_at(N[i], D, 1)
            check(dv == formula[i], "endpoint derivative formula mismatch")
            nverified += 1
            t1 = (p0[i] - p1[i]) * (pinf[i] - p1[i]) > 0
            t2 = dv * (pinf[i] - p1[i]) < 0
            t3 = dv * (p1[i] - p0[i]) < 0
            if t1 or t2 or t3:
                check(i in bad, "endpoint test flagged a monotone vertex")
            if i in bad:
                detect["total_nonmonotone"] += 1
                detect["test_same_side"] += int(t1)
                detect["test_derivative_vs_infinity"] += int(t2)
                detect["test_derivative_vs_zero"] += int(t3)
                detect["any_test"] += int(t1 or t2 or t3)
    out["endpoint_derivative_formula_vertices_checked"] = nverified
    out["endpoint_tests_eps_1_2"] = detect
    # two-step criterion at eps = 1/10 versus census
    tw = {"vertices": 0, "profile_nonmonotone": 0, "profile_nonmonotone_and_census_nonmonotone": 0,
          "census_nonmonotone": 0, "mixed_sigma": 0, "profile_nonmonotone_with_uniform_sigma": 0}
    for r in recs:
        n, edges = r["n"], [tuple(e) for e in r["edges"]]
        bad = set(r["per_eps"]["1/10"]["nonmonotone"])
        for k in range(n):
            terms, sigma = two_step_profile(n, edges, k)
            nm = profile_nonmonotone(terms)
            mixed = any(s > 0 for s in sigma) and any(s < 0 for s in sigma)
            tw["vertices"] += 1
            tw["profile_nonmonotone"] += int(nm)
            tw["census_nonmonotone"] += int(k in bad)
            tw["profile_nonmonotone_and_census_nonmonotone"] += int(nm and k in bad)
            tw["mixed_sigma"] += int(mixed)
            tw["profile_nonmonotone_with_uniform_sigma"] += int(nm and not mixed)
    check(tw["profile_nonmonotone_with_uniform_sigma"] == 0, "sigma criterion violated")
    out["two_step_vs_census_eps_1_10"] = tw
    # explicit profiles: house roof and gem vertex 1
    house = [(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3)]
    _, sr = two_step_profile(5, house, 4)
    gem = [(0, 1), (1, 2), (2, 3), (0, 4), (1, 4), (2, 4), (3, 4)]
    tg, sg = two_step_profile(5, gem, 1)
    check(all(s > 0 for s in sr), "house roof sigma")
    check(profile_nonmonotone(tg), "gem profile")
    out["sigma_house_roof"] = [str(s) for s in sr]
    out["sigma_gem_vertex1"] = [str(s) for s in sg]
    return out


# ------------------------------------------------------------ spread clauses
def spread(p):
    return max(p) - min(p)


def run_spread(census_result):
    out = {}
    house = [(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3)]
    mus = [F(0), F(1, 10), F(1, 2), F(1), F(2), F(5), F(20)]
    sp = [spread(B.solve_point(5, house, F(1, 2), m)[0]) for m in mus]
    spi = spread(B.pi_inf(5, house, F(1, 2)))
    check(all(a > b for a, b in zip(sp, sp[1:] + [spi])), "house spread not strictly decreasing")
    out["house_eps_1_2"] = {"mu": [str(m) for m in mus] + ["inf"], "spread": [str(s) for s in sp + [spi]]}
    # reversed clause: spread nonincreasing -- counterexample (atlas index 129)
    rec = next(r for r in census_result["graphs"] if r["atlas"] == 129)
    n, edges = rec["n"], [tuple(e) for e in rec["edges"]]
    eps = F(17, 20)
    s0 = spread(B.solve_point(n, edges, eps, 0)[0])
    s1 = spread(B.solve_point(n, edges, eps, F(1, 10))[0])
    check(s1 > s0, "reversed spread counterexample")
    check(rec["min_degree"] >= 2, "reversed spread graph min degree")
    out["reversed_clause_counterexample"] = {"atlas": 129, "n": n, "edges": [list(e) for e in edges],
                                              "eps": str(eps), "spread_mu_0": str(s0),
                                              "spread_mu_1_10": str(s1)}
    return out


# ------------------------------------------------------------ rank reversal
def les_mis_two_core():
    import networkx as nx
    G = nx.Graph(nx.les_miserables_graph().edges())
    C = nx.k_core(G, 2)
    C = C.subgraph(max(nx.connected_components(C), key=len))
    names = sorted(C.nodes())
    ix = {v: t for t, v in enumerate(names)}
    edges = sorted(tuple(sorted((ix[a], ix[b]))) for a, b in C.edges())
    return names, edges


def run_rank_reversal():
    names, edges = les_mis_two_core()
    n = len(names)
    eps = F(1, 2)
    i, j = names.index("Fauchelevent"), names.index("Pontmercy")
    p0 = B.solve_point(n, edges, eps, 0)[0]
    p1 = B.solve_point(n, edges, eps, 1)[0]
    pinf = B.pi_inf(n, edges, eps)
    d0, d1, di = p0[i] - p0[j], p1[i] - p1[j], pinf[i] - pinf[j]
    check(d0 < 0 and d1 > 0 and di < 0, "rank reversal certificate")
    return {"graph": "Les Miserables co-occurrence graph, unweighted, largest component of 2-core",
            "n": n, "m": len(edges), "eps": str(eps), "pair": ["Fauchelevent", "Pontmercy"],
            "diff_mu_0": str(d0), "diff_mu_1": str(d1), "diff_mu_inf": str(di)}


# ------------------------------------------------------------ degree tails
def clique_chain(L, beta):
    edges, portals, sizes, blocks = [], [], [], []
    n = 0
    for m in range(3, L + 1):
        for _ in range(L ** beta // m ** beta):
            portals.append(n)
            sizes.append(m)
            blocks.append(list(range(n, n + m)))
            edges.extend((a, b) for a in range(n, n + m) for b in range(a + 1, n + m))
            n += m
    edges.extend(zip(portals, portals[1:]))
    return n, edges, portals, sizes, blocks


def run_tails():
    out = []
    for beta, Ls in ((4, (4, 5, 8, 12, 16)), (5, (4, 6, 9))):
        for L in Ls:
            n, edges, portals, sizes, blocks = clique_chain(L, beta)
            nb = B.adjacency(n, edges)
            d = [len(x) for x in nb]
            h = [sum(F(1, d[j]) for j in nb[i]) for i in range(n)]
            K = len(sizes)
            # discarded cliques (1-based): 1, 2, K-1, K and both cliques at each size change
            disc = {1, 2, K - 1, K}
            for t in range(1, K):
                if sizes[t - 1] != sizes[t]:
                    disc.update({t, t + 1})
            surv = 0
            for t in range(1, K + 1):
                if t in disc:
                    continue
                m = sizes[t - 1]
                for vtx in blocks[t - 1]:
                    for eps in (F(1, 2), F(17, 20)):
                        score = (1 + eps * h[vtx]) / (1 + eps)
                        if vtx == portals[t - 1]:
                            want = 1 + 2 * eps / ((1 + eps) * (m + 1))
                        else:
                            want = 1 - 2 * eps / ((1 + eps) * (m * m - 1))
                        check(score == want, "generic clique score formula")
                surv += m
            check(max(h) <= F(5, 3), "reciprocal degree bound")
            excluded = n - surv
            check(excluded <= (4 + 2 * L) * L, "excluded vertex count bound")
            # two-sided local bound under min degree two
            for eps in (F(1, 2),):
                for i in range(n):
                    check((1 + eps * h[i]) / (1 + eps) <= (1 + eps * F(d[i], 2)) / (1 + eps), "degree bound")
            out.append({"beta": beta, "L": L, "n": n, "m": len(edges), "cliques": K,
                        "surviving_vertices": surv, "excluded_vertices": excluded})
    # regression: boundary clique (Theorem 11.2 proof); L = 4 chain K3-K3-K3-K4, vertex 3
    n, edges, portals, sizes, blocks = clique_chain(4, 4)
    check(sizes == [3, 3, 3, 4], "L=4 chain sizes")
    nb = B.adjacency(n, edges)
    d = [len(x) for x in nb]
    h3 = sum(F(1, d[j]) for j in nb[3])
    eps = F(1, 2)
    check((1 + eps * h3) / (1 + eps) == F(43, 36), "regression: endpoint-neighbour clique")
    return {"chains": out, "regression_L4_vertex3_score_eps_1_2": "43/36"}


def run_degree_bound(census_result):
    cnt = 0
    for r in census_result["graphs"]:
        if r["min_degree"] < 2:
            continue
        n, edges = r["n"], [tuple(e) for e in r["edges"]]
        nb = B.adjacency(n, edges)
        d = [len(x) for x in nb]
        for eps in (F(1, 2),):
            pinf = B.pi_inf(n, edges, eps)
            for i in range(n):
                check(n * pinf[i] <= (1 + eps * F(d[i], 2)) / (1 + eps), "n pi_inf degree bound")
                cnt += 1
    return cnt


# ------------------------------------------------------------ root-counter regressions
def run_root_regressions(census_result):
    out = B.known_root_regressions()
    # Atlas 399, eps = 1/2, coordinate 4: nonmonotone inside (0, 1)
    edges = [(0, 1), (0, 3), (0, 4), (1, 3), (1, 4), (1, 6), (2, 6), (4, 5)]
    rec = next(r for r in census_result["graphs"] if r["atlas"] == 399)
    check([tuple(e) for e in rec["edges"]] == edges, "atlas 399 edges")
    N, D = B.rational_response(7, edges, F(1, 2))
    r = B.sign_change_analysis(N[4], D)
    check(r["roots_pos"] == 2 and r["roots_01"] == 1, "atlas 399 coordinate 4 root counts")
    d0, d1 = deriv_at(N[4], D, 0), deriv_at(N[4], D, 1)
    check(d0 < 0 and d1 > 0, "atlas 399 coordinate 4 derivative signs")
    check(4 in rec["per_eps"]["1/2"]["nonmonotone_on_01"], "atlas 399 census entry")
    out["atlas_399_eps_1_2_vertex_4"] = {"roots_pos": r["roots_pos"], "roots_01": r["roots_01"],
                                         "derivative_mu_0": str(d0), "derivative_mu_1": str(d1)}
    return out


# ------------------------------------------------------------ tree counterexample
TREE_EDGES = [(0, 4), (1, 3), (2, 3), (3, 4)]


def run_tree():
    import sympy as sp
    eps = F(1, 2)
    n = 5
    N, D = B.rational_response(n, TREE_EDGES, eps)
    # claimed closed form pi_mu(1) = (21 mu^2 + 95 mu + 84) / (15 (9 mu^2 + 42 mu + 37))
    num = flint.fmpz_poly([84, 95, 21])
    den = flint.fmpz_poly([555, 630, 135])
    check(N[1] * den == D * num, "tree leaf closed form")
    x = sp.symbols("x")
    f = (21 * x ** 2 + 95 * x + 84) / (15 * (9 * x ** 2 + 42 * x + 37))
    fd = sp.cancel(sp.diff(f, x) - (27 * x ** 2 + 42 * x - 13) / (15 * (9 * x ** 2 + 42 * x + 37) ** 2))
    check(fd == 0, "tree leaf derivative")
    root = (2 * sp.sqrt(22) - 7) / 9
    check(sp.expand(27 * root ** 2 + 42 * root - 13) == 0 and root > 0 and root < 1, "tree turning point")
    r = B.sign_change_analysis(N[1], D)
    check(r["roots_pos"] == 1 and r["roots_01"] == 1, "tree root counts")
    mus = [F(1, 10), F(1, 4), F(1, 2)]
    vals = [B.solve_full(n, TREE_EDGES, eps, m)[1] for m in mus]
    check(vals == [F(9371, 61935), F(349, 2307), F(547, 3615)], "tree three-point values")
    check(vals[1] < vals[0] and vals[1] < vals[2], "tree three-point certificate")
    return {"edges": [list(e) for e in TREE_EDGES], "eps": "1/2", "vertex": 1,
            "mu": [str(m) for m in mus], "values": [str(v) for v in vals],
            "differences": [str(vals[1] - vals[0]), str(vals[2] - vals[1])],
            "turning_point": "(2*sqrt(22)-7)/9"}


# ------------------------------------------------------------ unique-extrema spread counterexample
UNIQ_EDGES = [(0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (5, 6)]


def _positive_coeffs(p):
    return all(int(c) > 0 for c in p.coeffs())


def run_unique_spread():
    n, eps = 7, F(17, 20)
    N, D = B.rational_response(n, UNIQ_EDGES, eps)
    pinf = B.pi_inf(n, UNIQ_EDGES, eps)
    gaps = 0
    for top, j in [(2, j) for j in range(n) if j != 2] + [(j, 0) for j in range(n) if j != 0]:
        if top == j:
            continue
        g = N[top] - N[j]
        c = g.gcd(D)
        gn, r1 = divmod(g, c)
        gd, r2 = divmod(D, c)
        check(r1.is_zero() and r2.is_zero(), "gap reduction")
        if int(gd.coeffs()[-1]) < 0:
            gn, gd = -gn, -gd
        check(_positive_coeffs(gn) and _positive_coeffs(gd), "unique extrema: positive-coefficient gap")
        check(gn.degree() == gd.degree(), "gap degrees")
        check(pinf[top] > pinf[j], "unique extrema at infinity")
        gaps += 1
    check(gaps == 12, "twelve gaps")
    s0 = B.solve_full(n, UNIQ_EDGES, eps, 0)
    s1 = B.solve_full(n, UNIQ_EDGES, eps, F(1, 10))
    d0, d1 = s0[2] - s0[0], s1[2] - s1[0]
    check(d0 == spread(s0) and d1 == spread(s1), "spread uses unique extrema")
    check(d1 > d0, "spread increases although extrema are unique")
    return {"edges": [list(e) for e in UNIQ_EDGES], "eps": str(eps), "max_vertex": 2, "min_vertex": 0,
            "gaps_with_positive_coefficients": gaps, "spread_mu_0": str(d0), "spread_mu_1_10": str(d1),
            "increment": str(d1 - d0)}


# ------------------------------------------------------------ finite exact rank-reversal search
RANK_GRID = [F(1, 4), F(1, 2), F(1), F(2), F(4), F(8), F(16)]


def run_rank_search(census_result):
    """Exact search: pairs ordered strictly the same way at mu = 0 and mu = inf whose
    order is strictly reversed at some grid point.  Finite grid only."""
    recs = [r for r in census_result["graphs"] if r["n"] == 7 and r["min_degree"] >= 2]
    cases = comparisons = reversals = ties = 0
    for r in recs:
        n, edges = r["n"], [tuple(e) for e in r["edges"]]
        for eps in (F(1, 2), F(17, 20)):
            cases += 1
            p0 = B.solve_point(n, edges, eps, 0)[0]
            pinf = B.pi_inf(n, edges, eps)
            grid = [B.solve_point(n, edges, eps, m)[0] for m in RANK_GRID]
            for i in range(n):
                for j in range(i + 1, n):
                    a, b = p0[i] - p0[j], pinf[i] - pinf[j]
                    if a * b <= 0:
                        continue
                    for pm in grid:
                        c = pm[i] - pm[j]
                        comparisons += 1
                        reversals += int(c * a < 0)
                        ties += int(c == 0)
    check(len(recs) == 507, "507 graphs")
    check(reversals == 0, "rank search found a reversal")
    return {"graphs": len(recs), "graph_damping_cases": cases, "eps": ["1/2", "17/20"],
            "mu_grid": [str(m) for m in RANK_GRID], "eligible_pair_grid_comparisons": comparisons,
            "strict_reversals": reversals, "ties": ties}
