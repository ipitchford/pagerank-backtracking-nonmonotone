# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Exact bivariate (mu, epsilon) certificates for P3, the paw and the gem.

Scores are obtained by solving the full directed-edge system over the field
Q(mu, epsilon) (sympy DomainMatrix LU, exact).  Sign claims on intervals of
epsilon are certified by exact Sturm root counts on open intervals
(bt_exact.roots_open) and exact rational evaluation.
"""
from fractions import Fraction as F
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import flint
import bt_exact as B

check = B.check
mu, e = sp.symbols("mu epsilon")
K = sp.QQ.frac_field(mu, e)

P3_EDGES = [(0, 1), (1, 2)]
PAW_EDGES = [(0, 3), (1, 2), (1, 3), (2, 3)]          # triangle 1,2,3; pendant 0 at 3
GEM_EDGES = [(0, 1), (1, 2), (2, 3), (0, 4), (1, 4), (2, 4), (3, 4)]  # path 0-1-2-3, apex 4
DIAMOND_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]  # hubs 0,1; tips 2,3


def symbolic_scores(n, edges):
    nb = B.adjacency(n, edges)
    d = [len(x) for x in nb]
    arcs = [(i, j) for i in range(n) for j in nb[i]]
    idx = {a: t for t, a in enumerate(arcs)}
    N = len(arcs)
    E, MU, one, zero = K.from_sympy(e), K.from_sympy(mu), K.one, K.zero
    M = [[one if r == c else zero for c in range(N)] for r in range(N)]
    b = [[zero] for _ in range(N)]
    for (i, j) in arcs:
        b[idx[(i, j)]][0] = (one - E) * K.from_sympy(sp.Rational(1, n * d[i]))
        for k in nb[j]:
            if d[j] == 1:
                w = one
            else:
                w = (one + (MU - one) * (one if k == i else zero)) / (K.from_sympy(sp.Integer(d[j] - 1)) + MU)
            M[idx[(j, k)]][idx[(i, j)]] -= E * w
    q = DomainMatrix(M, (N, N), K).lu_solve(DomainMatrix(b, (N, 1), K)).to_Matrix()
    return [sp.cancel(sum(q[idx[(i, j)]] for j in nb[i])) for i in range(n)]


def upoly_e(expr):
    """sympy expression in epsilon -> integer fmpz_poly (positive multiple)."""
    P = sp.Poly(sp.together(expr), e)
    coeffs = [sp.Rational(c) for c in reversed(P.all_coeffs())]
    den = 1
    for c in coeffs:
        den = sp.ilcm(den, c.q)
    return flint.fmpz_poly([int(c * den) for c in coeffs])


def roots_in_open_interval(p, a, b):
    """Number of distinct real roots of fmpz_poly p in (a, b), a < b rational (exact Sturm count)."""
    return B.roots_open(p, F(a), F(b))


def peval_e(p, x):
    return B.peval(p, x)


def bracket(p, lo, hi):
    """Exact: p has opposite strict signs at lo, hi and exactly one root in (lo, hi)."""
    vlo, vhi = peval_e(p, lo), peval_e(p, hi)
    check(vlo * vhi < 0, "bracket signs")
    check(roots_in_open_interval(p, lo, hi) == 1, "bracket uniqueness")
    return [str(F(lo)), str(F(hi))]


def num_deriv(s):
    return sp.factor(sp.numer(sp.together(sp.diff(s, mu))))


def run():
    out = {}
    # ---------------- P3: every coordinate is independent of mu
    sP3 = symbolic_scores(3, P3_EDGES)
    check(all(sp.simplify(sp.diff(s, mu)) == 0 for s in sP3), "P3 not constant")
    out["P3_scores"] = [str(sp.factor(s)) for s in sP3]

    # ---------------- paw
    sp_paw = symbolic_scores(4, PAW_EDGES)
    paw = {}
    for v, s in enumerate(sp_paw):
        dn = sp.expand(sp.numer(sp.together(sp.diff(s, mu))))
        dd = sp.denom(sp.together(sp.diff(s, mu)))
        paw[v] = (dn, dd)
    # derivative numerators (up to the positive factor e^2 (1-e)/den^2); signs:
    # d/dmu pi = dn/dd with dd a square times a positive constant -> check sign of dd
    res = {}
    for v in range(4):
        dn, dd = paw[v]
        Q = sp.Poly(sp.cancel(dn / (e ** 2 * (1 - e))), mu)
        check(Q.degree() == 2, "paw quadratic")
        a2, a1, a0 = [sp.factor(c) for c in Q.all_coeffs()]
        res[v] = (a2, a1, a0, dd)
    # leading coefficient (1-e^2)^2 times a positive constant for every vertex
    # denominator of the derivative is c * Den^2 with c > 0: check by sign at a point
    paw_cert = {}
    # Vertex 0 (leaf): a1 > 0 and a0 > 0 on (0,1) -> derivative of fixed sign
    for v in range(4):
        a2, a1, a0, dd = res[v]
        sgn_dd = sp.sign(dd.subs({mu: 1, e: sp.Rational(1, 2)}))
        paw_cert[v] = {"a2": str(a2), "a1": str(a1), "a0": str(a0)}
        res[v] = (a2, a1, a0, sgn_dd)
    def sign_on_01(coef):
        """+1/-1 if the polynomial in epsilon has constant sign on (0,1), else 0."""
        p = upoly_e(coef)
        if roots_in_open_interval(p, 0, 1) != 0:
            return 0
        val = sp.Rational(sp.Poly(sp.together(coef), e).eval(sp.Rational(1, 2)))
        return 1 if val > 0 else -1
    for v in range(4):
        check(res[v][3] > 0, "paw derivative denominator sign")
        s2, s1 = sign_on_01(res[v][0]), sign_on_01(res[v][1])
        check(s2 != 0 and s1 == s2, "paw: a2 and a1 of one constant sign on (0,1)")
    # leaf 0: a0 has the sign of a2 on (0,1) -> no positive root -> monotone
    check(sign_on_01(res[0][2]) == sign_on_01(res[0][0]), "paw leaf monotone")
    # vertices 1,2 and hub 3: a0 has exactly one root in (0,1)
    p_a0_1 = upoly_e(res[1][2])
    check(roots_in_open_interval(p_a0_1, 0, 1) == 1, "paw v1 a0 roots")
    paw_eps1 = bracket(p_a0_1, F(783254, 10000000), F(783255, 10000000))
    # a0*a2 < 0 (one sign change, nonmonotone) exactly for eps < eps1 at vertex 1
    check(sp.Poly(res[1][2], e).eval(sp.Rational(1, 20)) * sp.Poly(res[1][0], e).eval(sp.Rational(1, 20)) < 0, "paw v1 small eps")
    p_a0_3 = upoly_e(res[3][2])
    check(roots_in_open_interval(p_a0_3, 0, 1) == 1, "paw hub a0 roots")
    paw_eps2 = bracket(p_a0_3, F(4355202, 10000000), F(4355203, 10000000))
    check(sp.Poly(res[3][2], e).eval(sp.Rational(1, 2)) * sp.Poly(res[3][0], e).eval(sp.Rational(1, 2)) < 0, "paw hub large eps")
    out["paw"] = {"derivative_quadratic_coefficients": {str(k): v for k, v in paw_cert.items()},
                  "eps1_bracket": paw_eps1, "eps2_bracket": paw_eps2,
                  "scores": [str(sp.factor(s)) for s in sp_paw]}
    # exact three-point certificates for the paw
    pts = {}
    for eps, v, mus in [(F(1, 2), 3, [F(1, 10), F(1, 2), F(1)]), (F(1, 20), 1, [F(0), F(1, 10), F(1)])]:
        vals = [B.solve_full(4, PAW_EDGES, eps, x)[v] for x in mus]
        check((vals[1] > vals[0] and vals[1] > vals[2]) or (vals[1] < vals[0] and vals[1] < vals[2]),
              "paw three-point")
        pts[str(eps)] = {"vertex": v, "mu": [str(x) for x in mus], "values": [str(x) for x in vals]}
    out["paw"]["three_point"] = pts

    # ---------------- gem
    sg = symbolic_scores(5, GEM_EDGES)
    gem = {}
    # symmetry 0<->3, 1<->2
    check(sp.cancel(sg[1] - sg[2]) == 0 and sp.cancel(sg[0] - sg[3]) == 0, "gem symmetry")
    dn = sp.expand(sp.numer(sp.together(sp.diff(sg[1], mu))))
    dd = sp.denom(sp.together(sp.diff(sg[1], mu)))
    Rg = sp.cancel(dn / (e ** 2 * (e - 1)))
    Rg = sp.Poly(sp.expand(Rg), mu)
    check(Rg.degree() == 6, "gem degree")
    lead = sp.factor(Rg.LC())
    check(sp.cancel(lead / (1 - e ** 2) ** 4) == 4, "gem leading coefficient")
    R0 = sp.expand(Rg.eval(0))
    p_R0 = upoly_e(R0)
    # R0 < 0 on (0, 1/3]: no root in (0,1/3) and value at 1/3 negative
    check(roots_in_open_interval(p_R0, 0, F(1, 3)) == 0, "gem R0 roots")
    val13 = sp.Rational(sp.Poly(R0, e).eval(sp.Rational(1, 3)))
    check(val13 < 0 and sp.Poly(R0, e).eval(0) < 0, "gem R0 sign")
    eps_g = bracket(p_R0, F(3565987, 10000000), F(3565988, 10000000))
    check(roots_in_open_interval(p_R0, 0, F(3565987, 10000000)) == 0, "gem: first root")
    # sign of d/dmu pi(1) = e^2 (e-1) R / dd ; dd > 0 check at sample, dd is a square times const
    sgn = sp.sign(dd.subs({mu: 0, e: sp.Rational(1, 5)}))
    check(sgn > 0, "gem derivative denominator sign")
    ddf = sp.factor_list(dd)
    check(all(m % 2 == 0 for (_f, m) in ddf[1]) and ddf[0] > 0, "gem derivative denominator is a positive square")
    R_eps0 = sp.factor(sp.Poly(Rg.as_expr().subs(e, 0), mu).as_expr())
    gem["R_mu_eps0"] = str(R_eps0)
    gem["R_at_mu0"] = str(sp.Poly(R0, e).as_expr())
    gem["leading_mu_coefficient"] = str(lead)
    gem["eps_g_bracket"] = eps_g
    gem["value_R0_at_one_third"] = str(val13)
    pts = {}
    for eps, mus in [(F(1, 10), [F(0), F(1, 3), F(1)]), (F(1, 4), [F(0), F(1, 10), F(1)])]:
        vals = [B.solve_full(5, GEM_EDGES, eps, x)[1] for x in mus]
        check(vals[1] > vals[0] and vals[1] > vals[2], "gem three-point")
        pts[str(eps)] = {"vertex": 1, "mu": [str(x) for x in mus], "values": [str(x) for x in vals]}
    gem["three_point"] = pts
    out["gem"] = gem

    # ---------------- diamond spread (exact): tips increase, hubs decrease
    sd = symbolic_scores(4, DIAMOND_EDGES)
    z = sd[2]
    check(sp.cancel(sd[0] - (sp.Rational(1, 2) - z)) == 0, "diamond hub")
    check(sp.cancel(sd[2] - sd[3]) == 0 and sp.cancel(sd[0] - sd[1]) == 0, "diamond symmetry")
    # hub - tip = 1/2 - 2z > 0 at a point and z' > 0 (proved in legacy identity)
    check(sp.Rational((sp.Rational(1, 2) - 2 * z).subs({mu: 1, e: sp.Rational(1, 2)})) > 0, "diamond order")
    hub_tip = sp.cancel(sp.Rational(1, 2) - 2 * z)
    nz = sp.Poly(sp.numer(sp.together(hub_tip)), mu)
    out["diamond_spread"] = str(sp.factor(hub_tip))

    # ---------------- every connected graph of order four (uniform restart)
    # star, C4, K4: mu-independent; diamond: tips strictly increasing, hubs = 1/2 - tip;
    # P4: ends strictly increasing, inner vertices strictly decreasing; paw: above.
    order4 = {}
    for name, n4, ed in [("star", 4, [(0, 1), (0, 2), (0, 3)]),
                         ("C4", 4, [(0, 1), (1, 2), (2, 3), (0, 3)]),
                         ("K4", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])]:
        sc = symbolic_scores(n4, ed)
        check(all(sp.cancel(sp.diff(x, mu)) == 0 for x in sc), name + " not invariant")
        order4[name] = "mu-independent"
    sP4 = symbolic_scores(4, [(0, 1), (1, 2), (2, 3)])
    lin = (e ** 2 - 1) * mu - e ** 2 - e - 1
    dP4 = e ** 2 * (1 - e) / (8 * lin ** 2)
    check(sp.cancel(sp.diff(sP4[0], mu) - dP4) == 0, "P4 end derivative")
    check(sp.cancel(sp.diff(sP4[3], mu) - dP4) == 0, "P4 end derivative (3)")
    check(sp.cancel(sp.diff(sP4[1], mu) + dP4) == 0 and sp.cancel(sp.diff(sP4[2], mu) + dP4) == 0,
          "P4 inner derivative")
    # lin < 0 for mu >= 0, 0 < eps < 1 (both coefficients negative), so dP4 > 0
    order4["P4"] = {"end_derivative": "eps^2(1-eps)/(8[(eps^2-1)mu-eps^2-eps-1]^2) > 0",
                    "inner_derivative": "negative of the end derivative"}
    sdm = symbolic_scores(4, DIAMOND_EDGES)
    zd = sp.cancel(sp.diff(sdm[2], mu))
    zwant = e ** 2 * (1 - e) / (12 * ((1 - e ** 2) * mu + e ** 2 + 2 * e + 2) ** 2)
    check(sp.cancel(zd - zwant) == 0, "diamond tip derivative")
    order4["diamond"] = "tip derivative eps^2(1-eps)/(12[(1-eps^2)mu+eps^2+2eps+2]^2) > 0; hubs 1/2 - tip"
    order4["paw"] = "fails iff eps in (0, eps1) or (eps2, 1)"
    out["order4_all_connected"] = order4
    return out
