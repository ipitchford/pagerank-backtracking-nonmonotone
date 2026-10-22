# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Exact arithmetic for backtracking-weighted PageRank.

Model: finite connected simple graph; directed-edge states (i, j);
transition (i, j) -> (j, k) with weight (1 + (mu - 1)[k == i]) / (d_j - 1 + mu);
a degree-one vertex always returns (the continuous extension to mu = 0 of the
source kernel, which equals the source kernel for every mu > 0);
restart u(i, j) = v_i / d_i; q = (1 - eps) u + eps Q_mu q; pi(i) = sum_j q(i, j).

All arithmetic is exact (fractions.Fraction, python-flint fmpz_poly).
No assert statements are used, so the module behaves identically under -O.
"""
from fractions import Fraction as F
from math import lcm
import flint


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


# ---------------------------------------------------------------- linear algebra
def solve_det(A, b):
    """Gaussian elimination over Q. Returns (x, det A)."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    det = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        check(p is not None, "singular system")
        if p != c:
            M[c], M[p] = M[p], M[c]
            det = -det
        piv = M[c][c]
        det *= piv
        inv = 1 / piv
        for r in range(c + 1, n):
            f = M[r][c]
            if f:
                f *= inv
                Mr, Mc = M[r], M[c]
                for k in range(c, n + 1):
                    Mr[k] -= f * Mc[k]
    x = [F(0)] * n
    for r in range(n - 1, -1, -1):
        s = M[r][n] - sum(M[r][k] * x[k] for k in range(r + 1, n))
        x[r] = s / M[r][r]
    return x, det


# ---------------------------------------------------------------- graph helpers
def adjacency(n, edges):
    nb = [[] for _ in range(n)]
    for a, b in edges:
        nb[a].append(b)
        nb[b].append(a)
    return [sorted(x) for x in nb]


def solve_point(n, edges, eps, mu, v=None):
    """Exact vertex scores and the polynomial-denominator value Den(mu).

    Den(mu) = det M_full(mu), where M_full is the 2m x 2m edge system whose
    row (j,k) is multiplied by d_j - 1 + mu (non-leaf tail) so that every
    entry has degree <= 1 in mu.  Computed as det(X) det(I + Z X^{-1} Y)
    (block elimination of the head masses H_j = sum_i q(i, j)).
    """
    eps, mu = F(eps), F(mu)
    nb = adjacency(n, edges)
    d = [len(x) for x in nb]
    if v is None:
        v = [F(1, n)] * n
    a = [(1 - eps) * v[j] / d[j] for j in range(n)]

    def rowdata(j):
        # (S_j, off_j, y_j, r_j): row (j,k): S q(j,k) + off q(k,j) + y H_j = r
        if d[j] == 1:
            return F(1), -eps, F(0), a[j]
        S = d[j] - 1 + mu
        return S, -eps * (mu - 1), -eps, a[j] * S

    rd = [rowdata(j) for j in range(n)]
    A = [[F(int(i == k)) for k in range(n)] for i in range(n)]
    c = [F(0)] * n
    detX = F(1)
    blocks = []
    for (j, k) in edges:
        Sj, oj, yj, rj = rd[j]
        Sk, ok, yk, rk = rd[k]
        det = Sj * Sk - oj * ok
        check(det != 0, "zero block")
        detX *= det
        # q(j,k) = [Sk (rj - yj Hj) - oj (rk - yk Hk)] / det
        # q(k,j) = [Sj (rk - yk Hk) - ok (rj - yj Hj)] / det
        blocks.append((j, k, Sj, oj, yj, rj, Sk, ok, yk, rk, det))
        # H_k gets q(j,k); H_j gets q(k,j).  Equation: H - sum(q) = 0
        c[k] += (Sk * rj - oj * rk) / det
        A[k][j] += Sk * yj / det
        A[k][k] -= oj * yk / det
        c[j] += (Sj * rk - ok * rj) / det
        A[j][k] += Sj * yk / det
        A[j][j] -= ok * yj / det
    H, detA = solve_det(A, c)
    pi = [F(0)] * n
    for (j, k, Sj, oj, yj, rj, Sk, ok, yk, rk, det) in blocks:
        qjk = (Sk * (rj - yj * H[j]) - oj * (rk - yk * H[k])) / det
        qkj = (Sj * (rk - yk * H[k]) - ok * (rj - yj * H[j])) / det
        pi[j] += qjk
        pi[k] += qkj
    check(sum(pi) == 1, "mass not one")
    return pi, detX * detA


def solve_full(n, edges, eps, mu, v=None):
    """Independent check: full 2m-state edge system, transition from (eq. kernel)."""
    eps, mu = F(eps), F(mu)
    nb = adjacency(n, edges)
    d = [len(x) for x in nb]
    if v is None:
        v = [F(1, n)] * n
    arcs = [(i, j) for i in range(n) for j in nb[i]]
    idx = {e: t for t, e in enumerate(arcs)}
    N = len(arcs)
    M = [[F(int(r == s)) for s in range(N)] for r in range(N)]
    b = [F(0)] * N
    for (i, j) in arcs:
        b[idx[(i, j)]] = (1 - eps) * v[i] / d[i]
        for k in nb[j]:
            if d[j] == 1:
                w = F(1)
            else:
                w = (1 + (mu - 1) * (k == i)) / (d[j] - 1 + mu)
            M[idx[(j, k)]][idx[(i, j)]] -= eps * w
    q, _ = solve_det(M, b)
    pi = [F(0)] * n
    for (i, j) in arcs:
        pi[i] += q[idx[(i, j)]]
    return pi


def pi_one(n, edges, eps, v=None):
    return solve_point(n, edges, eps, 1, v)[0]


def pi_inf(n, edges, eps, v=None):
    eps = F(eps)
    nb = adjacency(n, edges)
    d = [len(x) for x in nb]
    if v is None:
        v = [F(1, n)] * n
    Pv = [sum(v[j] / d[j] for j in nb[i]) for i in range(n)]
    return [(v[i] + eps * Pv[i]) / (1 + eps) for i in range(n)]


# ---------------------------------------------------------------- polynomials
def interpolate(xs, ys):
    """Newton interpolation over Q; returns low->high coefficient list."""
    m = len(xs)
    coef = list(ys)
    for j in range(1, m):
        for i in range(m - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    poly = [F(0)] * m
    # Horner on Newton form
    out = [coef[m - 1]]
    for i in range(m - 2, -1, -1):
        # out = out*(x - xs[i]) + coef[i]
        new = [F(0)] * (len(out) + 1)
        for t, cc in enumerate(out):
            new[t + 1] += cc
            new[t] -= cc * xs[i]
        new[0] += coef[i]
        out = new
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    poly[:len(out)] = out
    return out


def to_fmpz(coeffs):
    den = 1
    for c in coeffs:
        den = lcm(den, F(c).denominator)
    return flint.fmpz_poly([int(F(c) * den) for c in coeffs])


def rational_response(n, edges, eps, v=None):
    """Exact (N_i, Den) as fmpz_poly with pi_mu(i) = N_i(mu)/Den(mu), mu >= 0.

    Degree bound 2m is proved (Cramer's rule on the degree-one-entry system);
    an extra exact evaluation at two further points is checked.
    """
    m = len(edges)
    xs = [F(t) for t in range(2 * m + 1)]
    vals = [solve_point(n, edges, eps, x, v) for x in xs]
    Den = interpolate(xs, [dv for (_, dv) in vals])
    Ns = [interpolate(xs, [p[i] * dv for (p, dv) in vals]) for i in range(n)]
    # independent consistency checks at non-node points
    for x in (F(2 * m + 3), F(1, 3)):
        p, dv = solve_point(n, edges, eps, x, v)
        dval = sum(c * x ** t for t, c in enumerate(Den))
        check(dval == dv, "Den interpolation mismatch")
        for i in range(n):
            nval = sum(c * x ** t for t, c in enumerate(Ns[i]))
            check(nval == p[i] * dv, "N interpolation mismatch")
    # put on a common integer scale
    allc = Den + [c for N in Ns for c in N]
    L = 1
    for c in allc:
        L = lcm(L, c.denominator)
    D = flint.fmpz_poly([int(c * L) for c in Den])
    Nz = [flint.fmpz_poly([int(c * L) for c in N]) for N in Ns]
    return Nz, D


def peval(p, x):
    x = F(x)
    s = F(0)
    for c in reversed([int(c) for c in p.coeffs()]):
        s = s * x + c
    return s


# ---------------------------------------------------------------- root counting
# Exact Sturm counting on a requested open interval.  Roots at the finite
# endpoints are first removed algebraically (exact division), so the Sturm
# count V(a) - V(b) counts exactly the roots in the open interval (a, b).
def _q(x):
    x = F(x)
    return flint.fmpq(x.numerator, x.denominator)


def _to_fmpq_poly(p):
    if isinstance(p, flint.fmpq_poly):
        return p
    return flint.fmpq_poly([int(c) for c in p.coeffs()])


def squarefree_part(p):
    """Product of the distinct irreducible-up-to-multiplicity factors of p (fmpq_poly)."""
    p = _to_fmpq_poly(p)
    check(not p.is_zero(), "zero polynomial")
    if p.degree() <= 0:
        return flint.fmpq_poly([1])
    g = p.gcd(p.derivative())
    s, r = divmod(p, g)
    check(r.is_zero(), "squarefree division")
    return s


def _sign(v):
    return (v > 0) - (v < 0)


def _sturm_chain(p):
    chain = [p, p.derivative()]
    while not chain[-1].is_zero() and chain[-1].degree() > 0:
        _, r = divmod(chain[-2], chain[-1])
        if r.is_zero():
            break
        chain.append(-r)
    return [c for c in chain if not c.is_zero()]


def _variations(signs):
    s = [x for x in signs if x != 0]
    return sum(1 for a, b in zip(s, s[1:]) if a != b)


def _signs_at(chain, x):
    if x is None:  # +infinity: sign of leading coefficient
        return [_sign(F(int(c.leading_coefficient().p), int(c.leading_coefficient().q))) for c in chain]
    xq = _q(x)
    return [_sign(F(int(c(xq).p), int(c(xq).q))) for c in chain]


def strip_root(p, a):
    """Divide out (x - a) from a squarefree p if p(a) == 0. Returns (p', removed)."""
    a = _q(a)
    if p(a) != 0:
        return p, False
    lin = flint.fmpq_poly([-a, 1])
    s, r = divmod(p, lin)
    check(r.is_zero(), "root division")
    check(s(a) != 0, "root of multiplicity > 1 in squarefree polynomial")
    return s, True


def roots_open(p, a, b=None):
    """Number of distinct real roots of p in the open interval (a, b); b=None is +inf."""
    a = F(a)
    check(b is None or F(b) > a, "empty interval")
    s = squarefree_part(p)
    s, _ = strip_root(s, a)
    if b is not None:
        s, _ = strip_root(s, F(b))
    if s.degree() <= 0:
        return 0
    chain = _sturm_chain(s)
    va = _variations(_signs_at(chain, a))
    vb = _variations(_signs_at(chain, None if b is None else F(b)))
    k = va - vb
    check(k >= 0, "negative Sturm count")
    return k


def odd_part(W):
    """Squarefree product of the factors of W with odd multiplicity."""
    _, facs = W.factor_squarefree()
    O = flint.fmpz_poly([1])
    for g, e in facs:
        if e % 2 == 1:
            O = O * g
    return O


def sign_change_analysis(N, D):
    """Decide exactly whether N/D (D > 0 on [0, inf)) is monotone on (0, inf), (0, 1).

    W = N' D - N D'.  pi' has the sign of W.  pi is non-monotone on an open
    interval iff W has a root of odd multiplicity in it.  Odd-multiplicity
    roots are the roots of the squarefree odd part O of W; they are counted
    by an exact Sturm sequence on each open interval after any root at 0 and
    at 1 has been divided out exactly.  A root at 1 counts for (0, inf) and
    not for (0, 1); a root at 0 counts for neither.
    """
    W = N.derivative() * D - N * D.derivative()
    if W.is_zero():
        return {"constant": True, "roots_pos": 0, "roots_01": 0, "root_at_1": False}
    O = _to_fmpq_poly(odd_part(W))
    pos = roots_open(O, 0, None)
    r01 = roots_open(O, 0, 1)
    at1 = O.degree() > 0 and O(_q(1)) == 0
    return {"constant": False, "roots_pos": pos, "roots_01": r01, "root_at_1": bool(at1),
            "W_degree": W.degree()}


def known_root_regressions():
    """Exact regressions for the root counter (returned for the verification record)."""
    from itertools import combinations
    out = {}
    # N = 4 mu^3 - 15 mu^2 + 12 mu, D = 1: N' = 6 (2mu - 1)(mu - 2)
    r = sign_change_analysis(flint.fmpz_poly([0, 12, -15, 4]), flint.fmpz_poly([1]))
    check(r["roots_pos"] == 2 and r["roots_01"] == 1, "regression 4mu^3-15mu^2+12mu")
    out["4mu^3-15mu^2+12mu"] = {"roots_pos": r["roots_pos"], "roots_01": r["roots_01"]}
    # N = 4 mu^3 - 3 mu^2, D = 1: N' = 6 mu (2 mu - 1), root at 0 excluded
    r = sign_change_analysis(flint.fmpz_poly([0, 0, -3, 4]), flint.fmpz_poly([1]))
    check(r["roots_pos"] == 1 and r["roots_01"] == 1, "regression 4mu^3-3mu^2")
    out["4mu^3-3mu^2"] = {"roots_pos": r["roots_pos"], "roots_01": r["roots_01"]}
    # every subset of known simple roots (including 0 and 1 and negatives)
    pool = [F(-2), F(-1, 3), F(0), F(1, 2), F(1), F(2), F(5, 2)]
    n = 0
    for k in range(len(pool) + 1):
        for S in combinations(pool, k):
            p = flint.fmpq_poly([1])
            for x in S:
                p = p * flint.fmpq_poly([-_q(x), 1])
            for mult in (1, 3):
                pp = p ** mult
                e_pos = sum(1 for x in S if x > 0)
                e_01 = sum(1 for x in S if 0 < x < 1)
                check(roots_open(pp, 0, None) == e_pos, "known-root subset (0,inf)")
                check(roots_open(pp, 0, 1) == e_01, "known-root subset (0,1)")
                check(roots_open(pp, F(-1), F(1, 2)) == sum(1 for x in S if -1 < x < F(1, 2)),
                      "known-root subset (-1,1/2)")
                n += 1
    out["known_root_subsets_checked"] = n
    return out
