# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Exact certificates for the house graph, the diamond, covers, an invariance
sanity check and clique chains.  Assertions were
replaced by explicit checks so that the module is valid under python -O."""
from fractions import Fraction as F
from itertools import permutations


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


class Poly(dict):
    """Sparse polynomial over Q in (mu, epsilon)."""
    def __init__(self, value=0):
        if isinstance(value, dict):
            super().__init__({k: F(v) for k, v in value.items() if v})
        else:
            super().__init__({(0, 0): F(value)} if value else {})
    def __add__(self, other):
        result = dict(self)
        for key, value in Poly(other).items():
            result[key] = result.get(key, F(0)) + value
        return Poly(result)
    __radd__ = __add__
    def __neg__(self): return Poly({key: -v for key, v in self.items()})
    def __sub__(self, other): return self + -Poly(other)
    def __rsub__(self, other): return Poly(other) + -self
    def __mul__(self, other):
        result = {}
        for (i, j), a in self.items():
            for (k, l), b in Poly(other).items():
                key = (i+k, j+l)
                result[key] = result.get(key, F(0)) + a*b
        return Poly(result)
    __rmul__ = __mul__
    def __pow__(self, n):
        result = Poly(1)
        for _ in range(n): result = result*self
        return result
    def dx(self): return Poly({(i-1, j): i*v for (i, j), v in self.items() if i})
    def eval(self, mu, eps):
        return sum((v*mu**i*eps**j for (i, j), v in self.items()), F(0))
    def at_mu(self, mu):
        result = Poly(0)
        for (i, j), v in self.items(): result += Poly({(0, j): v*mu**i})
        return result

def det(matrix):
    """Leibniz determinant, used only on 3x3 and 6x6 sparse matrices."""
    n = len(matrix)
    result = Poly(0)
    for perm in permutations(range(n)):
        if any(not matrix[i][perm[i]] for i in range(n)): continue
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n))
        term = Poly(sign)
        for i in range(n): term *= matrix[i][perm[i]]
        result += term
    return result

def solve(matrix, rhs):
    n = len(rhs)
    aug = [[F(x) for x in row]+[F(b)] for row, b in zip(matrix, rhs)]
    for i in range(n):
        pivot = next(k for k in range(i, n) if aug[k][i])
        aug[i], aug[pivot] = aug[pivot], aug[i]
        a = aug[i][i]
        aug[i] = [x/a for x in aug[i]]
        for k in range(n):
            if k != i:
                a = aug[k][i]
                aug[k] = [x-a*y for x, y in zip(aug[k], aug[i])]
    return [row[-1] for row in aug]

def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for i, j in edges:
        check(i != j and j not in adj[i], "legacy check failed")
        adj[i].append(j); adj[j].append(i)
    return [sorted(row) for row in adj]

HOUSE_EDGES = [(0,1), (0,3), (0,4), (1,2), (1,4), (2,3)]
HOUSE = adjacency(5, HOUSE_EDGES)

def edge_pagerank(adj, mu, eps, restart=None):
    """Definition-level construction, independent of all quotient formulas."""
    mu, eps = F(mu), F(eps)
    n = len(adj)
    if restart is None: restart = [F(1, n)]*n
    check(sum(restart) == 1, "legacy check failed")
    arcs = [(i,j) for i in range(n) for j in adj[i]]
    index = {arc: e for e, arc in enumerate(arcs)}
    m = len(arcs)
    transition = [[F(0)]*m for _ in range(m)]
    u = [restart[i]/len(adj[i]) for i,j in arcs]
    for e,(i,j) in enumerate(arcs):
        for k in adj[j]:
            transition[index[j,k]][e] = (mu if i == k else F(1))/(len(adj[j])-1+mu)
    check(all(sum(transition[i][j] for i in range(m)) == 1 for j in range(m)), "legacy check failed")
    matrix = [[F(i == j)-eps*transition[i][j] for j in range(m)] for i in range(m)]
    rhs = [(1-eps)*x for x in u]
    q = solve(matrix, rhs)
    check(all(sum(matrix[i][j]*q[j] for j in range(m)) == rhs[i] for i in range(m)), "legacy check failed")
    check(sum(q) == 1 and min(q) >= 0, "legacy check failed")
    p = [F(0)]*n
    for (i,j), val in zip(arcs, q): p[i] += val
    check(sum(p) == 1, "legacy check failed")
    return p, q

def coefficients(e):
    n = [2*e**5+4*e**4+9*e**3+10*e**2+8*e+12,
         6*e**4+13*e**3+25*e**2+22*e+24,
         6*e**3+14*e**2+19*e+15,
         (1+e)*(3+2*e)]
    d = [e**5+2*e**4+3*e**3+4*e**2+4*e+4,
         3*e**4+6*e**3+10*e**2+10*e+8,
         3*e**3+6*e**2+8*e+5,
         (1+e)**2]
    return n,d

def house_formula(mu, eps):
    n,d = coefficients(eps)
    s = (1-eps)*mu
    return sum(n[i]*s**i for i in range(4))/(15*sum(d[i]*s**i for i in range(4)))

def cover(k):
    edges = []
    for i,j in HOUSE_EDGES:
        for t in range(k): edges.append((k*i+t, k*j+(t+int((i,j)==(0,1)))%k))
    return adjacency(5*k, edges)

def clique_chain(L):
    """Deterministic assortative construction in the degree-tail theorem."""
    edges, portals, sizes = [], [], []
    n = 0
    for m in range(3, L+1):
        for _ in range(L**4//m**4):
            portals.append(n); sizes.append(m)
            edges.extend((i,j) for i in range(n,n+m) for j in range(i+1,n+m))
            n += m
    edges.extend(zip(portals,portals[1:]))
    return adjacency(n,edges), portals, sizes

def run_legacy():
    checks = {}
    x = Poly({(1,0):1}); e = Poly({(0,1):1}); t = 1-e
    ra = x+2; rb = x+1
    # Arc orbits a,b,c,d,h,f; see manuscript equation (quotient).
    M = [[ra-e*x,0,-e,0,0,-e],[-e,ra,-e*x,0,0,-e],
         [0,-e*x,rb,-e,0,0],[0,-e,0,rb-e*x,0,0],
         [-e,0,-e,0,ra,-e*x],[0,0,0,0,-e,1]]
    M = [[Poly(v) for v in row] for row in M]
    rhs30 = [2*t*ra,2*t*ra,3*t*rb,3*t*rb,2*t*ra,3*t]
    determinant = det(M)
    cramer_numerator = det([row[:5]+[b] for row,b in zip(M,rhs30)])
    nc,dc = coefficients(e)
    N = sum((nc[i]*(t*x)**i for i in range(4)), Poly(0))
    D = sum((dc[i]*(t*x)**i for i in range(4)), Poly(0))
    factor = (x+1)*(x+2)*t
    check(determinant == factor*D, "legacy check failed")
    check(cramer_numerator == factor*N, "legacy check failed")
    checks['bivariate_cramer_identity'] = True
    A = e**6-2*e**5-9*e**4-11*e**3-16*e**2-12*e+4
    B = 4*e**5-4*e**4-14*e**3-18*e**2-16*e+12
    C = 6*e**4-e+13
    s = t*x
    R = A+B*s+C*s**2+2*(1+e)*(2*e**2+3)*s**3+(1+e)**2*s**4
    check(N.dx()*D-N*D.dx() == e**2*t*R, "legacy check failed")
    checks['bivariate_derivative_identity'] = True
    check((N.dx()*D-N*D.dx()).at_mu(1) == 9*e**2*t*(4-9*e), "legacy check failed")
    checks['turning_location_at_mu_one'] = True
    Fpoly = 2*e**4+e**3+6*e**2+10*e-4
    check(N.at_mu(0)*(6+e-e**2)-3*D.at_mu(0)*(6-e) == e**2*t*Fpoly, "legacy check failed")
    check((3+2*e)*(6+e-e**2)-3*(1+e)*(6-e) == 2*e**2*t, "legacy check failed")
    checks['endpoint_sign_identities'] = True
    # Diamond formula and derivative certificate.
    Md = [[ra-e*x,0,-2*e],[-e,ra,-e*(x+1)],[0,-e,1]]
    Md = [[Poly(v) for v in row] for row in Md]
    rd = [2*t*ra,2*t*ra,3*t]
    Nd = (3-e-2*e**2)*x+2*e**2+4*e+6
    Dd = (1-e**2)*x+e**2+2*e+2
    check(det(Md) == ra*t*Dd, "legacy check failed")
    check(det([row[:2]+[b] for row,b in zip(Md,rd)]) == ra*t*Nd, "legacy check failed")
    check(Nd.dx()*Dd-Nd*Dd.dx() == e**2*t, "legacy check failed")
    checks['diamond_monotonicity_identity'] = True
    sample_count = 0
    for eps in [F(1,10),F(1,3),F(1,2),F(17,20),F(9,10)]:
        for mu in [F(0),F(1,4),F(1,2),F(1),F(10)]:
            p,q = edge_pagerank(HOUSE,mu,eps)
            check(p[4] == house_formula(mu,eps), "legacy check failed")
            sample_count += 1
    checks['independent_exact_house_solves'] = sample_count
    certificate = {}
    for eps in [F(1,3),F(1,2),F(17,20)]:
        vals = [house_formula(mu,eps) for mu in [F(0),F(1,2),F(1)]]
        inf = (3+2*eps)/(15*(1+eps))
        certificate[str(eps)] = {'mu_0':str(vals[0]),'mu_half':str(vals[1]),
                                'mu_1':str(vals[2]),'mu_infinity':str(inf),
                                'mu_0_minus_mu_1':str(vals[0]-vals[2]),
                                'infinity_minus_mu_1':str(inf-vals[2])}
    vals = [house_formula(mu,F(1,3)) for mu in [F(0),F(1,2),F(1)]]
    check(vals[0] > vals[1] < vals[2], "legacy check failed")
    certificate['unit_interval_gaps'] = [str(vals[0]-vals[1]), str(vals[2]-vals[1])]
    checks['strict_unit_interval_counterexample'] = True
    for k in [2,3]:
        adj = cover(k)
        reached={0}; stack=[0]
        while stack:
            for j in adj[stack.pop()]:
                if j not in reached: reached.add(j);stack.append(j)
        check(len(reached)==5*k and min(map(len,adj))==2 and max(map(len,adj))==3, "legacy check failed")
        for mu in [F(0),F(1,2),F(1)]:
            p,_ = edge_pagerank(adj,mu,F(1,3))
            base,_ = edge_pagerank(HOUSE,mu,F(1,3))
            check(all(p[k*i+j] == base[i]/k for i in range(5) for j in range(k)), "legacy check failed")
    checks['independent_exact_cover_solves'] = 6
    deg = [F(len(row),12) for row in HOUSE]
    for mu in [F(0),F(1,3),F(1),F(10)]:
        p,_=edge_pagerank(HOUSE,mu,F(2,3),deg)
        check(p == deg, "legacy check failed")
    checks['degree_restart_invariance_sanity_check'] = True
    # Exact rational brackets; uniqueness follows the proof in the manuscript.
    check(A.eval(0,F(2406813,10000000)) > 0 > A.eval(0,F(2406814,10000000)), "legacy check failed")
    check(Fpoly.eval(0,F(3291037,10000000)) < 0 < Fpoly.eval(0,F(3291038,10000000)), "legacy check failed")
    checks['exact_threshold_brackets'] = True
    clique_certificates = []
    for L in [4,6,10,14]:
        adj,portals,sizes = clique_chain(L)
        n = len(adj); degrees = list(map(len,adj))
        check(n == sum(m*(L**4//m**4) for m in range(3,L+1)), "legacy check failed")
        check(sum(degrees)//2 == sum(m*(m-1)//2*(L**4//m**4) for m in range(3,L+1))+len(sizes)-1, "legacy check failed")
        reached={0}; stack=[0]
        while stack:
            for j in adj[stack.pop()]:
                if j not in reached: reached.add(j); stack.append(j)
        check(len(reached) == n and min(degrees) >= 2, "legacy check failed")
        expected = {}
        for t,m in enumerate(sizes):
            b = int(t>0)+int(t+1<len(sizes))
            expected[m-1] = expected.get(m-1,0)+m-1
            expected[m-1+b] = expected.get(m-1+b,0)+1
        observed = {d:degrees.count(d) for d in set(degrees)}
        check(observed == expected, "legacy check failed")
        h = [sum((F(1,degrees[j]) for j in row),F(0)) for row in adj]
        check(max(h) <= F(5,3), "legacy check failed")
        check(all(h[i] <= 1 for i in range(n) if i not in set(portals)), "legacy check failed")
        for eps in [F(1,3),F(1,2),F(17,20)]:
            scaled = [(1+eps*x)/(1+eps) for x in h]
            check(sum(scaled) == n, "legacy check failed")
            check(max(scaled) <= (1+F(5,3)*eps)/(1+eps) < F(4,3), "legacy check failed")
        clique_certificates.append({'L':L,'vertices':n,'edges':sum(degrees)//2,
                                   'cliques':len(sizes),'max_reciprocal_degree_sum':str(max(h))})
    checks['exact_clique_chain_graphs'] = len(clique_certificates)
    return checks, certificate, clique_certificates
