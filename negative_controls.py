#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Semantic negative controls for this package.

A replay that cannot fail is not evidence.  Each control below copies code/ and results/
into a fresh temporary directory, corrupts one claim, one line of code or one output, and
requires an exact checker to REJECT it.  An unmodified baseline must pass first.  Nothing
is written beside this script or anywhere in the package directory.

The checker (run as a child process `python -B -O negative_controls.py --child ...`) uses the
package's own exact code:
  * legacy_checks.run_legacy() must reproduce results/verification.json "legacy_checks" and
    "house_certificates" exactly;
  * census.analyse_graph() on every connected graph of order <= 6 for the baseline (143 graphs,
    five damping values; order <= 5, 30 graphs, for the corrupted copies) must reproduce the corresponding records of results/census_full.json exactly;
  * the three-number house certificate (Section 1.2 of the paper) is re-derived with the
    independent full 2m-state solver bt_exact.solve_full and must show a strict dip;
  * bt_exact.known_root_regressions() must pass and reproduce the corresponding entries of
    results/verification.json "root_counter_regressions";
  * a list of exact monotonicity claims (house roof, K4, C5, diamond, paw, the five-vertex tree
    of Proposition 5.5, Atlas graph 399 on (0,1)) must agree with the exact decision procedure
    (interpolation + squarefree odd part + Sturm counts on open intervals).

Controls:
  1. baseline: the unmodified package passes the checker;
  2. certificate corrupted: house pi_0(4) at eps=1/2 claimed to be 22/125 (= pi_1(4)) -> rejected;
  3. claim flips: "K4 vertex nonmonotone", "C5 vertex nonmonotone", "diamond tip nonmonotone",
     "house roof monotone at eps=1/2", "tree leaf monotone at eps=1/2", and "Atlas 399 vertex 4
     monotone on (0,1) at eps=1/2" (the value printed before the root-counting repair) -> each rejected;
  4. wrong damping: "house roof nonmonotone at eps=1/10" (eps < eps_c ~ 0.2407) -> rejected;
  5. code mutations in bt_exact.py: sign error in W = N'Den - N Den'; keeping even- instead of
     odd-multiplicity factors; sign error in the backtracking weight (mu - 1 -> 1 - mu); Sturm
     count on (0,1) evaluated at +inf instead of 1; root at the right endpoint not divided out
     -> each rejected (mismatch with shipped results or an explicit RuntimeError);
  6. tampered outputs: the house record in results/census_full.json at eps=1/2, and one
     legacy-check count in results/verification.json -> each rejected;
  7. a tampered file is detected by SHA256SUMS (shipped sums if present, else sums of the copy).

Usage:  python -B -O negative_controls.py      Exit status 0 iff every control behaves.
Requires python-flint, sympy, networkx (see requirements.txt).  Runtime about a minute.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.abspath(__file__))

HOUSE = [(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3)]          # roof = vertex 4
K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
C5 = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]
DIAMOND = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]                 # tips 2, 3
PAW = [(0, 1), (0, 2), (1, 2), (2, 3)]
TREE = [(0, 4), (1, 3), (2, 3), (3, 4)]                            # Proposition 5.5, leaf 1
ATLAS399 = [(0, 1), (0, 3), (0, 4), (1, 3), (1, 4), (1, 6), (2, 6), (4, 5)]

# (label, n, edges, eps, vertex, claimed_nonmonotone, range) ; range "pos" = (0,inf), "01" = (0,1)
TRUE_CLAIMS = [
    ("house roof nonmonotone at eps=1/2", 5, HOUSE, "1/2", 4, True, "pos"),
    ("house roof nonmonotone at eps=1/3", 5, HOUSE, "1/3", 4, True, "pos"),
    ("house roof monotone at eps=1/10", 5, HOUSE, "1/10", 4, False, "pos"),
    ("K4 vertex monotone at eps=1/2", 4, K4, "1/2", 0, False, "pos"),
    ("C5 vertex monotone at eps=1/2", 5, C5, "1/2", 0, False, "pos"),
    ("diamond tip monotone at eps=1/2", 4, DIAMOND, "1/2", 2, False, "pos"),
    ("paw nonmonotone somewhere at eps=1/2", 4, PAW, "1/2", None, True, "pos"),
    ("tree leaf nonmonotone on (0,1) at eps=1/2", 5, TREE, "1/2", 1, True, "01"),
    ("Atlas 399 vertex 4 nonmonotone on (0,1) at eps=1/2", 7, ATLAS399, "1/2", 4, True, "01"),
]
HOUSE_CERT = {"eps": "1/2", "mu_0": "638/3615", "mu_1": "22/125", "mu_infinity": "8/45"}


# ------------------------------------------------------------------ child: the exact checker
def child(pkg, claims_json, cert_json, max_order):
    sys.path.insert(0, os.path.join(pkg, "code"))
    from fractions import Fraction as F
    import bt_exact as B
    import census
    import legacy_checks

    errs = []
    res = os.path.join(pkg, "results")
    with open(os.path.join(res, "verification.json")) as fh:
        ver = json.load(fh)
    with open(os.path.join(res, "census_full.json")) as fh:
        cen = json.load(fh)

    checks, certs, _ = legacy_checks.run_legacy()
    if checks != ver.get("legacy_checks"):
        errs.append("MISMATCH: legacy_checks")
    if certs != ver.get("house_certificates"):
        errs.append("MISMATCH: house_certificates")

    shipped = {r["atlas"]: r for r in cen["graphs"]}
    for item in census.atlas_graphs():
        if item[1] > max_order:
            continue
        rec = census.analyse_graph(item)
        if shipped.get(rec["atlas"]) != rec:
            errs.append("MISMATCH: census_full.json atlas %d" % rec["atlas"])
            break

    # three-number certificate, independent full-state solver
    cert = json.loads(cert_json)
    eps = F(cert["eps"])
    p0 = B.solve_full(5, HOUSE, eps, 0)[4]
    p1 = B.solve_full(5, HOUSE, eps, 1)[4]
    pinf = B.pi_inf(5, HOUSE, eps)[4]
    c0, c1, ci = F(cert["mu_0"]), F(cert["mu_1"]), F(cert["mu_infinity"])
    if (c0, c1, ci) != (p0, p1, pinf):
        errs.append("REJECT: house certificate values differ from exact solve")
    if not (c0 > c1 and ci > c1):
        errs.append("REJECT: claimed certificate does not show a strict dip")

    try:
        reg = B.known_root_regressions()
        shipped_reg = ver.get("root_counter_regressions") or {}
        if any(shipped_reg.get(k) != v for k, v in reg.items()):
            errs.append("MISMATCH: root_counter_regressions")
    except RuntimeError as exc:
        errs.append("REJECT: root-counter regression failed: %s" % exc)

    for label, n, edges, e, vtx, claim, rng in json.loads(claims_json):
        N, D = B.rational_response(n, [tuple(x) for x in edges], F(e))
        key = "roots_pos" if rng == "pos" else "roots_01"
        verdicts = [B.sign_change_analysis(N[i], D)[key] > 0 for i in range(n)]
        actual = any(verdicts) if vtx is None else verdicts[vtx]
        if actual != claim:
            errs.append("REJECT claim: " + label)

    for e in errs:
        print(e, file=sys.stderr)
    print("CHECKER " + ("PASS" if not errs else "FAIL"))
    return 0 if not errs else 1


# ------------------------------------------------------------------ parent
def fresh_copy():
    tmp = tempfile.TemporaryDirectory(prefix="pr-negctl-")
    dst = os.path.join(tmp.name, "pkg")
    os.makedirs(dst)
    for sub in ("code", "results"):
        shutil.copytree(os.path.join(ROOT, sub), os.path.join(dst, sub),
                        ignore=shutil.ignore_patterns("__pycache__"))
    return tmp, dst


def run_checker(pkg, claims=TRUE_CLAIMS, cert=HOUSE_CERT, max_order=5):
    cmd = [sys.executable, "-B", "-O", os.path.abspath(__file__), "--child", pkg,
           json.dumps(claims), json.dumps(cert), str(max_order)]
    return subprocess.run(cmd, cwd=os.path.dirname(pkg), capture_output=True, text=True)


def mutate(path, old, new):
    with open(path) as fh:
        text = fh.read()
    if text.count(old) != 1:
        raise SystemExit("control setup failed: pattern not unique in %s: %r" % (path, old))
    with open(path, "w") as fh:
        fh.write(text.replace(old, new))


failures = []


def expect(label, ok, detail=""):
    print("[%s] %s%s" % ("ok" if ok else "FAIL", label, (": " + detail) if detail else ""), flush=True)
    if not ok:
        failures.append(label)


def last(p):
    lines = (p.stderr or p.stdout).strip().splitlines()
    if lines and lines[0].startswith("Traceback"):
        return lines[-1]
    return lines[0] if lines else ""


def rejected(p):
    return p.returncode != 0


def main():
    # 1. baseline
    tmp, pkg = fresh_copy()
    p = run_checker(pkg, max_order=6)
    expect("baseline: unmodified package passes the exact checker (census orders <= 6)",
           p.returncode == 0 and "CHECKER PASS" in p.stdout, last(p) if p.returncode else "")
    tmp.cleanup()
    if failures:
        print("BASELINE FAILED; controls not meaningful")
        sys.exit(1)

    # 2. corrupted certificate
    tmp, pkg = fresh_copy()
    bad = dict(HOUSE_CERT, mu_0="22/125")
    p = run_checker(pkg, cert=bad)
    expect("certificate pi_0(4)=22/125 at eps=1/2 rejected", rejected(p) and "REJECT: house certificate" in p.stderr,
           last(p))
    tmp.cleanup()

    # 3-4. flipped claims and wrong damping
    flips = [
        ("claim 'K4 vertex nonmonotone' rejected", 3, None),
        ("claim 'C5 vertex nonmonotone' rejected", 4, None),
        ("claim 'diamond tip nonmonotone' rejected", 5, None),
        ("claim 'house roof monotone at eps=1/2' rejected", 0, None),
        ("wrong damping: 'house roof nonmonotone at eps=1/10' rejected", 2, None),
        ("claim 'tree leaf monotone' (Proposition 5.5) rejected", 7, None),
        ("claim 'Atlas 399 vertex 4 monotone on (0,1)' (pre-repair census value) rejected", 8, None),
    ]
    tmp, pkg = fresh_copy()
    for label, idx, _ in flips:
        claims = [list(c) for c in TRUE_CLAIMS]
        claims[idx][5] = not claims[idx][5]
        p = run_checker(pkg, claims=claims)
        expect(label, rejected(p) and ("REJECT claim: " + TRUE_CLAIMS[idx][0]) in p.stderr, last(p))
    tmp.cleanup()

    # 5. code mutations
    for label, old, new in [
        ("sign error in W = N'Den - N Den' rejected",
         "W = N.derivative() * D - N * D.derivative()", "W = N.derivative() * D + N * D.derivative()"),
        ("even- instead of odd-multiplicity factors rejected", "if e % 2 == 1:", "if e % 2 == 0:"),
        ("sign error in backtracking weight (mu - 1 -> 1 - mu) rejected",
         "return S, -eps * (mu - 1), -eps, a[j] * S", "return S, -eps * (1 - mu), -eps, a[j] * S"),
        ("Sturm count on (0,1) evaluated at +inf instead of 1 rejected",
         "vb = _variations(_signs_at(chain, None if b is None else F(b)))",
         "vb = _variations(_signs_at(chain, None))"),
        ("root at the right endpoint not divided out rejected",
         "        s, _ = strip_root(s, F(b))", "        s, _ = s, False"),
    ]:
        tmp, pkg = fresh_copy()
        mutate(os.path.join(pkg, "code", "bt_exact.py"), old, new)
        p = run_checker(pkg)
        expect(label, rejected(p), last(p))
        tmp.cleanup()

    # 6. tampered outputs
    tmp, pkg = fresh_copy()
    path = os.path.join(pkg, "results", "census_full.json")
    with open(path) as fh:
        cen = json.load(fh)
    house = [r for r in cen["graphs"] if r["atlas"] == 43]
    if len(house) != 1 or house[0]["per_eps"]["1/2"]["nonmonotone"] != [4]:
        raise SystemExit("control setup failed: house record (atlas 43) not as expected")
    house[0]["per_eps"]["1/2"]["nonmonotone"] = []
    with open(path, "w") as fh:
        json.dump(cen, fh, sort_keys=True, indent=1)
        fh.write("\n")
    p = run_checker(pkg)
    expect("tampered census_full.json (house, eps=1/2) rejected",
           rejected(p) and "MISMATCH: census_full.json atlas 43" in p.stderr, last(p))
    tmp.cleanup()

    tmp, pkg = fresh_copy()
    mutate(os.path.join(pkg, "results", "verification.json"),
           '"independent_exact_house_solves": 25', '"independent_exact_house_solves": 26')
    p = run_checker(pkg)
    expect("tampered verification.json (legacy count) rejected",
           rejected(p) and "MISMATCH: legacy_checks" in p.stderr, last(p))
    tmp.cleanup()

    # 7. SHA256SUMS
    tmp = tempfile.TemporaryDirectory(prefix="pr-negctl-sums-")
    dst = os.path.join(tmp.name, "pkg")
    shutil.copytree(ROOT, dst, ignore=shutil.ignore_patterns(".git", "__pycache__"))
    sums = os.path.join(dst, "SHA256SUMS")
    if os.path.exists(sums):
        with open(sums) as fh:
            entries = [ln.split(None, 1) for ln in fh.read().splitlines() if ln.strip()]
        entries = [(d, n.lstrip("*")) for d, n in entries]
        source = "shipped SHA256SUMS"
    else:
        entries = []
        for base, _, files in os.walk(dst):
            for f in sorted(files):
                full = os.path.join(base, f)
                with open(full, "rb") as fh:
                    entries.append((hashlib.sha256(fh.read()).hexdigest(), os.path.relpath(full, dst)))
        source = "sums of the copy (no SHA256SUMS shipped)"
    target = os.path.join("results", "verification.json")
    with open(os.path.join(dst, target), "ab") as fh:
        fh.write(b" ")
    bad = []
    for digest, name in entries:
        with open(os.path.join(dst, name), "rb") as fh:
            if hashlib.sha256(fh.read()).hexdigest() != digest:
                bad.append(name)
    expect("tampered results/verification.json detected by " + source, bad == [target], ",".join(bad))
    tmp.cleanup()

    print("ALL NEGATIVE CONTROLS BEHAVED" if not failures else "NEGATIVE CONTROLS FAILED: %s" % failures)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    if len(sys.argv) == 6 and sys.argv[1] == "--child":
        sys.exit(child(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])))
    main()
