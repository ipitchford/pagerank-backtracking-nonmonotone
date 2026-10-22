#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Replay every exact certificate and the certified census.

Usage:  python3 -B -O verify_all.py OUTDIR [--jobs J]

Writes OUTDIR/verification.json (sorted keys, no timestamps, no absolute
paths).  Never writes beside itself.  Every check is an explicit
`if not cond: raise RuntimeError`, so -O does not disable anything.
Requires: python-flint, sympy, networkx.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True

import bt_exact as B            # noqa: E402
import census                    # noqa: E402
import legacy_checks             # noqa: E402
import symbolic_checks           # noqa: E402
import criteria_checks as C      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.abspath(args.outdir)
    if outdir == here:
        raise RuntimeError("refusing to write beside the verifier")
    os.makedirs(outdir, exist_ok=True)

    res = {"status": "PASS"}
    checks, certs, cliques = legacy_checks.run_legacy()
    res["legacy_checks"] = checks
    res["house_certificates"] = certs
    res["symbolic"] = symbolic_checks.run()
    cen = census.run(jobs=args.jobs)
    # sanity: the house and the gem are the only order-5 min-degree-2 failures
    five = [r for r in cen["graphs"] if r["n"] == 5 and r["min_degree"] >= 2
            and any(r["per_eps"][e]["nonmonotone"] for e in cen["damping"])]
    B.check(len(five) == 2, "order-5 list")
    for r in cen["graphs"]:
        if r["n"] <= 3 or (r["n"] == 4 and r["min_degree"] >= 2):
            B.check(not any(r["per_eps"][e]["nonmonotone"] for e in cen["damping"]), "small graph")
    res["census"] = {"summary": cen["summary"], "damping": cen["damping"], "method": cen["method"],
                     "order5_mindeg2_nonmonotone": [{"atlas": r["atlas"], "edges": r["edges"],
                                                     "per_eps": r["per_eps"]} for r in five],
                     "order4_nonmonotone": [{"atlas": r["atlas"], "edges": r["edges"], "per_eps": r["per_eps"]}
                                            for r in cen["graphs"] if r["n"] == 4 and
                                            any(r["per_eps"][e]["nonmonotone"] for e in cen["damping"])]}
    res["endpoint_and_two_step"] = C.run_endpoint_and_twostep(cen)
    res["spread"] = C.run_spread(cen)
    res["rank_reversal"] = C.run_rank_reversal()
    res["rank_search_order7"] = C.run_rank_search(cen)
    res["root_counter_regressions"] = C.run_root_regressions(cen)
    res["tree_counterexample"] = C.run_tree()
    res["unique_extrema_spread_counterexample"] = C.run_unique_spread()
    res["tails"] = C.run_tails()
    res["tails"]["degree_bound_vertices_checked"] = C.run_degree_bound(cen)
    res["arithmetic"] = "exact rational (fractions, python-flint fmpz_poly, sympy QQ(mu,eps))"
    path = os.path.join(outdir, "verification.json")
    with open(path, "w") as fh:
        json.dump(res, fh, sort_keys=True, indent=1)
        fh.write("\n")
    with open(os.path.join(outdir, "census_full.json"), "w") as fh:
        json.dump(cen, fh, sort_keys=True, indent=1)
        fh.write("\n")
    print(json.dumps({"status": res["status"], "legacy": checks,
                      "census": {k: v["graphs"] for k, v in cen["summary"].items()}}, sort_keys=True))


if __name__ == "__main__":
    main()
