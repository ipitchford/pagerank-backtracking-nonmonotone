# Assurance

This file says what has been checked and what has not. Each row is a separate kind of assurance, and no row implies another.

| Dimension | State | Evidence |
|---|---|---|
| Availability | Public GitHub release and Zenodo DOI 10.5281/zenodo.23580381 | Release assets and [SHA256SUMS](SHA256SUMS) |
| Written proofs | Complete proofs in the manuscript for every result marked P in Table 1 | [paper/paper.pdf](paper/paper.pdf). Covers, the period-class argument (Theorem 9.1) and the asymptotic statements (Theorem 11.2) are proved in prose only. |
| Exact certificates | Exact rational arithmetic throughout | [code/verify_all.py](code/verify_all.py) uses Python fractions, python-flint `fmpz_poly`, and sympy exact fields over Q(μ, ε). Every check is an explicit `RuntimeError`, so `python -O` disables nothing. |
| Certified census | Exact, for five damping values only | [results/census_full.json](results/census_full.json): all 995 connected graphs on 2–7 vertices at ε ∈ {1/10, 1/4, 1/3, 1/2, 17/20}. Interpolation with a proved degree bound (two extra check points), squarefree decomposition of W = N′Den − N Den′, and exact Sturm counts of odd-multiplicity roots on open intervals, with roots at 0 and 1 divided out first. Known-root regressions guard the counter: two synthetic responses, Atlas graph 399, and every subset of a known root set (256 cases, with multiplicities 1 and 3). Behaviour between those ε values is not certified. An earlier root-isolation routine misordered interval endpoints and misclassified 27 coordinates on (0, 1). The post-revision adversarial audit found this. The routine was replaced, and the final census agrees with that audit's independent Sturm implementation on all 33,900 coordinates. |
| Numerical results | Floating point, labelled as such | [code/prevalence.py](code/prevalence.py), [results/prevalence.json](results/prevalence.json) (Section 10, Table 4). Exact grid searches certify nothing between grid points. |
| Internal replay | Passed | After the final fix round, the producer's replays of `verify_all.py` (all certificates and the census) were byte-identical. `--jobs 1` and `--jobs 4` give identical bytes. An earlier clean replay from a separate copy also reproduced `prevalence.json`, which did not change. The frozen-archive replay is recorded in [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). |
| Semantic negative controls | Passed | [negative_controls.py](negative_controls.py) runs an unmodified baseline, then checks that each of the following is rejected: a corrupted house certificate; flipped monotonicity claims on K4, C5, the diamond, the house, the tree of Proposition 5.5 and Atlas graph 399; a wrong damping value; five plausible code bugs, two of them in the Sturm root counter; two tampered outputs; and a tampered file. |
| Implementation diversity | Partial | The first adversarial audit rebuilt the full 12-arc house system symbolically and ran independent exact attacks (covers, restarts, invariance, clique chains) without importing package code. The second audit's independent Sturm counter reused only the package's rational-response reconstruction, and checked every census coordinate against the final census. It also derived the tree and unique-extrema certificates with its own solvers. The re-review recomputed the certificates with its own exact solver. This is diversity inside one producer workflow, not independent reimplementation. |
| Internal AI-assisted review | Done | First-round multi-perspective review; AI-assisted adversarial audits by a different model (configured as OpenAI gpt-6-astra) of the first draft and of the revised version; and a verification re-review. Summary in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). This is not peer review. |
| Independent rerun or reimplementation by an unaffiliated party | Not established | |
| Specialist (human) review | Not established | |
| Editorial peer review | Not established | |
| Formal verification | Not attempted | No proof assistant was used. |
| Source fidelity | arXiv v2 only | Conjectures 3.1 and 3.2 are quoted verbatim from arXiv:2211.13353v2. The journal text (doi:10.1016/j.amc.2026.130182) was inaccessible (HTTP 403), so its wording is unverified. |
| Environment | Pinned Python packages; deterministic JSON | [ENVIRONMENT.txt](ENVIRONMENT.txt), [requirements.txt](requirements.txt) |
| Novelty | Bounded search only | [RESEARCH_GATES.json](RESEARCH_GATES.json), [SOURCES.md](SOURCES.md) |
| Citations | Metadata checked | [CITATION_AUDIT.md](CITATION_AUDIT.md) |

## Scope and limits

* The census is certified at five damping values. It says nothing certified about other values of ε, or about graphs on eight or more vertices.
* No proof assistant was used. The finite computations support the finite statements and check exact identities used in the proofs. The statements for all ε, for all orders and for asymptotic regimes rest on the written proofs.
* The conjecture wording of the journal version is unverified because of a paywall; arXiv v2 was checked verbatim.
* The configuration-model tail transfer (Section 11) is stated as an open problem with a heuristic only.
* Prevalence, excursion, shortcut-failure and ranking statistics are numerical. The only exact ranking result is the single Les Misérables certificate (Proposition 10.1).
* With leaves, results concern μ ∈ (0, ∞). At μ = 0 the source kernel is undefined at a leaf, and the package uses forced return (Section 2).

Byte-identical replay shows that the shipped outputs are exactly what the shipped programs produce. It does not by itself show that the programs implement the manuscript's definitions. That rests on the method description (Sections 7 and 12), the negative controls, and the independent exact recomputations in the internal reviews.
