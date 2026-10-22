# Agent navigation

## Identity

- Release slug: `pagerank-backtracking-nonmonotone`. Version `1.0.0-candidate`, tag `v1.0.0-candidate`. DOI 10.5281/zenodo.23580381; concept DOI 10.5281/zenodo.23580380.
- Repository: https://github.com/ipitchford/pagerank-backtracking-nonmonotone. Release page: https://evidencepress.org/releases/pagerank-backtracking-nonmonotone/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](paper/paper.pdf) (17 pages) and [LaTeX source](paper/paper.tex). The figure is [paper/house_response.pdf](paper/house_response.pdf) (illustration only).
- Machine-readable claims: [CLAIMS.json](CLAIMS.json). Overview for humans: [README.md](README.md).
- Assurance limits: [ASSURANCE.md](ASSURANCE.md). Internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
- Provenance: [PROVENANCE.md](PROVENANCE.md). Sources: [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [RESEARCH_GATES.json](RESEARCH_GATES.json).
- Citation: [CITATION.cff](CITATION.cff) and [Zenodo metadata](.zenodo.json). Licences: [LICENSES.md](LICENSES.md), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

Object: on a finite connected simple graph, a random walk on directed edges. The immediate return has weight μ ≥ 0 and the other continuations weight 1. Restart has probability 1 − ε and is spread as u_ij = v_i/d_i. Vertex score π_μ(i) = Σ_j q_μ(i, j) (tail projection). μ = 0 is non-backtracking, μ = 1 is ordinary PageRank, and μ → ∞ gives π_∞ = (v + εPv)/(1 + ε). Source: Glover, Jones, Kempton and Oveson, Appl. Math. Comput. 531 (2026) 130182; arXiv:2211.13353v2, Conjectures 3.1 and 3.2.

- Conjecture 3.1 (each coordinate monotone in μ, and the max − min spread nondecreasing) is **false**. Each claim below gives its status and its location in [paper/paper.pdf](paper/paper.pdf).
  - House graph, ε = 1/2: π_0(4) = 638/3615, π_1(4) = 22/125, π_∞(4) = 8/45, so both endpoints lie above the middle. Section 1.2, exact.
  - Roof classification: monotone iff ε ≤ ε_c ≈ 0.2406813. Otherwise one strict minimum, at μ = 1 iff ε = 4/9. Theorems 3.1 and 4.1, **proved**.
  - The endpoint shortcut is wrong iff ε > ε_0 ≈ 0.3291038. Proposition 4.3, **proved**.
- Minimality.
  - Minimum degree ≥ 2: order 5 is smallest for every ε (house ∪ gem). Theorem 5.2, **proved**.
  - Leaves allowed: the paw (order 4) fails iff ε ∈ (0, ε_1) ∪ (ε_2, 1), with ε_1 ≈ 0.0783 and ε_2 ≈ 0.4355. The smallest order is 4 there and 5 on [ε_1, ε_2]. Proposition 5.3 and Corollary 5.4, **proved**.
  - A five-vertex tree has a nonmonotone leaf at ε = 1/2. Proposition 5.5, exact certificate.
  - Covers give failures at every order 5k. Lemma 5.6 and Corollary 5.7, **proved**.
- Why vertices dip: a small-damping two-step criterion (Theorem 6.1), the exact slope dπ/dμ at μ = 1 (Proposition 6.2), and three one-sided endpoint tests (Corollary 6.3). All **proved**.
- Census: all 995 connected graphs on 2–7 vertices, ε ∈ {1/10, 1/4, 1/3, 1/2, 17/20}. Section 7, Table 3, **certified exact** at those five ε only. Data in [results/census_full.json](results/census_full.json); summary in [results/verification.json](results/verification.json).
- Spread clause: fails both as printed and when reversed, even with unique extreme vertices. Propositions 8.1 and 8.2, exact certificates.
- Invariance: the full vector is μ-invariant iff P²v = v. Theorem 9.1, **proved**.
- Prevalence: 10–41% of vertices on standard networks. Section 10, Table 4, **numerical**; [results/prevalence.json](results/prevalence.json). One exact transient rank reversal is given in Proposition 10.1.
- Conjecture 3.2 (degree tail implies a scale-free π_∞) fails for every tail index α > 1. Theorem 11.2, **proved**; the mechanism is adapted from Henning, van der Hofstad and Litvak. Proposition 11.1 gives a one-sided bound.

Do not read:
- the census as certified for damping values other than the five listed;
- grid searches as certificates of absence;
- the prevalence numbers as exact;
- "four vertices with leaves" as holding for every ε;
- any result as covering directed or weighted graphs;
- the configuration-model tail transfer as proved;
- internal AI-assisted review as peer review.

## Code

- [code/bt_exact.py](code/bt_exact.py): exact solver (block elimination, plus an independent full 2m-state solver), interpolation, squarefree odd part, exact Sturm root counts on open intervals, and known-root regressions.
- [code/census.py](code/census.py): the certified census. [code/symbolic_checks.py](code/symbolic_checks.py): P3, paw, gem and diamond over Q(μ, ε).
- [code/criteria_checks.py](code/criteria_checks.py): slope formula, endpoint and two-step tests, spread (including the unique-extrema witness), the tree counterexample, the rank reversal and the finite exact rank search, degree tails. [code/legacy_checks.py](code/legacy_checks.py): house, diamond, cover and clique-chain identities.
- [code/verify_all.py](code/verify_all.py): runs all of the above. [code/prevalence.py](code/prevalence.py): the numerical study. [code/make_figure.py](code/make_figure.py): the figure.

## Replay

Python 3.13 with the pins in [requirements.txt](requirements.txt):

```sh
sha256sum -c SHA256SUMS                         # or: shasum -a 256 -c SHA256SUMS
python -B -O code/verify_all.py OUT --jobs 4    # all exact certificates + census, a few minutes
cmp OUT/verification.json results/verification.json && cmp OUT/census_full.json results/census_full.json
python -B -O code/prevalence.py OUT             # numerical study, about 1 minute
python -B -O negative_controls.py               # semantic negative controls, under a minute
```

Results of the release replay: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). Environment: [ENVIRONMENT.txt](ENVIRONMENT.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml).

## Useful next research question

Which graphs have all coordinates monotone for every ε? Trees are not such a class (Proposition 5.5). A certified extension of the census to order 8, or to the whole damping interval, would be the cheapest way to sharpen the picture (Section 13, open problems 1 and 4).
