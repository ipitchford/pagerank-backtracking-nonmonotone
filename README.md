# Backtracking-weighted PageRank is not coordinatewise monotone

**Release 1.0.0-candidate (23 October 2026). Unrefereed research candidate.**
DOI: [10.5281/zenodo.23580381](https://doi.org/10.5281/zenodo.23580381) (this version) ·
concept DOI [10.5281/zenodo.23580380](https://doi.org/10.5281/zenodo.23580380) (all versions) ·
repository <https://github.com/ipitchford/pagerank-backtracking-nonmonotone> (tag `v1.0.0-candidate`) ·
Evidence Press page <https://evidencepress.org/releases/pagerank-backtracking-nonmonotone/>.

Checked by internal AI-assisted review only (multi-perspective review, adversarial audits by a different model
before and after revision, and a verification re-review): not peer review, not independent reproduction, not formal verification.
See [STATUS.md](STATUS.md) and [ASSURANCE.md](ASSURANCE.md).

Author: Anonymous. Text CC0-1.0; code MIT (see `LICENSE`).

This package concerns Conjectures 3.1 and 3.2 of Glover, Jones, Kempton and Oveson,
*Effects of backtracking on PageRank* (Appl. Math. Comput. 531 (2026) 130182;
arXiv:2211.13353v2). The paper is `paper/paper.pdf` (17 pages); its source is `paper/paper.tex`.

## Main results (status in brackets)

- House graph (5 vertices): exact cubic roof response; monotone iff ε ≤ ε_c ≈ 0.2406813;
  one strict minimum otherwise, at μ = 1 iff ε = 4/9 [proved].
- The endpoint shortcut (π_∞ vs π_1 predicts the side of π_0) fails at the roof iff ε > ε_0 ≈ 0.3291038 [proved].
- Minimum degree ≥ 2: order 5 is minimal for every ε ∈ (0,1); house for ε > ε_c, gem for ε < ε_g ≈ 0.3565988 [proved].
- Leaves allowed: the paw (4 vertices) fails iff ε ∈ (0, ε_1) ∪ (ε_2, 1), ε_1 ≈ 0.0783255, ε_2 ≈ 0.4355202;
  orders ≤ 3 are μ-invariant; the smallest order is 4 for ε ∉ [ε_1, ε_2] and 5 on [ε_1, ε_2] [proved].
- A five-vertex tree has a nonmonotone leaf on μ ∈ (0, 1) at ε = 1/2 [exact certificate].
- Certified census of all 995 connected graphs on 2–7 vertices at ε ∈ {1/10, 1/4, 1/3, 1/2, 17/20}:
  e.g. 353 of 507 seven-vertex graphs of minimum degree 2 have a nonmonotone vertex for some listed ε [certified, exact].
- Small-damping two-step criterion and an exact formula for dπ/dμ at μ = 1 with three endpoint tests [proved].
- Both readings of the spread clause of Conjecture 3.1 fail, also with unique extreme vertices [exact certificates].
- Exact invariance iff P²v = v [proved]; degree-tail Conjecture 3.2 fails for every tail index α > 1 [proved].
- Prevalence on standard networks: 10–41% of vertices nonmonotone, small excursions, stable rankings [numerical].

## Reproduce

Requires Python ≥ 3.11 with `python-flint`, `sympy`, `networkx` (tested: python-flint 0.9.0, sympy 1.14.0,
networkx 3.6.1); the prevalence script and figure also need `numpy`, `scipy`, `matplotlib`.

```bash
python3 -B -O code/verify_all.py OUTDIR --jobs 4     # all exact certificates + census (~1 min with 4 jobs)
python3 -B -O code/prevalence.py OUTDIR              # numerical prevalence study (~1 min)
python3 -B code/make_figure.py OUTDIR                # figure (illustration only)
cd paper && pdflatex paper.tex && pdflatex paper.tex && pdflatex paper.tex
```

The scripts write only into `OUTDIR`; JSON is written with sorted keys and no timestamps or paths.
`verify_all.py` uses explicit `RuntimeError` checks, so `-O` disables nothing. Running it with
`--jobs 1` and `--jobs 4` gives byte-identical output.

## Files

- `code/bt_exact.py`: exact solver, interpolation, squarefree decomposition and exact Sturm root counts on open intervals
  (roots at 0 and 1 divided out first), with known-root regressions.
- `code/census.py`: certified census; `code/symbolic_checks.py`: P3, paw, gem and diamond over Q(μ, ε).
- `code/criteria_checks.py`: slope formula, endpoint and two-step tests, spread (including the unique-extrema witness),
  the tree counterexample, rank reversal and the finite exact rank search, degree tails.
- `code/legacy_checks.py`: house, diamond, cover and clique-chain identities.
- `results/verification.json`, `results/census_full.json` (exact), `results/prevalence.json` (numerical).

## Limits

The census is certified at five damping values only. Covers, the period-class argument and the asymptotic
statements are proved in prose; no proof assistant was used. The journal version of the source was not
accessible, so the conjecture wording was checked against arXiv v2 only. The configuration-model tail
transfer is stated as an open problem.

## Release records

- [STATUS.md](STATUS.md): what is proved, certified or numerical, and what is open.
- [ASSURANCE.md](ASSURANCE.md): what has and has not been checked.
- [CLAIMS.json](CLAIMS.json): machine-readable claims, with status and paper location.
- [AI_INDEX.md](AI_INDEX.md): a navigation page for agents.
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): substantive objections from the internal AI-assisted reviews, and their resolution.
- [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md): replay commands and results. [ENVIRONMENT.txt](ENVIRONMENT.txt) and [requirements.txt](requirements.txt): the environment and pinned packages.
- [negative_controls.py](negative_controls.py): semantic negative controls. Run `python -B -O negative_controls.py`; it should reject corrupted claims, code and outputs.
- [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [RESEARCH_GATES.json](RESEARCH_GATES.json): sources, citation checks and the bounded novelty search.
- [PROVENANCE.md](PROVENANCE.md), [LICENSES.md](LICENSES.md), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md): provenance and the licence map.
- `SHA256SUMS`: hashes of every file. Check with `sha256sum -c SHA256SUMS`.

## Citation

Anonymous (2026). *Backtracking-weighted PageRank is not coordinatewise monotone: minimal counterexamples, a certified
census and endpoint tests* (Version 1.0.0-candidate) [Unrefereed research candidate]. Evidence Press.
https://doi.org/10.5281/zenodo.23580381

Machine-readable metadata: [CITATION.cff](CITATION.cff) and [.zenodo.json](.zenodo.json). Scholarly creator: Anonymous.
Ian Pitchford maintains the repository and publishes it; that role is not authorship.
