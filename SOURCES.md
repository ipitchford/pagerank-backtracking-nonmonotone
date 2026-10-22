# Sources and their roles

The manuscript's reference list has 22 entries: 17 with DOIs, 2 arXiv-only, 2 books and 1 conference paper without a DOI. Reading depth is recorded below. It follows the manuscript (Section 12), the novelty search of 20 October 2026, and the internal review records. Bibliographic checks are in [CITATION_AUDIT.md](CITATION_AUDIT.md).

| Source | Role in this work | Reading |
|---|---|---|
| Glover, Jones, Kempton & Oveson, *Effects of backtracking on PageRank*, arXiv:2211.13353v2 (9 Feb 2026); journal version Appl. Math. Comput. 531 (2026) 130182, doi:10.1016/j.amc.2026.130182 | The target. It supplies Conjectures 3.1 and 3.2 (quoted verbatim in §1.1), the model (Definition 2.4), the sufficient invariance cases (Theorems 3.1–3.2) and the π_∞ formula (Proposition 3.3) | arXiv v2 read in full. The journal text was inaccessible (HTTP 403), so numbering and wording follow arXiv v2 |
| Henning, van der Hofstad & Litvak, *Power-law hypothesis for PageRank on undirected graphs*, Electron. J. Probab. 30 (2025), doi:10.1214/25-EJP1415; arXiv:2407.13730 | The assortative-component mechanism adapted in Theorem 11.2, and the local-weak-convergence reading of Conjecture 3.2 | Section 3 (Theorem 3.6), arXiv version |
| Aleja, Criado, García del Amo, Pérez & Romance (2019), doi:10.1016/j.chaos.2019.06.017 | The same kernel on μ ∈ [0, 1] (Definitions 3.5, 3.7); different leaf conventions | Full published paper |
| Arrigo, Higham & Noferini (2019), doi:10.1007/s10915-019-00981-8 | Sufficiency of invariance at μ = 0 (Theorem 4) | Full paper |
| Criado, Flores, García, García del Amo, Pérez & Romance (2019), doi:10.1016/j.cam.2018.09.048 | A related return-weighted centrality (unnormalised Perron vector) | Sections 2–4 of the author preprint |
| Meng & Masuda (2020), doi:10.1098/rspa.2020.0447; Grover & Leskovec (2016), doi:10.1145/2939672.2939754 | node2vec kernel translation; leaf convention; stationary law without restart | Meng–Masuda relevant sections; node2vec for the kernel |
| Avena, Bet, Schroeder & Stegehuis, arXiv:2604.13681 (2026) | A later paper citing Glover et al.; the unrestarted α, β, γ walk | Full v1 (Definition 2.1, Theorem 3.6) |
| Boldi, Lonati, Santini & Vigna (2006), doi:10.1051/ita:2006004; Kemeny & Snell (1960) | PageRank invariance under fibrations and lumpability: the basis of Lemma 5.6 | Cited for the method |
| Garavaglia, van der Hofstad & Litvak (2020), doi:10.1214/19-AAP1494 | Local weak convergence of PageRank | Cited for context |
| Zhang & Wu (2025), doi:10.1007/s10915-024-02735-7 | Recent computational work on non-backtracking PageRank | Abstract and metadata |
| Knuth, *The Stanford GraphBase* (1993); Hagberg, Schult & Swart, NetworkX (SciPy 2008) | Data source of the Les Misérables co-occurrence graph (Proposition 10.1) and the graph atlas and test graphs used by the code (networkx 3.6.1) | Cited as data and software sources |
| Boldi, Santini & Vigna (2009); Gleich (2015); Gleich, Lim & Yu (2015); Benson, Gleich & Lim (2017); Hashimoto (1989); Alon et al. (2007); Kempton (2016) | Context: dependence on the damping parameter, the PageRank survey, higher-order (nonlinear) PageRank, non-backtracking walks and edge lifts | Cited for context; metadata verified |

The manuscript's literature boundary (Section 12): within the search, no prior counterexample to Conjecture 3.1 and no census was found. No exhaustive priority claim is made. See also [RESEARCH_GATES.json](RESEARCH_GATES.json).
