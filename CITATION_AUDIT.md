# Citation audit

This file states only what the retained review and revision records support. All checks were made inside the Evidence Press workflow by AI agents. They are bibliographic checks, not a judgement of the cited works.

## Identifier resolution

| Check | Date | Result | Record |
|---|---|---|---|
| The source's journal version, doi:10.1016/j.amc.2026.130182 | 22 Oct 2026 | Crossref returns *Effects of backtracking on PageRank*, Applied Mathematics and Computation 531, article 130182 (2026), with authors Cory Glover, Tyler Jones, Mark Kempton and Alice Oveson. The full text returned HTTP 403, so the conjecture wording was not compared with the journal version. | First-round domain review; re-review; publication check on 22 Oct 2026 |
| Six DOIs of the first draft's direct predecessors | 22 Oct 2026 | All resolved via Crossref to the stated title, journal, volume and authors | First-round domain review |
| All 17 DOIs in the revised reference list | 22 Oct 2026 | All resolved via Crossref. Titles, venues, volumes, pages and authors match the bibliography. | Verification re-review |
| arXiv:2211.13353v2 and arXiv:2604.13681 | 22 Oct 2026 | Resolved on arxiv.org. Title and authors match. The v2 PDF carries the 9 February 2026 stamp. | Verification re-review |
| Kemeny & Snell, *Finite Markov Chains* (1960; Springer reprint 1976, ISBN 978-0-387-90192-3) | 22 Oct 2026 | Has no DOI. Checked against the Open Library catalogue. | Revision record |
| Knuth, *The Stanford GraphBase* (1993); Hagberg, Schult & Swart, NetworkX, Proc. SciPy 2008, 11–15 | 22–23 Oct 2026 | Added in the final fix round as data and software sources. Neither entry has a DOI. No separate catalogue check is recorded. | Final fix round |

Corrections made during revision: the Arrigo et al. issue number and the Henning et al. page range. Zhang & Wu is dated by its 2025 volume year; it was online on 16 November 2024.

## Source-to-claim checks

* Conjectures 3.1 and 3.2 are quoted verbatim from arXiv:2211.13353v2 (p. 10 and p. 11; §3.3 and §3.4). The re-review checked them word for word against an independently fetched copy, and the adversarial audit inspected them directly.
* Model fidelity: the kernel, edge restart u = Tᵀ(TTᵀ)⁻¹v and tail projection agree with the source's Definition 2.4. The π_∞ formula is the source's Proposition 3.3. Both the adversarial audit and the re-review checked this. The source's transpose inconsistency in some proofs is noted in Remark 2.1.
* Attributions: the re-review checked these against the cited works' content —
  * Boldi et al. (fibrations) and Kemeny–Snell (lumpability) for Lemma 5.6;
  * Henning et al. Section 3 for the mechanism of Theorem 11.2;
  * Arrigo et al. Theorem 4 for invariance at μ = 0;
  * the context citations.
* Aleja et al.'s leaf conventions could not be rechecked at re-review because of a paywall. The full published paper was read during the 20 October novelty search.

## Not done

There was no systematic citation export from MathSciNet, zbMATH, Scopus or Web of Science, no exhaustive search in languages other than English, and no enquiry to the source authors.
