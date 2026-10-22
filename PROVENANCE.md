# Provenance

* Scholarly creator: **Anonymous**, the standing Evidence Press attribution. No human research contribution is claimed or inferred. There is no ORCID and no affiliation.
* Publisher: Evidence Press (https://evidencepress.org; release page https://evidencepress.org/releases/pagerank-backtracking-nonmonotone/). Ian Pitchford maintains the repository and carries out publication operations under the repository owner's account. That role is not authorship.
* Production: AI-assisted research. AI research agents produced the manuscript, proofs and computations. Errors were sought in four ways: an internal multi-perspective AI-assisted review of the first draft; an AI-assisted adversarial audit of that draft by a different model (configured as OpenAI gpt-6-astra, 22 October 2026); a verification re-review of the revised version, which recomputed the certificates with its own exact solver; and a second adversarial audit of the revised version by the same external model. That second audit found a census root-counting error and two new counterexamples (a five-vertex tree, and a unique-extrema spread witness). A final fix round followed. Their substantive points and how each was resolved are in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). The original reports are kept privately.
* Target: Glover, Jones, Kempton and Oveson, *Effects of backtracking on PageRank*, Applied Mathematics and Computation 531 (2026) Article 130182, doi:10.1016/j.amc.2026.130182. Open version: arXiv:2211.13353v2 (revised 9 February 2026), whose numbering is used throughout. Section 1.1 of the manuscript quotes Conjectures 3.1 and 3.2 verbatim from arXiv v2, with attribution.
* History:
  * An 11-page first draft, dated 20 October 2026, had a bounded novelty search run on the same day.
  * The draft was reviewed and audited on 22 October 2026.
  * It was then revised and extended the same day: retitled, with a paw classification, a gem proposition, a certified census, a criterion with endpoint tests, the spread clause, a prevalence study, and a degree-tail theorem for every tail index.
  * The revised package was replayed from a clean copy, re-reviewed and audited again.
  * In a final fix round the census root counter was replaced by exact Sturm counts and the census regenerated. The paper gained Corollary 5.4, Proposition 5.5 and Proposition 8.2, and the re-review's textual corrections were made. The package was then replayed again (17 pages).
  * Only the revised version is published.
* Publication edits:
  * the reserved DOI and repository URL inserted in the manuscript front matter, and the PDF rebuilt (see the publication notes);
  * release identity added to README.md;
  * the release records listed in README.md added, together with `negative_controls.py` and `.github/workflows/replay.yml`;
  * SHA256SUMS generated.
  No computation script or results file was changed by publication.
