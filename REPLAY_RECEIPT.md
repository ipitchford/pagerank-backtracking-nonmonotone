# Replay receipt

## Commands

From the root of a fresh extraction, using Python 3.13 with the pins in [requirements.txt](requirements.txt):

```sh
sha256sum -c SHA256SUMS                                   # macOS: shasum -a 256 -c SHA256SUMS
python -B -O code/verify_all.py OUT --jobs 4              # exact certificates + census
cmp OUT/verification.json results/verification.json
cmp OUT/census_full.json  results/census_full.json
python -B -O code/prevalence.py OUT                       # numerical study
cmp OUT/prevalence.json   results/prevalence.json         # same host and pins; floating point
python -B -O negative_controls.py                         # semantic negative controls
```

`verify_all.py` and `prevalence.py` write only into `OUT` and refuse to write beside themselves. `negative_controls.py` works only in temporary directories.

## Runs before publication

| Run | Command | Result |
|---|---|---|
| Clean replay of the revised package before the final fix round, 22 October 2026, CPython 3.13.5, macOS arm64 | `code/verify_all.py OUT` and `code/prevalence.py OUT` (the retained log records the output, not the interpreter flags) | Status PASS; byte-identical to the results of that version. Superseded for the exact results by the runs below; `prevalence.json` and `prevalence.py` did not change in the final fix round |
| Producer replays of the final package, final fix round, 22–23 October 2026 | `python -B -O code/verify_all.py OUT` | Three replays byte-identical, according to the final assurance record |
| Release preparation, final package, 23 October 2026, CPython 3.13.5, macOS arm64, from a fresh copy | `python -B -O code/verify_all.py OUT --jobs 2` (about 300 s) | Status PASS. `verification.json` and `census_full.json` byte-identical to `results/` |
| Negative controls, release preparation, final package, 23 October 2026 | `python -B -O negative_controls.py` (about 40–60 s) | All 17 controls behaved. The baseline passed. Rejected: the corrupted certificate; six flipped claims (including the Proposition 5.5 tree and the pre-repair Atlas 399 value); the wrong damping value; five code mutations (two in the Sturm root counter); two tampered outputs; and a tampered file |

## Release-payload replay

This file cannot record a hash of the archive that contains it, so the replay below was run on the assembled release tree. That tree is the final package plus these release records and the DOI-bearing PDF. Its code and results files are the ones listed in `SHA256SUMS`.

| Item | Value |
|---|---|
| Date and platform | 2026-10-22T23:38Z, CPython 3.13.5, macOS arm64 (Darwin) |
| `python -B -O code/verify_all.py OUT --jobs 4` | exit 0; `verification.json` and `census_full.json` byte-identical to `results/` |
| `python -B code/verify_all.py OUT --jobs 4` (assertions enabled) | exit 0; both files byte-identical to `results/` |
| `python -B -O code/prevalence.py OUT` | exit 0; `prevalence.json` byte-identical to `results/` (same platform and pins) |
| `python -B -O negative_controls.py` | ALL NEGATIVE CONTROLS BEHAVED |

The frozen archive was replayed again at publication, from a fresh extraction (`-B -O`, both exact result files byte-identical, negative controls). That run is part of the publisher's gate record. The repository's CI workflow repeats the SHA-256 check, the exact replay and the negative controls on GitHub for the tagged commit.

## What this establishes

Byte-identical replay shows that the shipped outputs are exactly what the shipped programs produce on the stated environment. It does not show that the programs implement the manuscript's definitions correctly. That rests on three things: the method description in Sections 7 and 12; the negative controls; and the independent exact recomputations in the internal reviews (a full 12-arc symbolic solve of the house, and the re-review's own exact solver). None of these runs is an independent reproduction by an unaffiliated party. The prevalence comparison is floating point. Byte identity of `prevalence.json` is expected only on the same platform and pins.
