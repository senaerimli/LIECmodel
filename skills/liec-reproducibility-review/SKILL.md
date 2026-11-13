---
name: liec-reproducibility-review
description: Review the LIEC GitHub repository, MATLAB/COBRA scripts, or a proposed upload for reproducibility and publication readiness. Use whenever the user shares code, a folder tree, a repo link, or asks about GitHub, file organisation, versions, solvers, large files, licenses, or MEMOTE.
---

# Reproducibility review

Check, and report as a short checklist with fixes ordered by priority:
- Hard-coded absolute paths (e.g., C:\Users\...); replace with relative paths.
- Files over 50 MB (warning) or 100 MB (rejected by GitHub); move to Zenodo.
- Scripts run in part-number order from a clean clone; outputs regenerate figures and tables.
- Versions recorded in environment.md (MATLAB, COBRA commit, each solver used per result).
- Commercial-solver dependence: say which results need it, and test with a free solver.
- Random seeds fixed where relevant.
- Raw GEO data downloaded by script, not stored.
- Models exported in SBML plus .mat; MEMOTE report in docs/.
- No license files, credentials, personal data, or private paths committed.
- README, LICENSE, CITATION.cff, release tag, and Zenodo DOI present before submission.
- Credits and license compatibility for reused task lists and models (Harvey/Harvetta, iColon, etc.).
Explain Git steps for a Windows beginner using GitHub Desktop, in plain language.
