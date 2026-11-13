---
name: liec-evidence-check
description: Keep every scientific statement in the LIEC / ulcerative colitis genome-scale modeling project evidence-backed. Use whenever the user asks for scientific explanations, manuscript text, literature claims, method justifications, dataset descriptions, or interpretation of model results, even if they do not ask for references.
---

# Evidence check for LIEC research

1. Separate claims into: (a) established in literature, (b) shown by the user's own results, (c) hypothesis or interpretation. Label them.
2. Never invent references, DOIs, GEO details, or numbers. If unsure, say so and search the web or ask the user to supply the paper.
3. For dataset facts (platform, sample groups, timepoints, response definitions), verify against the GEO record or the original paper before using them.
4. For claims about tools or algorithms (COBRA, GIMME, E-Flux, MOMA, FVA, dFBA, RAVEN, MEMOTE), check the original paper or current documentation.
5. For results, point to the exact script/table/figure in the repository. If the result cannot be traced to a script, flag it as a reproducibility risk.
6. In manuscript drafting, mark any statement that lacks a citation with [CITATION NEEDED] instead of filling it in from memory.
7. State limitations plainly (e.g., microarray vs RNA-seq, cross-platform batch effects, bound artefacts) rather than smoothing them over.
8. Keep quotations short and paraphrase sources.
