# E. coli MIC prediction: curated resistance features versus Evo 2

This repository contains the completed comparison of curated AMRFinderPlus
features and frozen Evo 2 resistance-locus embeddings for antibiotic MIC
prediction in *Escherichia coli*.

Curated allele, gene-family, and drug-class features outperformed the tested
Evo 2 representation in controlled variant-withholding experiments. Adding Evo 2
did not demonstrate improvement in the original locked external evaluation.
The external family-feature follow-up is exploratory.

- [Project overview and reproduction instructions](amr-ecoli-project/README.md)
- [Paper PDF](amr-ecoli-project/paper/manuscript.pdf)
- [Supplement](amr-ecoli-project/paper/supplement.pdf)
- [Code, data, embeddings, and submission source release](https://github.com/Sekiromike/AMR-Resistance-Gene-Detection/releases/tag/v1.0.0)
- [Current study protocol](amr-ecoli-project/docs/STUDY_PROTOCOL.md)
- [Verified results](amr-ecoli-project/paper/verification/audit.json)

**Author:** Rushath Rajeev, independent researcher. **Computing:** UMass Unity
HPC. Claude Code and Codex are acknowledged as AI assistance.

Commands in the project documentation run from `amr-ecoli-project/`.
Original code is [MIT licensed](amr-ecoli-project/LICENSE); upstream materials
retain their [own terms](amr-ecoli-project/THIRD_PARTY_NOTICES.md).
