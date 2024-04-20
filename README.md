# ktest_experiment_genome_biology_2024

Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper "ktest : Kernel-Based Testing for Single-Cell Differential Analysis"

The datasets can be loaded here : https://zenodo.org/records/11001869

For the analysis corresponding to Section 2.4 "Kernel testing is calibrated and powerful on simulated data", see `simulated_data/`

For the analysis corresponding to Section 2.5 "Challenging DEA metods on experimental scRNA-Seq data", see `experimental_scrnaseq_data/`

For the analysis corresponding to Section 2.6 "Kernel testing reveals the heterogeneity of reverting cells", see `reversion/`

For the analysis corresponding to Section 3 "Towards a new testing framework for differential binding analysis in single-cell ChIP-Seq data", see `scchipseq/`

Scripts corresponding to Sections 2.4 and 2.5 are based on the branch [publication_genome_biology](https://github.com/LMJL-Alea/ktest/tree/publication_genome_biology) of the ktest package : 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology#subdirectory=python
```

Scripts corresponding to Sections 2.6 and 3 are based on the branch[publication_genome_biology_reversion](https://github.com/LMJL-Alea/ktest/tree/publication_genome_biology_reversion) of the ktest package :
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology_reversion#subdirectory=python
```
