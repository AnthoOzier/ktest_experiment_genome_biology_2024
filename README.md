# ktest_experiment_genome_biology_2024

Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper "ktest : Kernel-Based Testing for Single-Cell Differential Analysis"

## Data availability 

- The reversion and scChIP-Seq datasets can be loaded here : https://zenodo.org/records/11001869
- The experimental scRNAseq datasets associated to [Squair, J.W., Gautier, M., Kathe, C. et al. Confronting false discoveries in single-cell differential expression. Nat Commun 12, 5692 (2021)](https://doi.org/10.1038/s41467-021-25960-2) can be loaded here : https://zenodo.org/records/5048449

## Folders 

- `simulated_data/` analyses corresponding to Section 2.4 "Kernel testing is calibrated and powerful on simulated data"
- `experimental_scrnaseq_data/` analyses corresponding to Section 2.5 "Challenging DEA metods on experimental scRNA-Seq data"
- `reversion/` analyses corresponding to Section 2.6 "Kernel testing reveals the heterogeneity of reverting cells"
- `scchipseq/` analyses corresponding to Section 3 "Towards a new testing framework for differential binding analysis in single-cell ChIP-Seq data"


## ktest versions 

The analyses are performed with different versions of the `ktest` package. 

- For scripts in `simulated_data/` and `experimental_scrnaseq_data/` :
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology#subdirectory=python
```
- For scipts in `reversion/` and `scchipseq/` : 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology_reversion#subdirectory=python
```
