Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper 
Section 3 "Towards a new testing framework for differential binding analysis in single-cell ChIP-Seq data"

**Install the adapted version of the ktest package** 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology_reversion#subdirectory=python
```

**Load the datasets** : https://zenodo.org/records/11001869

1/ init.R to load the package and use it in R
2/ analysis_persistent_vs_untreated.R to perform differential analysis based on kernels
3/ get_sub_populations.R to get the subpopulations of the untreated cells and to identify the corresponding cells
4/ analysis_untreated_x_vs_untreated_y.R with x=1,2,3, y = 1,2,3 to perform pairwise comparisons between subpopulations of the untreated population
