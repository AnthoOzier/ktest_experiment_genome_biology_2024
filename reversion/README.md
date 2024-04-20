Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper 
Section 2.6 "Kernel testing reveals the heterogeneity of reverting cells"


**Install the adapted version of the ktest package** 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology_reversion#subdirectory=python
```

**Load the datasets** : https://zenodo.org/records/11001869

### Notebooks : 
- `reversion_RTqPCR.ipynb` analyses of the RTsPCR data 
- `reversion_scRNAseq.ipynb` analyses of the scRNAseq data

### Scripts : 
- `utils_reversion.py` functions of interest for the analyses

### DE-tables : 
- `.csv` files containing the asymptotic p-values associated to each gene for each pairwise comparison
