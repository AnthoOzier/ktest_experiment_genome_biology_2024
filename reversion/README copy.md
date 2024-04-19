Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper 
Section 2.6 "Kernel testing reveals the heterogeneity of reverting cells"


Install the ktest package 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology_reversion#subdirectory=python
```
Files description : 
- `reversion_RTqPCR.ipynb` analyses of the RTsPCR data 
- `reversion_scRNAseq.ipynb` analyses of the scRNAseq data
- `utils_reversion.py` functions of interest for the analyses
- `data/RTqPCR_reversion_logx.csv` : log(x+1) normalized RTqPCR data
- `data/RTqPCR_reversion_logcentered.csv` : log(x+1) normalized RTqPCR data centered by batch 
- `data/reversion_SCT_residuals.csv`: Person residuals of the scRNAseq data obtained through the sctransform R package.
- `data/reversion_SCT_residuals_batch_corrected.csv`: Person residuals of the scRNAseq data obtained through the sctransform R package centered by batch.

