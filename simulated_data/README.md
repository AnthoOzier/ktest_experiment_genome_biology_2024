Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper 
Section 2.4 "Kernel testing is calibrated and powerful on simulated data"

Install the adapted version of the ktest package 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology#subdirectory=python
```

### Data generation : 
- `data_generation.R` to generate the data (this code is the property of the authors of [**Distribution-free complex hypothesis testing for single-cell RNA-seq differential expression analysis** - 
Marine Gauthier, Denis Agniel, Rodolphe Thiébaut, Boris P. Hejblum
 ](https://doi.org/10.1101/2021.05.21.445165))

### Data analysis : 
- We propose a script separated in three files `data_analysis0.sh` `data_analysis1.sh` and `data_analysis.py` to parallelize the data analysis on a computer cluster based on SLURM. The `data_analysis0.sh` file sends requests to the SLURM job scheduler as `data_analysis1.sh` files with specified parameters. Then the `data_analysis1.sh` file launch the `data_analysis.py` file that performs the ktest data analysis. 

### Post analysis : 
- `post_analysis.py` computes the Benjamini-Hochberg correction of the p-values and aggregates the results in one file. 
- `figure_generation.ipynb` is a Jupiter Notebook that generates the figures of the paper. 

### Others : 
- `get_params.py` is imported in `data_analysis.py` to load the parameters specified in `data_analysis1.sh` 