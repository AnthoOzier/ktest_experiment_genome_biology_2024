Use this code to reproduce the analyses of the Ozier-Lafontaine et al. 2024 Genome Biology paper Section 2.5 "Challenging DEA metods on experimental scRNA-Seq data"

Install the adapted version of the ktest package 
```
pip  install ktest@git+https://github.com/LMJL-Alea/ktest@publication_genome_biology#subdirectory=python
```

Data preparation : 
- Load the data on https://zenodo.org/records/7659806
- `data_preparation.R` to normalize the data and run the other DE tests. 

Data analysis : 
- We propose a script separated in three files `data_analysis0.sh` `data_analysis1.sh` and `data_analysis.py` to parallelize the data analysis on a computer cluster based on SLURM. The `data_analysis0.sh` file sends requests to the SLURM job scheduler as `data_analysis1.sh` files with specified parameters. Then the `data_analysis1.sh` file launch the `data_analysis.py` file that performs the ktest data analysis. 

Post analysis : 
- `results_aggregation.py` to compute the Benjamini-Hochberg correction of the p-values and aggregate the ktest results in one file per dataset.
- `AUCC_computation.R` to compute the AUCC scores based on functions in `utils_aucc.R`. 
- `figure_generation.ipynb` is a Jupiter Notebook that generates the figures of the paper. 

Others : 
- `get_params.py` is imported in `data_analysis.py` to load the parameters specified in `data_analysis1.sh` 