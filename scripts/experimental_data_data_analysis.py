
import pandas as pd
from time import time
import os 
from ktest.tester import Ktest
import numpy as np
from get_params import get_params


if __name__ == '__main__':

    key = get_params()
    nystrom = key['nystrom']
    n = key['n']
    n_jobs = key['n_jobs']
    
    s0 = key['s0']
    s1 = key['s1']
    seeds= list(range(s0,s1)) if s1-s0>1 else [s0]
    stat = key['stat']
    diroutput = key['diroutput']
    dirdata =  key['dirdata']
    filename = key['filename']
    lots = key['lots']
    permutation = key['permutation']
    n_permutations = key['n_permutations']
    ignore_zeros = key['ignore_zeros'] 
    kernel_initial = key['kernel']
        

    if stat == 'kfda':
        if (nystrom and n>40) or not nystrom:
            for seed in seeds:
                print(f'n={n} s={seed}',end=' ')
                name=f'n{n}_s{seed}'
                
                dfx = pd.read_csv(f'{dirdata}X{n}_s{seed}.csv',index_col=0)
                dfy = pd.read_csv(f'{dirdata}Y{n}_s{seed}.csv',index_col=0).T
                dfx.index = dfy.index
                t = Ktest(dfy,dfx,condition='x',data_name=name,nystrom=nystrom,
                          permutation=permutation,n_permutations=n_permutations)
                if kernel_initial=='linear':
                    kernel={'function':'linear'}
                    kernel_info = None    
                elif kernel_initial == 'fisher_zero_inflated_gaussian':
                    kernel = kernel_initial
                    propZeros = 1 - pd.read_csv(f'{dirdata}propNonZeros{n}_s{seed}.csv',index_col=0)
                    propZeros.columns = ['propZeros']
                    propZeros.index = t.get_var().index
                    t.update_var_from_dataframe(propZeros)
                    kernel_info = ['propZeros','propZeros']
                elif kernel_initial == 'naive_fisher_zero_inflated_gaussian':
                    t.add_zero_proportions_to_var()
                    kernel = 'fisher_zero_inflated_gaussian'
                    kernel_info = [f'{k}_pz' for k in t.get_samples_list()]  
                elif kernel_initial == 'gauss_quantile':
                    quantile = key['quantile']/100
                    kernel = {'function':'gauss','bandwidth':'quantile','median_coef':quantile,'kernel_name':f'gq_{quantile}'}
                    kernel_info=None
                else:
                    kernel={'function': 'gauss', 'bandwidth': 'median', 'median_coef': 1, 'kernel_name': None}
                    kernel_info = None    
                t.univariate_test(n_jobs=n_jobs,lots=lots,save_path=diroutput,
                                  kernel=kernel,kernel_info=kernel_info,
                                  ignore_zeros=ignore_zeros,verbose=1)
     