### version ccipl 


import pandas as pd 
import os
from ktest.pvalues import correct_BenjaminiHochberg_pval_of_dfcolumn
from ktest.utils_simu_ccdf import post_traitement_ccdf_of_column
from itertools import product
from joblib import parallel_backend
from joblib import Parallel, delayed


data="2023_06_ccdf"
log= False
log_str='_lognormalized' if log else ''
data_str=f'{data}{log_str}'


stat="kfda"
kernel='naive_fisher_zero_inflated_gaussian'
kernel='linear'
kernel= 'gauss_kernel_q80'
ignore_zeros=False
permutation=False
nystrom=False


if kernel=='linear':
    kernel_str='_klin_'
elif kernel == 'fisher_zero_inflated_gaussian':
    kernel_str='_kfzig_'
elif kernel =='naive_fisher_zero_inflated_gaussian':
    kernel_str='_naivekfzig_'
elif kernel == 'gauss_kernel_q80':
    kernel_str='_kgq80_'
else:
    kernel_str=''
perm_str = '_perm_' if permutation else ''
ny_str = '_ny_' if nystrom else ''
ignore_zeros_str = '_iz_' if ignore_zeros else ''

model_str=f'{stat}{kernel_str}{ny_str}{perm_str}{ignore_zeros_str}'

path = f'/home/data/experimental_data/{data_str}_{model_str}/'
path_output = f'/home/data/experimental_data/scores_{data_str}_{model_str}/'
tmax = 30
print(path)
print(path_output)
seeds = range(501)
ns = [20,40,60,80,100,160,200,500]
n_jobs = 64


def assert_output_exists(path_output,file_output):    
    output_exists = os.path.isfile(path_output+file_output)
    if output_exists:
        df = pd.read_csv(path_output+file_output,index_col=0)
        if df.empty:
            output_exists=False
    return(output_exists)

def post_analysis(s,n,ny,path,path_output):
    print(s,n,ny)
    if (ny and n>40) or not ny:
        name = f'_ny_lmrandom_m{int(n/10)}_basisw_n{n}_s{s}x' if ny else f'_standard_n{n}_s{s}x'
        file = name+'_univariate.csv'
        file_output= name+f'_scores.csv'

        input_exists = os.path.isfile(path+file) and os.path.getsize(path+file)!=0
        output_exists = assert_output_exists(path_output,file_output)
        print(file,input_exists,output_exists)
        if input_exists and not output_exists:
            df = pd.read_csv(path+file,index_col=0) 
            print(df.columns)
            if len(df) == 10000:
                ts = []
                for t in range(1,tmax):
                    # print(f't={t}')
                    if f'{name}_t{t}_pval' in df:
                        pval = df[f'{name}_t{t}_pval']
                        print('na pval',pval.isna().sum())
                        if pval.isna().sum()==0:
                            ts += [t]
                for t in ts:
                    pval = df[f'{name}_t{t}_pval']
                    pvalc = correct_BenjaminiHochberg_pval_of_dfcolumn(pval)
                    df[f'{name}_t{t}_pvalBH'] = pvalc
                res = []        
                for t in ts:
                    dict_hyperparams = {'m':'ny' if ny else 'st','n':n,'s':s,'k':'g','t':t}
                    pval = df[f'{name}_t{t}_pval']
                    pvalc = df[f'{name}_t{t}_pvalBH']
                    res += post_traitement_ccdf_of_column(pval,pvalc,dict_hyperparams)
                dfres = pd.DataFrame(res)
                dfres.to_csv(path_output+file_output)
                print(f'{file_output} saved')
if n_jobs >1:
    with parallel_backend('loky'):
        a = Parallel(n_jobs=n_jobs)(delayed(post_analysis)(s,n,ny,path,path_output) \
                                    for s,n,ny in product(seeds[::-1],ns,[True,False]))
else:
    for s,n,ny in product(seeds[::-1],ns,[True,False]):
        post_analysis(s,n,ny,path,path_output)

        
