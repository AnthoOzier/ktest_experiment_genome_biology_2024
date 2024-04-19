from ktest.tester import Ktest
import pandas as pd
import os
from get_params import get_params

def center_in_input_space(self,center_by=None,center_non_zero_only=False):
    if center_by == 'replicate':

        data_r = self.get_data(condition='replicate',dataframe=True,in_dict=True)
        all_data_ = self.get_data(condition='replicate',dataframe=True,in_dict=False)
        dfc = all_data_[all_data_!=0].mean() if center_non_zero_only else all_data_.mean()
        
        data_ ={}
        for k,v in data_r.items():
#             print(k)
            if center_non_zero_only:
                v[v!=0] = v[v!=0] - v[v!=0].mean() + dfc
                data_[k] = v
            else:
                data_[k] = v - v.mean() + dfc


    elif center_by == '#-replicate_+label':    
        data_r = self.get_data(condition='replicate',dataframe=True)
        data_l = self.get_data(condition='label',dataframe=True)
        if center_non_zero_only:
            mean_l = {k:v[v!=0].mean() for k,v in data_l.items()}
            mean_r = {k:v[v!=0].mean() for k,v in data_r.items()}
        else:
            mean_l = {k:v.mean() for k,v in data_l.items()}
            mean_r = {k:v.mean() for k,v in data_r.items()}

        metadata_l = self.get_metadata(condition='label',in_dict=True)
        r_in_l = {k:metadata_l[k]['replicate'].unique().to_list() for k in metadata_l.keys()}

        data_ = {}
        for l,ml in mean_l.items():
            for r,dr in data_r.items():
                if r in r_in_l[l]:
                    if center_non_zero_only:
                        dr[dr!=0] = dr[dr!=0] - mean_r[r] + ml
                        data_[r] = dr
                    else:
                        data_[r] =  dr - mean_r[r] + ml
                            
    new_df = pd.concat(data_.values())
    new_meta = self.get_metadata(condition='label',in_dict=False)
    new_meta = new_meta.iloc[new_meta.index.get_indexer(new_df.index)].copy()

    return(new_df,new_meta)



key = get_params()
n_jobs=key['n_jobs']
file=key['file']+'.csv'
diroutput=key['diroutput']
dirdata=key['dirdata']
filename=key['filename']
kernel_name=key['kernel'] if 'kernel' in key else 'gauss_median'
ignore_zeros=key['ignore_zeros'] 
binary=key['binary']
center_by=key['center_by']
multivariate=key['multivariate']
center_non_zero_only=key['center_non_zero_only']

print(key)
if 'v0' in key:
    v0 = key['v0']
    v1 = key['v1']
else:
    v0 = 0
    v1 = -1

lots = key['lots'] if 'lots' in key else 100

prefix = 'normalized' if 'normalized' in file else 'residuals'
dataset = file.replace(prefix,'').replace('.csv','')
data_type = 'sn_' if prefix == 'normalized' else ''
izstr = 'iz_' if ignore_zeros else ''
binstr = 'bin_' if binary else ''
# cbstr = f'cb{center_by}_' if center_by is not None else ''

if kernel_name=='linear':
    kernel={'function':'linear'}
    kernelstr = 'lin_' 
elif kernel_name=='gauss_quantile':
    quantile=key['quantile']/100
    kernelstr=f'kgq{quantile}'
    kernel={'function': 'gauss', 'bandwidth': 'quantile', 'median_coef':quantile, 'kernel_name': kernelstr}
    
else:
    kernel={'function': 'gauss', 'bandwidth': 'median', 'median_coef': 1, 'kernel_name': None}
    kernelstr = ''
data_name = f'{kernelstr}{binstr}{izstr}{data_type}{dataset}'


file_meta = file.replace(prefix,"meta")
df = pd.read_csv(f'{dirdata}{file}',index_col=0).T
if 'normalized' in file:
    df = df.T[df.mean()!=0].T
    df = df[df.columns[(df[df!=0]).count()>2]]
if binary:
    df = (df != 0).astype(int)
meta = pd.read_csv(f'{dirdata}{file_meta}',index_col=0)
print(file,data_name,meta.label.unique(),df.shape[0],' cells',df.shape[1],' variables')

print('df',df.shape,'meta',meta.shape)
nystrom = True if len(df)>4000 else False

center_by_in_input_space = None
if center_by is not None:
    if "#plus#" in center_by:
        center_by = center_by.replace('#plus#','+')
    if 'in_input_space' in center_by:
        center_by_in_input_space = center_by.replace('_in_input_space','')
        center_by = None

t = Ktest(
    data=df,
    metadata=meta.copy(),
    data_name=data_name,
    condition='label',
    nystrom=nystrom,
    kernel=kernel,
    center_by=center_by
    )

if center_by_in_input_space is not None:

    datac,metac = center_in_input_space(t,
                                        center_by=center_by_in_input_space,
                                        center_non_zero_only=center_non_zero_only)

    if center_by_in_input_space == 'replicate':
        data_name += f'cris_'
    elif center_by_in_input_space == '#-replicate_+label':
        data_name += f'crlis_'
    if center_non_zero_only :
        data_name+='cnz_'


    t = Ktest(data=datac,
                metadata=metac,
                condition='label',
                nystrom=nystrom,
                center_by=None,
                kernel=kernel,
                data_name=data_name)


if multivariate:
    t.multivariate_test()
    name = t.projections()
    t.df_proj_kfda[name].to_csv(f'{diroutput}{name}_proj_kfda.csv')
    t.df_proj_kpca[name].to_csv(f'{diroutput}{name}_proj_kpca.csv')
    t.get_diagnostics(diff=True,var_samples=True,var_within=True,kfdr=True).to_csv(f'{diroutput}{name}_diagnostics.csv')
    sp,ev=t.get_spev()
    name=t.get_kfdat_name()
    pd.DataFrame(sp,columns=['spectrum']).to_csv(f'{diroutput}{name}_sp.csv')
    pd.DataFrame(ev).to_csv(f'{diroutput}{name}_ev.csv')
    t.df_pval.to_csv(f'{diroutput}{name}_pval.csv')
    t.df_kfdat.to_csv(f'{diroutput}{name}_kfda.csv')

    t.add_zero_proportions_to_var()
    var = t.get_var()
    x = t.get_data(in_dict=False,dataframe=True)
    var['var'] = x.var()
    var['mean'] = x.mean()
    var.to_csv(f'{diroutput}{name}_var.csv')



else:
    if kernel_name == 'fisher_zero_inflated_gaussian':


        meta = t.get_metadata(samples='all',in_dict=False)
        df = t.get_data(samples='all',in_dict=False,dataframe=True)
        cd =  pd.concat([meta[['label','replicate']],df],axis=1)
        expr_per_replicate = cd[cd!=0].groupby(by=['label','replicate']).count()

        not_expressed_in_one_condition = expr_per_replicate.groupby(by='label').max().min() == 0
        expressed_less_than_once_per_replicate = expr_per_replicate.groupby(by='label').max().max() == 1

        # mask = (~(not_expressed_in_one_condition * expressed_less_than_once_per_replicate))
        mask = (~expressed_less_than_once_per_replicate)
        goi = mask[mask].index
        
        df = df[goi]

        t = Ktest(data=df,
                metadata=meta.copy(),
                condition='label',
                nystrom=nystrom,
                center_by=center_by,
                data_name=data_name)
        
        t.add_zero_proportions_to_var()
        kernel_ = 'fisher_zero_inflated_gaussian'
        kernel_info_ = [f'{k}_pz' for k in t.get_samples_list()]


    else:
        kernel_ = None
        kernel_info_ = None

    t.univariate_test(n_jobs=n_jobs,
            lots=lots,
            fromto=[v0,v1],
            save_path=diroutput,name='DEA',
            ignore_zeros=ignore_zeros,
            kernel=kernel_,
            kernel_info = kernel_info_,
            verbose=4)

