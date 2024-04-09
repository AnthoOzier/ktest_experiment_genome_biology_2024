from ktest.tester import Ktest
import pandas as pd
import os
from ktest.pvalues import correct_BenjaminiHochberg_pval_of_dfcolumn


path_squair = '/home/data/squair/ktest_results/'

folders = ['cbris_cbnz','cbrlis_cbnz',
          'cbris_cbnz_kfzig','cbrlis_cbnz_kfzig',
          'cbris_cbnz_klin','cbrlis_cbnz_klin',
          'cbris_cbnz_kgq20','cbrlis_cbnz_kgq20',
          'cbris_cbnz_kgq80','cbrlis_cbnz_kgq80',
          ]


datasets = {}
# load results  
for folder in folders: 
    print('\n',folder)
    path = f'{path_squair}Courtine_KFDA_{folder}_results/'

    datasets[folder] = {}
    for file in os.listdir(path):
        if file[0] == '0':
            dataset = file.split(sep='sn__')[1].split(sep='label_cb')[0]
            if dataset not in datasets[folder]:
                datasets[folder][dataset]={}
            else:
                print(dataset,file)
            spec = file.split(sep='DEA')[1].split(sep='.csv')[0]
            datasets[folder][dataset][spec] = file
    print(folder,'\n',len(list(datasets[folder].keys())),list(datasets[folder].keys()))

# aggregate results 
for folder in folders:
    path = f'{path_squair}Courtine_KFDA_{folder}_results/'
    path_output = f'{path_squair}Courtine_KFDA_{folder}_results_assembled/'
    print(folder)
    for dataset,specs in datasets[folder].items() : 
        print('\n',dataset,len(specs))    
        for spec,file0 in specs.items():
            print('\n',spec)
            df0 = pd.read_csv(path+file0,index_col=0)
            if df0.shape[0] == 1000:
                file = file0.replace('0_1000','')
            else:
                file = file0.replace(file0.split(sep='DEA')[0],'')
            if file in os.listdir(path_output):
                print(dataset,spec,'already assembled')
            else:
                if df0.shape[0] == 1000:
                    for s in range(1000,26000,1000):
                        file = file0.replace('0_1000',f'{s}_{s+1000}')
                        if file in os.listdir(path):
                            df = pd.read_csv(path+file,index_col=0)
                            print('loaded',df.shape[0],s)
                            df0 = pd.concat((df,df0))
                        else:
                            print(file)
                    file = file0.replace('0_1000','')
                    print(file)
                else:
                    file = file0.replace(file0.split(sep='DEA')[0],'')
                    print(file,df0.shape)
                df0 = df0[~df0.index.duplicated(keep='first')]
                df0.to_csv(path_output+file)



# p-value corrections 
for folder in folders:
    print(folder)
    path_output = f'{path_squair}Courtine_KFDA_{folder}_results_assembled/'
    path_adj = f'{path_squair}Courtine_KFDA_{folder}_results_adj/'
    for file in os.listdir(path_output):
        if file not in os.listdir(path_adj):
            print(file)
            df = pd.read_csv(path_output+file,index_col=0)
            for c in df.columns:
                for trunc in range(1,11):
                    if f't{trunc}_pval' in c and 'BH' not in c:
                        df[f'{c}BH'] = correct_BenjaminiHochberg_pval_of_dfcolumn(df[c])
            df['gene'] = df.index
            df.to_csv(path_adj+file)
        else:
            print('###',file,'already corrected')      


