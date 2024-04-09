library(tidyverse)
library(magrittr)
library(future)
library(progressr)
library(foreach)

# load scripts from the original study here : https://github.com/neurorestore/DE-analysis
source("/home/scripts/DE-analysis-master/R/functions/calculate_overlap.R") 
source("/home/scripts/utils_aucc.R") 
verbose=3

cbis_dir = '/home/data/squair/ktest_results/'
bulkdir = '/home/data/squair/bulk_results/'
scdir = '/home/data/squair/sc_results/'
path_aucc= '/home/data/squair/AUCC/couples/'

datasets = c(
    'Reyfman2020_pneumo','Reyfman2020_alvmac','CanoGamez2020_Naive-Th2', 'CanoGamez2020_Naive-iTreg',
    'CanoGamez2020_Naive-Th0', 'CanoGamez2020_Naive-Th17','CanoGamez2020_Memory-iTreg','CanoGamez2020_Memory-Th0',
    'CanoGamez2020_Memory-Th17','CanoGamez2020_Memory-Th2','Angelidis2019_alvmac','Angelidis2019_pneumo',
    'Hagai2018_mouse-pic', 'Hagai2018_rat-lps','Hagai2018_rabbit-lps', 'Hagai2018_rat-pic',
    'Hagai2018_mouse-lps', 'Hagai2018_pig-lps'
    )

datasets_w_Reyfman = datasets
datasets_wo_Reyfman = datasets[!grepl('Reyfman',datasets)]

de_methods1 = c(
    'kfda_cbisreplicate_sn_t4',
    'kfda_cbisreplicate_klin_sn_t1',
    'kfda_cbisreplicate_kfzig_sn_t4'
    )
de_methods2 = c("bimod","LR","negbinom","poisson","pseudobulk_DESeq2,test?LRT,replicate?cells",
    "pseudobulk_DESeq2,test?LRT","pseudobulk_DESeq2,test?Wald,replicate?cells",
    "pseudobulk_DESeq2,test?Wald","pseudobulk_edgeR,test?LRT,replicate?cells",
    "pseudobulk_edgeR,test?LRT","pseudobulk_edgeR,test?QLF,replicate?cells",
    "pseudobulk_edgeR,test?QLF","pseudobulk_limma,mode?trend,replicate?cells",
    "pseudobulk_limma,mode?trend","pseudobulk_limma,mode?voom,replicate?cells",
    "pseudobulk_limma,mode?voom","t","wilcox",'MAST'
)

bulk_methods = c(
    "bulk_DESeq2,test?LRT",  "bulk_DESeq2,test?Wald", "bulk_edgeR,test?LRT",  
    "bulk_edgeR,test?QLF",   "bulk_limma,mode?trend", "bulk_limma,mode?voom"
)

de_methods = c(de_methods1,de_methods2,bulk_methods)

get_adapted_dataset_list = function(de_method1,de_method2){
    datasets = datasets_w_Reyfman
    if (substr(de_method1,1,4)=='bulk'){
        if (de_method1 != "bulk_DESeq2,test?LRT"){
            datasets = datasets_wo_Reyfman
        }
    }
    if (substr(de_method2,1,4)=='bulk'){
        if (de_method2 != "bulk_DESeq2,test?LRT"){
            datasets = datasets_wo_Reyfman
        }
    }
    return(datasets)
}


aucc_of_two_methods = function(de_method1,de_method2,k){
    print(paste('de_method1 =',de_method1))
    print(paste('de_method2 =',de_method2))
    df_aucc = data_frame()
    df_aucc[1,]=NA
    auccs=list()
    aucc_name = paste('de1',de_method1,'de2',de_method2,paste0('k',k),sep='_')
    datasets = get_adapted_dataset_list(de_method1,de_method2)
    if (!file.exists(paste(path_aucc,aucc_name,'.csv',sep=''))){
        auccs = lapply(datasets,aucc_of_dataset,de_method1=de_method1,de_method2=de_method2,k=k,verbose=verbose)
        names(auccs)=datasets
        print(length(auccs))
        df_aucc = as.data.frame(auccs,row.names=aucc_name)
        print(df_aucc)
        fileaucc = paste(path_aucc,aucc_name,'.csv',sep='')
        write.csv(df_aucc,fileaucc)
    }
}


for(k in c(100,200,500)){
    for(de_method1 in de_methods){
        for (de_method2 in de_methods){
            if (de_method1>de_method2){
                print(paste('de_method1 =',de_method1))
                print(paste('de_method2 =',de_method2))
                df_aucc = data_frame()
                df_aucc[1,]=NA
                auccs=list()
                aucc_name = paste('de1',de_method1,'de2',de_method2,paste0('k',k),sep='_')

                datasets = get_adapted_dataset_list(de_method1,de_method2)

                if (!file.exists(paste(path_aucc,aucc_name,'.csv',sep=''))){
                    auccs = lapply(datasets,aucc_of_dataset,de_method1=de_method1,de_method2=de_method2,k=k,verbose=verbose)
                    names(auccs)=datasets
                    print(length(auccs))
                    df_aucc = as.data.frame(auccs,row.names=aucc_name)
                    print(df_aucc)
                    fileaucc = paste(path_aucc,aucc_name,'.csv',sep='')
                    write.csv(df_aucc,fileaucc)
                }
            }
        }
    }
}

