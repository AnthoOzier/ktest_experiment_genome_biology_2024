
rename_de_table = function(table){
    # fonctionne pour le bulk et pour le sc
    colnames(table) %<>%
        fct_recode('p_val' = 'p.value',  ## DESeq2
                   'p_val' = 'pvalue',  ## DESeq2
                   'p_val' = 'p.value',  ## t/wilcox
                   'p_val' = 'P.Value',  ## limma
                   'p_val' = 'PValue'  , ## edgeR
                   'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                   'p_val_adj' = 'adj.P.Val',      ## limma
                   'p_val_adj' = 'FDR',            ## edgeER
                   'avg_logFC' = 'log2FoldChange', ## DESEeq2
                   'avg_logFC' = 'logFC', ## limma/edgeR
                   'test_statistic' = 'stat', ## DESeq2
                   'test_statistic' = 'F', ## edgeR
                   'test_statistic' = 't', ## limma
                   'test_statistic' = 'LR', ## edgeR LRT
                   'test_statistic' = 'statistic' ## t
        ) %>%
        as.character()
    return(table)
}

rename_kfda_de_table = function(table,name,trunc,verbose=0){
    if (verbose >0){
        print(paste('- rename kfda de table'))
    }
    colnames(table) %<>%
        fct_recode('p_val' = paste(name,'_t',trunc,'_pval',sep=''),  ## DESeq2
                   'p_val_adj' = paste(name,'_t',trunc,'_pvalBH',sep=''),  ## DESeq2
                   'test_statistic' = paste(name,'_t',trunc,'_kfda',sep='') ## DESeq2
        ) %>%
        as.character()
    return(table)
}

pvalue_preparation = function(table){
    # remove NAs
    table %<>% filter(!is.na(p_val), !is.na(p_val_adj)) # , !is.na(test_statistic))
    # replace p=0 with minimum p-value
    min_pval_adj = min(table$p_val_adj[table$p_val_adj > 0])
    table %<>% 
        mutate(p_val_adj = ifelse(p_val_adj <= min_pval_adj, min_pval_adj, p_val_adj))
    ## repeat for raw p-values
    min_pval = min(table$p_val[table$p_val > 0])
    table %<>% 
        mutate(p_val = ifelse(p_val <= min_pval, min_pval, p_val))
    return(table)
}

get_prepared_pval_table = function(table,kfda=FALSE,name=NA,trunc=NA){
    if(kfda){
        table = rename_kfda_de_table(table,name,trunc)
    }else{
        table = rename_de_table(table)
    }
    table = pvalue_preparation(table)
    return(table[c('gene','p_val','p_val_adj')])
}

get_kfda_de = function(kfdadir,kfdafile,verbose=0){
    if (verbose>0){
        print('- Read kfda DE file')
        print(paste('file :',kfdafile))
    }
    df = read_csv(paste0(kfdadir,kfdafile),show_col_types = FALSE)
    return(df)
}

get_kfda_de_table_cbis_replicate = function(de_method,dataset,fisher_kernel=FALSE,verbose=0){
    if (grepl('kfzig',de_method)){
        kstr = 'kfzig_'
    }else if(grepl('kgq20',de_method)){
        kstr = 'kgq20_'
    }else if(grepl('kgq80',de_method)){
        kstr = 'kgq80_'
    }else if(grepl('klin',de_method)){
        kstr = 'klin_'
    }else{
        kstr = ''
    }
    if ((grepl('Angelidis',dataset) || grepl('Reyfman',dataset))){
        cbstr = 'crlis_cnz'
        dir=paste0(cbis_dir,'Courtine_KFDA_','cbrlis_cbnz_',kstr,'results_adj/')
    }else{
        cbstr = 'cris_cnz'
        dir=paste0(cbis_dir,'Courtine_KFDA_','cbris_cbnz_',kstr,'results_adj/')
    }
    for (kfdafile in list.files(dir)){
        if (grepl(dataset,kfdafile) && grepl(cbstr,kfdafile)){
            print(paste('kfdafile',kfdafile))
            table = get_kfda_de(dir,kfdafile)
            trunc = strtoi(sub(',','',strsplit(de_method,'_t')[[1]][2]))
            name = sub("_univariate.csv$", "",kfdafile)
            table = rename_kfda_de_table(table,name,trunc)
            table = pvalue_preparation(table)
        }
    }
    return(table)
}





get_kfda_de_table = function(de_method,dataset,verbose=0){
    if (verbose>0){
        print('- Get kfda de table ')
        print(paste('de_method:',de_method,'datset:',dataset))
    }
    if (grepl('kfda_cbisreplicate_',de_method)){
        print('center by replicate in input space')
        table = get_kfda_de_table_cbis_replicate(de_method=de_method,dataset=dataset,verbose=verbose)
    }else{
        print(paste('method',de_method,'not recognized'))
    }
    return(table[c('gene','p_val','p_val_adj')])
}

get_sc_de_table = function(de_method,dataset,verbose=0){
    file = paste0(dataset,'-de_test=',de_method,'.csv')
    file.exists.in.dir = file %in% list.files(scdir)
    print(paste(file,file.exists.in.dir))
    if (file.exists.in.dir){
        table = read_csv(paste0(scdir,file),show_col_types = FALSE)
        if (verbose>1){
            print(colnames(table))
        }
        table = rename_de_table(table)
        table = pvalue_preparation(table)
        table = table[c('gene','p_val','p_val_adj')]
    }
    return(table)
}

get_bulk_de_table = function(de_method,dataset,verbose=0){
    if (grepl('Reyfman',dataset) & de_method=='bulk_DESeq2,test?LRT'){
        print('Reyfman')
        file = paste0(dataset,'_results.tsv.gz')
        table = read_tsv(paste0(bulkdir,file),show_col_types=FALSE)
        table = rename_de_table(table)
        table = pvalue_preparation(table)
        table = table[c('gene','p_val','p_val_adj')]
        
    }else{
        file = paste0(dataset,'-de_test=',de_method,'.rds')
        file.exists.in.dir = file %in% list.files(bulkdir)
        print(paste(file,'exists ?',file.exists.in.dir))
        if (file.exists.in.dir){
            table = read_csv(paste0(bulkdir,file),show_col_types = FALSE)
            table = rename_de_table(table)
            table = pvalue_preparation(table)
            table = table[c('gene','p_val','p_val_adj')]
        }
    }
    return(table)
}

get_de_table = function(de_method,dataset,verbose=0){
    if (verbose >0){
        print('- Get de table')
    }
    if (grepl('kfda',de_method,fixed=TRUE)){
        table = get_kfda_de_table(de_method,dataset,verbose=verbose)
    }else if (substr(de_method,1,4)=='bulk'){
        if (verbose>0){print('- bulk data')}
        table = get_bulk_de_table(de_method,dataset,verbose=verbose)
    }else{
        table = get_sc_de_table(de_method,dataset,verbose=verbose)
    }
    return(table)
}

reduce_to_intersection = function(table1,table2,verbose=0){
    if (verbose>0){
        print('- Reduce to intersection')
    }
    # filter to genes detected in both tables
    genes = intersect(table1$gene, table2$gene)
    table1 %<>% filter(gene %in% genes) %>% arrange(gene)
    table2 %<>% filter(gene %in% genes) %>% arrange(gene)
    return(list(table1=table1,table2=table2))
}

compute_aucc = function(table1,table2,k,verbose=0){
    if (verbose>0){
        print('- compute aucc')
        print(paste('k =',k))
    }
    # area under the concordance curve
    k = as.integer(k)
    ## rank in descending order first by p_val
    ## break ties by the abs() of the test_statistic
    vec1 = table1 %>%
        arrange(p_val, desc(abs(p_val))) %>%
        # arrange(p_val, desc(abs(test_statistic))) %>%
        pull(gene) %>%
        head(k)
    vec2 = table2 %>%
        # arrange(p_val, desc(abs(test_statistic))) %>%
        arrange(p_val, desc(abs(p_val))) %>%
        pull(gene) %>%
        head(k)
    
    concordance_curve = map_dbl(seq_len(k), ~ {
        v1 = vec1[seq_len(.)]
        v2 = vec2[seq_len(.)]
        length(intersect(v1, v2))
    })
    denom = k * (k + 1) / 2
    aucc = sum(concordance_curve) / denom
    # print(aucc)
    return(aucc)
}

aucc_of_dataset = function(dataset,de_method1,de_method2,k,verbose=0){
    if (verbose>0){
        print('aucc of dataset')
    }
    print('- Load table 1')
    table1 = get_de_table(de_method1,dataset,verbose=verbose)
    print('- Load table 2')
    table2 = get_de_table(de_method2,dataset,verbose=verbose)
    tables = reduce_to_intersection(table1,table2,verbose=verbose)
    table1 = tables[[1]]
    table2 = tables[[2]]
    aucc = compute_aucc(table1,table2,k,verbose=verbose)
    auccs[dataset] = aucc
    print(paste('aucc =',aucc))
    return(aucc)
}

