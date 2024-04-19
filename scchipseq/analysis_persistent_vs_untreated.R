rm(list=ls())
source("init.R")

mincells = 5
tmax     = 20 
treat1   = "persister"
treat2   = "untreated"

sce        = qread("/data/MM468_H3K27me3_peaks_3000_10000_95_uncorrected_correlation_filtered.qs")
sce        = logNormCounts(sce)
tmp        = modelGeneVar(sce)
chosen     = getTopHVGs(tmp, prop = 1)
sce        = sce[chosen,]
sc_df      = as.data.frame(t(as.matrix(assay(sce,"counts"))))
meta_sc_df = as.data.frame(colData(sce))
sc_df      = log(sc_df+1)
ncells     = dim(sc_df)[1]

day   = sapply(rownames(meta_sc_df), FUN=function(s){strsplit(s,"_")[[1]][3]})

persister_cells = c(
  names(unlist(sapply(day,FUN=function(s){grep("day33",s)}))),
  names(unlist(sapply(day,FUN=function(s){grep("day67",s)})))
  )

untreated_cells = c(
  names(unlist(sapply(day,FUN=function(s){grep("day60",s)}))),
  names(unlist(sapply(day,FUN=function(s){grep("day77",s)}))),
  names(unlist(sapply(day,FUN=function(s){grep("day131",s)})))
)

untreated_cells = names(unlist(sapply(untreated_cells,FUN=function(s){grep("DMSO",s)})))

meta_sc_df[ rownames(meta_sc_df) %in% persister_cells, ]$cell_cluster = "persister"
meta_sc_df[ rownames(meta_sc_df) %in% untreated_cells, ]$cell_cluster = "untreated"


sc_df_treat      = sc_df[ ( meta_sc_df$cell_cluster =="untreated" | meta_sc_df$cell_cluster =="persister"),]
meta_sc_df_treat = meta_sc_df[ ( meta_sc_df$cell_cluster =="untreated" | meta_sc_df$cell_cluster =="persister"),]

keep_genes       = apply(sc_df_treat,2,FUN=function(z){
  keep = TRUE; 
  if (sum(z!=0)==0){
    keep = FALSE
  } else if (sum(z==0)>0){
    tt   = table(meta_sc_df_treat$cell_cluster,z!=0)
    if (sum(tt[,2])<mincells){
      keep = FALSE
    }  
  }
  return(keep)
})
sc_df_treat  = sc_df_treat[,keep_genes]
gene_list    = sort(colnames(head(sc_df_treat)))
gene_list_df = Reduce("rbind",mapply(gene_list,FUN=function(s){strsplit(s,"_")}))


kt = ktest(sc_df_treat, meta_sc_df_treat,condition='cell_cluster', samples=c(treat1,treat2), verbose=1)
kt$multivariate_test(verbose=1)
kt$projections()    
res = list(
  pv        = kt$get_pvalue()[1:tmax],        
  proj_kfda = kt$init_df_proj(proj='proj_kfda')[,1:tmax],
  proj_kpca = kt$init_df_proj(proj='proj_kpca')[,1:tmax]
  )
save(res,file=paste0("/results/multivariate_results_",treat1,"_",treat2,".RData"))


kt = ktest(sc_df_treat, meta_sc_df_treat,condition='cell_cluster', samples=c(treat1,treat2), verbose=1, nystrom = TRUE)
kt$univariate_test(
              n_jobs    = 24L, 
              save_path = "/results/",
              name      = paste0(treat1,"_",treat2),
              kernel    = list('function'='gauss','bandwidth'='median'),
              verbose   = 1L
)


x = CreateSeuratObject(counts = t(exp(sc_df_treat)))
x@meta.data$cluster      = meta_sc_df_treat$cell_cluster
cell_cluster             = meta_sc_df_treat$cell_cluster
Idents(x)                = cell_cluster
x@meta.data$cell_cluster = cell_cluster
x = SCTransform(x)
x = FindVariableFeatures(x, nfeatures = 660)

DE_genes = FindMarkers(x, 
            ident.1 = treat1, ident.2 = treat2, 
            min.cells.group = 0, 
            min.cells.feature = 0,
            min.pct = 0,
            logfc.threshold = -Inf,
            only.pos = FALSE,
            test.use = "wilcox"
)
DE_genes           = DE_genes[,c(1,2)]
colnames(DE_genes) = c("p_val","avg_log2FC")
DE_genes_wilcox    = DE_genes
gene_list_df       = Reduce("rbind",mapply(rownames(DE_genes_wilcox),FUN=function(s){str_replace_all(s, "-", "_")}))
DE_genes_wilcox    = data.frame(ave_log2FC = DE_genes_wilcox$avg_log2FC)
rownames(DE_genes_wilcox) = gene_list_df 

DE_genes_ktest            = read.table("/results/persister_untreated_ny_lmrandom_m499_basisw_datacell_clusterpersisteruntreated_univariate.csv",h=T,sep=",")
t                         = 5
gene_list                 = DE_genes_ktest[,1]
gene_list_coord           = Reduce("rbind",mapply(gene_list,FUN=function(s){strsplit(s,"_")}))
DE_genes_ktest            = DE_genes_ktest[,c(2*t,2*t+1)]
colnames(DE_genes_ktest)  = c("pv","kfda")
rownames(DE_genes_ktest)  = gene_list
DE_genes_ktest$chr        = gene_list_coord[,1]
DE_genes_ktest$start      = gene_list_coord[,2]
DE_genes_ktest$end        = gene_list_coord[,3]
DE_genes_ktest            = DE_genes_ktest[,c(3,4,5,1,2)]
DE_genes_ktest$adj_pval   = p.adjust(DE_genes_ktest$pv,"bonferroni")
DE_genes_ktest$rank       = rank(DE_genes_ktest$kfda)

df_mean           = (t(apply(sc_df_treat,2,FUN=function(x){tapply(x,meta_sc_df_treat$cell_cluster,mean)})))
colnames(df_mean) = paste("average",sort(c(treat1,treat2)),sep="_")
DE_genes_ktest    = merge(DE_genes_ktest,df_mean,by="row.names")
rownames(DE_genes_ktest) = DE_genes_ktest[,1]
DE_genes_ktest    = DE_genes_ktest[,-1]
DE_genes_ktest    = DE_genes_ktest[order(DE_genes_ktest$rank,decreasing=TRUE),]
DE_genes_ktest    = merge(DE_genes_ktest,DE_genes_wilcox, by ="row.names")
rownames(DE_genes_ktest) = DE_genes_ktest[,1]
DE_genes_ktest    = DE_genes_ktest[,-1]
o                 = order(DE_genes_ktest$rank,decreasing=TRUE)
DE_genes_ktest    = DE_genes_ktest[o,]


#### need bedtools to process the genomic data 
#### https://bedtools.readthedocs.io/en/latest/

file.out         = paste0("/results/",treat1,"_vs_",treat2,"_ktest_coordinates.bed")
file.out.sorted  = paste0("/results/",treat1,"_vs_",treat2,"_ktest_coordinates_sorted.bed")
file.out.merged  = paste0("/results/",treat1,"_vs_",treat2,"_ktest_coordinates_genes.bed")
write.table(DE_genes_ktest[,1:3], file = file.out,quote=F,col.names=F,row.names=F,sep="\t")
system(paste("bedtools sort -i ",file.out," > ",file.out.sorted,sep=""))
system(paste("mv ",file.out.sorted," ",file.out,sep=""))
gene.file = "/data/human_epdnew_qoQwV.bed"
system(paste("bedtools intersect -wa -wb -a ",file.out," -b ",gene.file," > ",file.out.merged,sep=""))
DE_genes_ktest_annotation = read.table(file.out.merged)
colnames(DE_genes_ktest_annotation) = c("chr","start","end","tss_chr","tss_start","tss_end","gene_id","gene_exon","tss_strand")
tt = merge(DE_genes_ktest,DE_genes_ktest_annotation, by=c("chr","start","end"))
o  = order(tt$rank,decreasing=TRUE)
tt = tt[o,]

write.table(tt, file = paste0("/results/DEG_",treat1,"_",treat2,".txt"),quote=F,col.names=T,row.names=F,sep="\t")

