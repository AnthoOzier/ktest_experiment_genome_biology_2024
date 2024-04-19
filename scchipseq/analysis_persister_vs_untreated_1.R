source("get_sub_populations.R")
treat1           = "persister"
treat2           = "untreated_1"
sc_df_treat      = rbind(sc_df_persister,  sc_df_untreated_1)
meta_sc_df_treat = rbind(meta_sc_df_persister,  meta_sc_df_untreated_1)


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


kt = ktest(sc_df_treat, meta_sc_df_treat,condition='new_pop', samples=c(treat1,treat2), verbose=1, nystrom = TRUE)
kt$univariate_test(
              n_jobs    = 24L, 
              save_path = "/results/",
              name      = paste0(treat1,"_",treat2),
              kernel    = list('function'='gauss','bandwidth'='median'),
              verbose   = 1L
)

x = CreateSeuratObject(counts = t(exp(sc_df_treat)))
x@meta.data$cluster      = meta_sc_df_treat$new_pop
new_pop                  = meta_sc_df_treat$new_pop
Idents(x)                = new_pop
x@meta.data$new_pop      = new_pop
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

ff  = list.files("/results/",full.names=T)
fff = names(unlist(sapply(ff,FUN=function(x){grep("ny_lmrandom_m",x)*grep(treat2,x,fixed=TRUE)})))
fff = names(unlist(sapply(fff,FUN=function(x){grep("ny_lmrandom_m",x)*grep(treat1,x,fixed=TRUE)})))
DE_genes_ktest            = read.table(fff,h=T,sep=",")
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

df_mean           = (t(apply(sc_df_treat,2,FUN=function(x){tapply(x,meta_sc_df_treat$new_pop,mean)})))
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

tt = read.table(file = paste0("/results/DEG_",treat1,"_",treat2,".txt"),h=T)
tt$locus_id = paste(tt[,1],tt[,2],tt[,3],sep="_")
ww = lapply(1:dim(tt)[1],FUN=function(i){tt[tt$locus_id == tt$locus_id[i],]$gene_id})
ww = lapply(ww,FUN=function(y){unique(unlist(lapply(strsplit(y,c("_")),FUN=function(x){x[[1]][1]})))})
ww = lapply(ww, FUN=function(x){paste(x,collapse="-")})
tt$gene_id = Reduce("c",ww)
tt = tt[!duplicated(tt$locus_id),]

write.table(tt, file = paste0("/results/DEG_",treat1,"_",treat2,".txt"),quote=F,col.names=T,row.names=F,sep="\t")


