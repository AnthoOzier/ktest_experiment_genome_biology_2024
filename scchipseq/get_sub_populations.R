library(Mclust)
mincells      = 5
tmax          = 20 
t             = 10
nb_genes      = 100

treat1 = "persister"
treat2 = "untreated"

load(file=paste0("/results/multivariate_results_",treat1,"_",treat2,".RData"))
t        = 5
dd       = meta_sc_df_treat
dd$index = rownames(dd)
res$proj_kfda$index = rownames(res$proj_kfda)
res$proj_kpca$index = rownames(res$proj_kpca)
cn       = c("index",colnames(dd)[1:13])
dd       = merge(dd,res$proj_kfda[,c("index",t)],by="index");colnames(dd) = c(cn,"proj_kfda")
dd       = merge(dd,res$proj_kpca[,c("index",t+1)],by="index"); colnames(dd)[16] = "proj_kpca"


x    = dd[dd$cell_cluster==treat2,]$proj_kfda
id   = dd[dd$cell_cluster==treat2,]$index
K    = 3
out  = Mclust(x,K,"E")
dens = lapply(1:K,FUN=function(k){out$parameters$pro[k]*dnorm(dd$proj_kfda, mean = out$parameters$mean[k], sd = sqrt(out$parameters$variance$sigmasq))})
lik  = apply(Reduce("rbind",dens),2,sum) 

dd$density_stage_1 = dens[[1]]
dd$density_stage_2 = dens[[2]] 
dd$density_stage_3 = dens[[3]] 

gg = ggplot(dd, aes(x = proj_kfda, color=cell_cluster)) +  
        geom_density()+
        geom_line(aes(x = proj_kfda,y = density_stage_1), color = brewer.pal(9, "Blues")[5], linetype="dashed") +
        geom_line(aes(x = proj_kfda,y = density_stage_2), color = brewer.pal(9, "Blues")[7], linetype="dashed") +
        geom_line(aes(x = proj_kfda,y = density_stage_3), color = brewer.pal(9, "Blues")[9], linetype="dashed") +
        my_theme() + #ggtitle("Persister vs. Untreated")+        
        scale_color_manual(values = c(brewer.pal(9, "YlOrRd")[7],"black"),labels=c("Persister","Untreated"))+
        annotate("text", x =(-6000), y = (2.5e-4), label = "Persister Cells", color = brewer.pal(9, "YlOrRd")[7], size=6)  +
        annotate("text", x =(-1107.159), y = (0.25e-4), label = "Persister-like Cells", color = brewer.pal(9, "Blues")[5], size=6)  +
        annotate("text", x =(6830.357) , y = (1.2e-4), label = "Intermediate Cells", color = brewer.pal(9, "Blues")[7], size=6)  +
        annotate("text", x =(12783.024), y = (0.9e-4), label = "Naive Cells", color = brewer.pal(9, "Blues")[9], size=6)  +
        theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank())

x    = dd[dd$cell_cluster==treat2,]$proj_kfda
id   = dd[dd$cell_cluster==treat2,]$index
K    = 3
out  = Mclust(x,K,"E")
z    = out$classification
dens = lapply(1:K,FUN=function(k){out$parameters$pro[k]*dnorm(x, mean = out$parameters$mean[k], sd = sqrt(out$parameters$variance$sigmasq))})
lik  = apply(Reduce("rbind",dens),2,sum) 

dd$density_stage_1[dd$cell_cluster==treat2] = dens[[1]]
dd$density_stage_2[dd$cell_cluster==treat2] = dens[[2]] 
dd$density_stage_3[dd$cell_cluster==treat2] = dens[[3]] 

cells_persister   = dd[dd$cell_cluster==treat1,]$index
cells_untreated_1 = id[z==1]
cells_untreated_2 = id[z==2]
cells_untreated_3 = id[z==3]

write.table(cells_untreated_1, file = "/results/cells_untreated_1.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(cells_untreated_2, file = "/results/cells_untreated_2.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(cells_untreated_3, file = "/results/cells_untreated_3.txt",quote=F,col.names=F,row.names=F,sep="\t")


sc_df_persister   = sc_df[ rownames(sc_df) %in% cells_persister, ]
sc_df_untreated_1 = sc_df[ rownames(sc_df) %in% cells_untreated_1, ]
sc_df_untreated_2 = sc_df[ rownames(sc_df) %in% cells_untreated_2, ]
sc_df_untreated_3 = sc_df[ rownames(sc_df) %in% cells_untreated_3, ]

meta_sc_df_persister           = meta_sc_df[ rownames(meta_sc_df) %in% cells_persister,]
meta_sc_df_persister$new_pop   = "persister"
meta_sc_df_untreated_1         = meta_sc_df[ rownames(meta_sc_df) %in% cells_untreated_1,]
meta_sc_df_untreated_1$new_pop = "untreated_1"
meta_sc_df_untreated_2         = meta_sc_df[ rownames(meta_sc_df) %in% cells_untreated_2,]
meta_sc_df_untreated_2$new_pop = "untreated_2"
meta_sc_df_untreated_3         = meta_sc_df[ rownames(meta_sc_df) %in% cells_untreated_3,]
meta_sc_df_untreated_3$new_pop = "untreated_3"

#########################################################
