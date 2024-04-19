
rm(list=ls())
library(tidyverse)
library(parallel)
library(sctransform)
library(magrittr)
library(Matrix)
library(peakRAM)
library(future)
suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(SingleR)
  library(SingleCellExperiment)
  library(BiocManager)
})
source("/home/scripts/DE-analysis-master/R/functions/get_comparisons.R") 
source("/home/scripts/DE-analysis-master/R/functions/run_DE.R") 




dirdata ="/home/data/squair/rds/"
path_output = "/home/data/squair/normalized/"


# all genes sct 
for (file in list.files(dirdata)){
    print(file)
    input_file=paste(dirdata,file,sep="")
    dataset = input_file %>% basename() %>% gsub("\\.rds$", "", .)
    output_filename = paste0("normalized_",dataset,".csv")
    output_metafilename = paste0("meta_",dataset,".csv")
    output_file = paste0(path_output,output_filename)
    output_metafile = paste0(path_output,output_metafilename)
    print(output_file)

    if (file.exists(output_file)){
        print("already normalized")
    }else{
        x = readRDS(paste(dirdata,file,sep=''));  
        x = NormalizeData(x)
        normx = x@assays$RNA@data
        meta = x@meta.data
        write.csv(normx,output_file) 
        write.csv(meta,output_metafile) 

        
    }
}



# DE analysis with existing methods

args = commandArgs(trailingOnly=TRUE)
file = args[1]
diroutput= args[2]
dirdata = args[3]

print(paste("file =",file,"diroutput =",diroutput,"dirdata =",dirdata))


de_tests = c(
    ## single-cell methods, implemented in Seurat
    "wilcox", "bimod", "t", "negbinom", "poisson", "LR", 
    "MAST",
    ## pseudobulk methods
    "pseudobulk_DESeq2,test?LRT", "pseudobulk_DESeq2,test?Wald",
    "pseudobulk_limma,mode?voom", "pseudobulk_limma,mode?trend",
    "pseudobulk_edgeR,test?QLF", "pseudobulk_edgeR,test?LRT"

    ## pseudobulk methods run without aggregation
    # "pseudobulk_DESeq2,test?LRT,replicate?cells",
    # "pseudobulk_DESeq2,test?Wald,replicate?cells",
    # "pseudobulk_limma,mode?voom,replicate?cells",
    # "pseudobulk_limma,mode?trend,replicate?cells",
    # "pseudobulk_edgeR,test?QLF,replicate?cells",
    # "pseudobulk_edgeR,test?LRT,replicate?cells"
  )

input_file=paste(dirdata,file,sep="")
for (de_test in de_tests){
   
    print(paste(input_file,de_test))
    
    dataset = input_file %>% basename() %>% gsub("\\.rds$", "", .)
    output_filename = paste0(dataset, "-de_test=", de_test,".csv")    
    output_file = file.path(diroutput,output_filename)
    
    if (file.exists(output_file)){
      print("already tested")
    }else{
      print('testing')
      
      sc = readRDS(input_file) 
      expr = GetAssayData(sc, slot = 'counts')
      meta = sc@meta.data
      
      DE = run_DE(sc, de_test = de_test)
      write_csv(DE,output_file)
    }
}


