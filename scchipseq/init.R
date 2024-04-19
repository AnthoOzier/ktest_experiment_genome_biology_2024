# https://github.com/LMJL-Alea/ktest
# load ktest package
library(ktest)

# load python environment
library(reticulate)
venv <- "ktest"
use_virtualenv(venv)
py_config()

# check ktest
check_ktest()


my_theme <- function(){
  theme(   
    axis.line        = element_line(colour = "black"),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    panel.background = element_rect(fill=NA),
    strip.text.x     = element_text(size = 10),
    strip.text.y     = element_text(size = 10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position  = "top",
    legend.title     = element_blank(),
    legend.text      = element_text(size=6),  
    aspect.ratio     = 1,
    axis.text.x      = element_text(size=5),
    axis.text.y      = element_text(size=5)
  )
}

suppressPackageStartupMessages({
  library(plyr)
  library(knitr)
  library(kableExtra)
  library(ggplot2)
  library(parallel)
  library(DESeq2)
  library(gridExtra)
  library(tidyverse)
  library(magrittr)
  library(Seurat)
  library(Matrix)
  library(peakRAM)
  library(future)
  library(xgboost)
  library(pROC)
  library(caTools)
  library(randomForest)
  library(e1071)
  library(ggridges)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(ggvenn)
  library(NMF)
	library(Matrix)
	library(harmony)
	library(qs)
	library(scran)	
	library(ggpubr)
	library(RColorBrewer)
	library(archive)
	library(readr)
	library(FactoMineR)
  library(Seurat)
  library(dplyr)
  library(sctransform)
  library(SingleR)
  library(patchwork)
  library(cowplot)
  library(stringr)
  library(cluster)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(viridis)
  library(gridExtra)
  library(RColorBrewer)
  library(xlsx)
  library(Nebulosa)
  library(slingshot, quietly = FALSE)
  library(SingleCellExperiment)
  library(destiny, quietly = TRUE)
  library(mclust, quietly = TRUE)
  library(cluster)
  library(enrichR)
  library(harmony)
  library(metap)
  library(BiocManager)
})





GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}







