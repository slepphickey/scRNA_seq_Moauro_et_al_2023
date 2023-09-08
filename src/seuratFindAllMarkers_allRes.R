#'args[1] path to seurat object
#'args[2] cluster assignment column prefix ie "SCT_snn_res"
#'args[3] adjusted pval filter
#'args[4] path to output dir
args <- commandArgs(TRUE)
library(Seurat)
library(dplyr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

print(paste0("path to seurat object: ", args[1]))
print(paste0("prefix: ", args[2]))
print(paste0("adjp = ", args[3]))
print(paste0("outdir: ", args[4]))



seuratOb = loadRData(args[1])
prefix = args[2]
adjp = as.numeric(args[3])
out_dir = args[4]

# if there is a res = 0 remove it
res_oi = colnames(seuratOb@meta.data)[grep(prefix, colnames(seuratOb@meta.data))]
res0 = paste0(prefix, ".0")
res_oi = res_oi[!res_oi %in% res0]

for(res in res_oi){
  
  # set the identity
  Idents(object = seuratOb) <- res
  
  # find all markers
  markers = FindAllMarkers(seuratOb, 
                           only.pos = TRUE,
                           logfc.threshold = 0,
                           min.pct = 0)
   
  # filter for p_val_adj 
  # markers = 
  # markers %>%
  # filter(p_val_adj < adjp) 
  
  # write to a csv
  write.csv(markers,
            file = paste0(out_dir,
                          "/",  
                          res, 
                          "_marker_genes.csv"))
}