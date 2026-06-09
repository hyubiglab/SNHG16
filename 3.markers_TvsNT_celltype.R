library(Seurat)
library(plyr)

runDeg <- function(input_object, outputPath, outputSuffix) {
  library(plyr)
  library(dplyr)
  library(ggplot2)
  
  degmarkers <- FindAllMarkers(object = input_object, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
  degmarkers$Symbol <- unlist(lapply(degmarkers$gene, function (x) unlist(strsplit(as.character(x), split = '-ENSG'))[1] ))
  write.table(degmarkers, paste0(c(outputPath, '/', 'markers.', outputSuffix, '.txt'), collapse = ''), sep = '\t', quote = F, col.names = T, row.names = F) 
}


myObj <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')
head(myObj@meta.data)

summary(myObj@meta.data)
levels(myObj@meta.data$Cell_type_v2)
dir.create('degs_tvsnt')

for (celltype in levels(myObj@meta.data$Cell_type_v2)) {
  if (celltype != 'Malignant cells') {
    tmp <- subset(myObj, subset = Cell_type_v2 == celltype)
    tmp@meta.data <- droplevels(tmp@meta.data)
    Idents(tmp) <- 'Condition'
    
    if (length(levels(tmp@meta.data$Condition)) != 1) {
      if (min(summary(tmp@meta.data$Condition)) >= 3) {
        print (celltype)
        celltype <- gsub(celltype, pattern = ' ', replacement = '')
        celltype <- gsub(celltype, pattern = '/', replacement = '_')
        
        runDeg(input_object = tmp, outputPath = 'degs_tvsnt', 
               outputSuffix = paste0(c('Cell_celltype.', celltype, '.T_', summary(tmp@meta.data$Condition)[1], '_vs_NT_', summary(tmp@meta.data$Condition)[2]), collapse = ''))
      }
    }
    remove(tmp)
  }
}



epinormal_bc <- rownames(subset(myObj@meta.data, Condition == 'Normal' & Cell_type_v2 == 'Epithelial cells'))
maligannt_bc <- rownames(subset(myObj@meta.data, Cell_type_v2 == 'Malignant cells'))
tmp <- subset(myObj, cells = c(epinormal_bc, maligannt_bc))
tmp@meta.data <- droplevels(tmp@meta.data)
summary(tmp@meta.data)
Idents(tmp) <- 'Condition'
runDeg(input_object = tmp, outputPath = 'degs_tvsnt', 
       outputSuffix = paste0(c('Cell_celltype.Malignantcells.T_', length(maligannt_bc), '_vs_NT_', length(epinormal_bc)), collapse = ''))


