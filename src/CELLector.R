### load packages ###
library(devtools)
library(CELLector)

options(header = T)
options(stringsAsFactors = F)

# setwd('/Volumes/home/workspace/CELLector_pancreatic/')

### load input ###
TGCA_mat <- read.csv('data/tcga_snv_matrix_drivers.csv')
inhouse_mat <- read.csv('data/inhouse_snv_matrix_drivers.csv')

convert_input_to_matrix <- function(df){
  
  rownames(df) <- df[,1]
  mat <- as.matrix(df[, -1])
  
  return(mat)
} 

TGCA_mat <- convert_input_to_matrix(TGCA_mat)
inhouse_mat <- convert_input_to_matrix(inhouse_mat)

inhouse_mat <- t(inhouse_mat)
inhouse_mat <- cbind(id = 1:nrow(inhouse_mat), CellLine=rownames(inhouse_mat), 
                     as.data.frame(inhouse_mat))
inhouse_mat$CellLine <- factor(inhouse_mat$CellLine)


### unicize the sample identifiers for the tumour data
TGCA_mat <- CELLector.unicizeSamples(TGCA_mat)

CSS <- CELLector.Build_Search_Space(ctumours = t(TGCA_mat),
                                    verbose = T,
                                    minGlobSupp = 0.03,
                                    cancerType = 'Pancreatic',
                                    mutOnly=T, 
                                    UD_genomics = T)


### take all the signatures from the searching space
Signatures <- CELLector.createAllSignatures(CSS$navTable)
ModelMat <- CELLector.buildModelMatrix(Signatures$ES, inhouse_mat, CSS$navTable)
CELLector.visualiseSearchingSpace_sunBurst(CSS)

### visualising the CELLector searching space as interactive collapsible binary tree
CELLector.visualiseSearchingSpace(searchSpace = CSS, CLdata = inhouse_mat)

### selecting 10 cell lines
selectedCellLines<-CELLector.makeSelection(modelMat = ModelMat,
                                           n=12,
                                           searchSpace = CSS$navTable)
knitr::kable(selectedCellLines,align = 'l')

### Scoring cell lines based on the searching space assembled in the previous examples
CSscores <- CELLector.Score(NavTab = CSS$navTable, CELLlineData = inhouse_mat)
knitr::kable(CSscores,align = 'l')
