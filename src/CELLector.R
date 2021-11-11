### load packages ###
library(devtools)
install_github("Francescojm/CELLector")
library(CELLector)
options(header = T)
options(stringsAsFactors = F)

### load input ###
TGCA_mat <- read.csv('data/tcga_snv_matrix_drivers.csv')
inhouse_mat <- read.csv('data/inhouse_snv_matrix_drivers.csv')

rownames(TGCA_mat) <- TGCA_mat$driver
TGCA_mat <- TGCA_mat %>% select(-driver)
TGCA_mat <- as.matrix(TGCA_mat)

### unicize the sample identifiers for the tumour data
tumours_BEM <- CELLector.unicizeSamples(TGCA_mat)

CSS <- CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                    verbose = T,
                                    minGlobSupp = 0.03,
                                    cancerType = 'Pancreatic',
                                    mutOnly=T, 
                                    UD_genomics = T)

### take all the signatures from the searching space
Signatures <- CELLector.createAllSignatures(CSS$navTable)
CELLector.visualiseSearchingSpace_sunBurst(CSS)


