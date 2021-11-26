# format data table to perform drug sensitivity analysis
library(readxl)
library(tidyverse)
library(dplyr)
library(data.table)

#### folder and input ####
setwd('/Volumes/home/workspace/CELLector_pancreatic/')
cmap_drug_file <- '/Volumes/iorio/CellModelPassports/drugs/GDSC2_fitted_dose_response_25Feb20.xlsx'
CSS_cmap_cl_BEM_file <- 'output/project_CSS_cmap_cl_BEM.tsv'
cmap_cl_BEM_file <- 'data/cmap_PA_cl_snv_matrix_drivers.RData'
cmap_model_file <- '/Volumes/iorio/CellModelPassports/models/model_list_20210611.csv'

##########################################
## convert IC50 in format for gdscTools ##
##########################################

cmap_drug <- read_excel(cmap_drug_file, sheet = 1)
cosmic_id <- unique(cmap_drug$COSMIC_ID)
drug_id <- unique(cmap_drug$DRUG_ID)

ic50_drug_mat <- matrix(NA, ncol = length(drug_id), nrow = length(cosmic_id))
colnames(ic50_drug_mat) <- paste0('Drug_', drug_id, '_IC50')

for(id_col in 1:length(drug_id)){
  
  tmp <- cmap_drug %>% filter(DRUG_ID == drug_id[id_col])
  
  common_cosmic_id <- intersect(tmp$COSMIC_ID, cosmic_id)
  id_row <- match(common_cosmic_id, cosmic_id)
  
  ic50_drug_mat[id_row, id_col] <- tmp$LN_IC50[match(common_cosmic_id, tmp$COSMIC_ID)]
  
}

ic50_drug <- cbind(data.frame(COSMIC_ID = cosmic_id), ic50_drug_mat)

# save
write.table(file = 'data/GDSC2_IC50.csv', x = ic50_drug, quote = F, sep = ',', 
            col.names = T, row.names = F)

#####################
## drug annotation ##
#####################

drug_ann <-  cmap_drug %>% 
  distinct(DRUG_ID, .keep_all= TRUE) %>% 
  select(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET) %>%
  rename(DRUG_TARGET = PUTATIVE_TARGET)

# save
write.table(file = 'data/GDSC2_drug_annotation.csv', x = drug_ann, quote = T, sep = ',', 
            col.names = T, row.names = F)


#########################################
## convert BEM in format for gdscTools ##
#########################################

cmap_model <- read_csv(cmap_model_file)
cmap_cl_BEM <- get(load(cmap_cl_BEM_file))
CSS_cmap_cl_BEM <- fread(CSS_cmap_cl_BEM_file, h=T, sep = '\t')
#
cl_id <- colnames(CSS_cmap_cl_BEM)[-1]
cosmic_id_cl <- cmap_model[match(cl_id, cmap_model$model_name),]

cl_BEM <-  cbind(data.frame(COSMIC_ID = cosmic_id_cl$COSMIC_ID, 
                      TISSUE_FACTOR = cosmic_id_cl$tissue, 
                      MSI_FACTOR = as.numeric(cosmic_id_cl$msi_status == "MSS")), 
                 t(CSS_cmap_cl_BEM[, -1]))
cl_BEM$MSI_FACTOR[is.na(cl_BEM$MSI_FACTOR)] <- 0
colnames(cl_BEM)[-(1:3)] <- CSS_cmap_cl_BEM$Signature

# # exclude MSI factor if less than 2 
# if(sum(cl_BEM$MSI_FACTOR == 0) < 2 | sum(cl_BEM$MSI_FACTOR == 1) < 2){
#   cl_BEM <- cl_BEM %>% select(-MSI_FACTOR)
# }
# # exclude signature if less than 2
# id_rm <- which(sapply(CSS_cmap_cl_BEM$Signature, function(x)
#                       sum(cl_BEM[, x]== 0)< 2 | sum(cl_BEM[,x]==1) < 2))
# sign_rm <- names(id_rm)
# cl_BEM <- cl_BEM %>% select(-(!!sign_rm))

# save
write.table(file = 'data/GDSC2_GF_CS_cmap_cl.csv', x = cl_BEM, quote = T, sep = ',', 
            col.names = T, row.names = F)

# 
# single mutations in CS
single_mut <- sapply(CSS_cmap_cl_BEM$Signature, function(x) strsplit(x, split = ',')[[1]])
single_mut <- sapply(single_mut, function(x) gsub(" ", "", x, fixed = TRUE))
single_mut <- sapply(single_mut, function(x) gsub("~", "", x, fixed = TRUE))
single_mut <- sapply(single_mut, function(x) gsub("mut", "", x, fixed = TRUE))
single_mut <- unique(unlist(single_mut))

cl_id <- cmap_cl_BEM$CellLine
cosmic_id_cl <- cmap_model[match(cl_id, cmap_model$model_name),]

cl_BEM <-  cbind(data.frame(COSMIC_ID = cosmic_id_cl$COSMIC_ID, 
                            TISSUE_FACTOR = cosmic_id_cl$tissue, 
                            MSI_FACTOR = as.numeric(cosmic_id_cl$msi_status == "MSS")), 
                 cmap_cl_BEM[, single_mut])
cl_BEM$MSI_FACTOR[is.na(cl_BEM$MSI_FACTOR)] <- 0
write.table(file = 'data/GDSC2_GF_CFE_cmap_cl.csv', x = cl_BEM, quote = T, sep = ',', 
            col.names = T, row.names = F)
