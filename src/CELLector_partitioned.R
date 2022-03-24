### load packages ###
library(cgdsr) # Cancer Genomics Data Server
library(tidyverse)
library(dplyr)
library(data.table)
library(CELLector)
library(survival)
library(survminer)
library(knitr)
library(RColorBrewer)
library(ggrepel)



source('/group/iorio/lucia/CELLector_partitioned/R/CELLector_partioned.R')

tot_mat <- get(load('data/combined_snv_matrix_drivers.RData'))

CSS_p <- CELLector.Build_Search_Space_Partioned(ctumours = t(tot_mat),
                                    verbose = F,
                                    minGlobSupp = 0.02,
                                    cancerType = 'Pancreatic',
                                    mutOnly=T,  
                                    UD_genomics = T) 

# save space
save(CSS_p, file = 'output/CSS_partitioned_output.RData')

##########################
## clinical differences ##
##########################

create_class_from_CSS_p <- function(CSS_navTable, GF_input){
   
  Signatures <- CSS_navTable$Signature
  samples_group_list <- lapply(CSS_navTable$Points, function(x) strsplit(x, split = '[,]')[[1]])
  class_df <- data.frame(id = colnames(GF_input), group = NA, group_id = NA)

  for(i in 1:nrow(CSS_navTable)){
    class_df$group[class_df$id %in% samples_group_list[[i]]] <- CSS_navTable$Idx[i]
    class_df$group_id[class_df$id %in% samples_group_list[[i]]] <- CSS_navTable$Signature[i]
  }
  class_df$group_id <- factor(class_df$group_id, levels = CSS_navTable$Signature)
  return(class_df)
  
}

# cox model
cox_regression <- function(class_df, df_clinical, time_var, status_var, save_file = NULL){
  
  tmp <- df_clinical %>% 
     mutate(CSS = class_df$group_id) %>%
     rename(TIME = (!!time_var), STATUS = (!!status_var))
  df <- as.data.frame(tmp)
  res_cox <- coxph(Surv(TIME, STATUS) ~ CSS + SEX + AGE, data = tmp)
  
  ylab_name <- ifelse(grepl('OS_', time_var), 
                      'Overal survival probability','Recurrence probability')
  
  pl <- ggsurvplot(
        fit = survfit(Surv(TIME, STATUS) ~ CSS, data = df),
        legend = "right", 
        xlab = "Months", 
        ylab = ylab_name)
  
  cox_pvalue <- summary(res_cox)$logtest['pvalue']
    
  annotateText <- paste0('Cox log Rank P-value= ', round(cox_pvalue, 4))
  pl$plot <- pl$plot+ 
      ggplot2::annotate("text", 
                        x = Inf, y = Inf, 
                        hjust=1,vjust=1, 
                        label = annotateText, size = 4)
    
  if(!is.null(save_file)){
    ggsave(filename = save_file, plot = pl$plot, width = 10, height = 5, dpi = 200)
  }

  return(list(plot = pl, regression = res_cox))
}



### create BEM 
class_CSSp <- create_class_from_CSS_p(CSS_p$partitioned, tot_mat)
samples_id <- class_CSSp$id

# consider only mutations in CS signature
single_mut <- sapply(CSS_p$partitioned$Signature, function(x) strsplit(x, split = ',')[[1]])
single_mut <- sapply(single_mut, function(x) gsub(" ", "", x, fixed = TRUE))
single_mut <- sapply(single_mut, function(x) gsub("~", "", x, fixed = TRUE))
single_mut <- unique(unlist(single_mut))
BEM_CFE_df <- as.data.frame(t(tot_mat[single_mut, ]))


### get survival data
cgds_db <- CGDS("http://www.cbioportal.org/")
# use only paad_tcga: contains survival info, 
# also overlap with paad_qcmg_uq_2016 and paad_icgc
case_list <- getCaseLists(cgds_db,'paad_tcga')[1,1]
clinical_data <- getClinicalData(cgds_db,case_list)
new_names <- sapply(rownames(clinical_data), function(x) substr(x, 1,  nchar(x)-3))
# remove duplicates 
rm_names <- new_names[duplicated(new_names)]
if(length(rm_names)>0){
  rm_id <- unlist(lapply(rm_names, function(x) which(grepl(x, rownames(clinical_data)))))
  clinical_data <- clinical_data[-rm_id, ]    
  rownames(clinical_data) <- new_names[-rm_id]
}else{
  rownames(clinical_data) <- new_names
}
df_clinical <- clinical_data[rownames(clinical_data) %in% samples_id,]

## filter samples in input
class_CSSp <- class_CSSp[match(rownames(df_clinical), class_CSSp$id),]
BEM_CFE_df <- BEM_CFE_df[match(rownames(df_clinical),rownames(BEM_CFE_df)),]

df_clinical_subset <- df_clinical %>% 
  select(OS_STATUS, OS_MONTHS, DFS_STATUS, DFS_MONTHS, AGE, SEX) %>%
  mutate(SAMPLE_ID = rownames(df_clinical), .before = 'OS_STATUS') %>%
  mutate(SEX = as.numeric(SEX == 'Female')) %>%
  mutate(OS_STATUS_LOG = OS_STATUS == '1:DECEASED', .after = 'OS_STATUS') %>%
  mutate(DFS_STATUS_LOG = DFS_STATUS == '1:Recurred/Progressed', .after = 'DFS_STATUS')
df_clinical_subset$DFS_STATUS_LOG[df_clinical$DFS_STATUS == ''] <- NA

# cox regression
cox_OS_CS <- cox_regression(class_df =  class_CSSp, df_clinical = df_clinical_subset, 
               time_var = 'OS_MONTHS', status_var = 'OS_STATUS_LOG', 
               save_file = 'output/plot/clinical/km_curve_survival_P.png')
cox_OS_CS$plot

cox_DFS_CS <- cox_regression(class_df =  class_CSSp, df_clinical = df_clinical_subset, 
                            time_var = 'DFS_MONTHS', status_var = 'DFS_STATUS_LOG', 
                            save_file = 'output/plot/clinical/km_curve_recurrence_P.png')
cox_DFS_CS$plot


