# survival analysis 
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


#### folder and input ####
setwd('/Volumes/home/workspace/CELLector_pancreatic/')
GF_input_file <- 'data/combined_snv_matrix_drivers.RData'
CSS_output_file <- 'output/CSS_output.RData'

CSS_output <- get(load(CSS_output_file))
GF_input <- get(load(GF_input_file))

#### functions ####
create_BEM_from_CSS <- function(CSS){
   
   Signatures <- CELLector.createAllSignatures(CSS$navTable)
   tmp <- cbind(data.frame(id = 1:ncol(GF_input), CellLine = colnames(GF_input)), t(GF_input))
   ModelMat <- CELLector.buildModelMatrix(Signatures$ES, tmp, CSS$navTable)
  
   
   return(ModelMat)
}

# cox model
cox_regression <- function(BEM_df, df_clinical, time_var, status_var){
  
  res_cox_out <- list()
  for(i in 1:ncol(BEM_df)){
    
    tmp <- df_clinical %>% 
      mutate(SIGNATURE = BEM_df[,i]) %>%
      rename(TIME = (!!time_var), STATUS = (!!status_var))
    
    res_cox <- coxph(Surv(TIME, STATUS) ~ SIGNATURE + SEX + AGE, data = tmp)
    res_cox_out[[i]] <-  as.data.frame(t(coef(summary(res_cox))[1,]))
    res_cox_out[[i]]$SIGNATURE <- colnames(BEM_df)[i]
  }
  
  res_cox_out <- do.call(rbind, res_cox_out)
  colnames(res_cox_out)[5] <- 'pvalue'
  res_cox_out <- res_cox_out %>%
    mutate(pvalue_adj = p.adjust(pvalue, method = 'BH'), .after = pvalue)
  
  return(res_cox_out)
}

plot_cox_summary <- function(cox_reg, title_pl = NULL, save_file = NULL){
  
  cox_reg <- cox_reg %>% mutate(pvalue_log = -log10(pvalue))
  cox_reg$label_sign <- cox_reg$SIGNATURE
  cox_reg$label_sign[cox_reg$pvalue_log <= 1.5] <- ''
  
  pl <- ggplot(cox_reg, aes(x = coef, 
                        y = pvalue_log, color = type, label = label_sign)) + 
    geom_point(size = 3, alpha = 0.8) + 
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'blue') + 
    geom_label_repel(min.segment.length = 0, box.padding = 0.5, 
                     force = 1) +
    theme_bw() + 
    theme(plot.title = element_text(hjust =0.5),
          axis.text = element_text(size = 12), 
          legend.title = element_blank(), 
          legend.position = 'bottom', legend.text = element_text(size = 11))+
    ylab('-log10(pvalue)') + xlab('Cox regression coefficient') +
    ggtitle(title_pl) +
    scale_color_manual(values = c('chocolate1', 'grey40'))
   
  
  if(!is.null(save_file)){
    ggsave(filename = save_file, plot = pl, width = 5, height = 5, dpi = 200)
  }
  return(pl)
  
}

plot_km_curve <- function(BEM_df, df_clinical, time_var, status_var, sign_var, 
                          cox_reg_df = NULL, save_file = NULL){
  
  df <- df_clinical %>% 
    mutate(Feature = BEM_df[,sign_var]) %>%
    rename(TIME = !!time_var, STATUS = !!status_var)
  df <- as.data.frame(df)
  
  ylab_name <- ifelse(grepl('OS_', time_var), 
                      'Overal survival probability','Recurrence probability')
  
  pl <- ggsurvplot(
        fit = survfit(Surv(TIME, STATUS) ~ Feature, data = df), 
        xlab = "Months", 
        ylab = ylab_name, 
        title = paste('Feature:', sign_var))
  
  if(!is.null(cox_reg_df)){
    id <- which(cox_reg_df$SIGNATURE == sign_var)
    cox_pvalue <- cox_reg_df$pvalue[id]
    cox_fdr <- cox_reg_df$pvalue_adj[id]
    
    annotateText <- paste0('Cox P-value= ', round(cox_pvalue, 4), 
                        '\n FDR=', round(cox_fdr*100, 1), '%')
    pl$plot <- pl$plot+ 
      ggplot2::annotate("text", 
                        x = Inf, y = Inf, 
                        hjust=1,vjust=1, 
                        label = annotateText, size = 4)
    
  }
  
  if(!is.null(save_file)){
    ggsave(filename = save_file, plot = pl$plot, width = 5, height = 5, dpi = 200)
  }
  
  return(pl)
  
}

#################
### load data ###
#################

Signatures <- CELLector.createAllSignatures(CSS$navTable)

BEM_CS <- create_BEM_from_CSS(CSS_output)
samples_id <- colnames(BEM_CS)
BEM_CS_df <- as.data.frame(t(BEM_CS))
colnames(BEM_CS_df) <- Signatures$ES

# consider only mutations in CS signature
single_mut <- sapply(Signatures$ES, function(x) strsplit(x, split = ',')[[1]])
single_mut <- sapply(single_mut, function(x) gsub(" ", "", x, fixed = TRUE))
single_mut <- sapply(single_mut, function(x) gsub("~", "", x, fixed = TRUE))
single_mut <- unique(unlist(single_mut))
BEM_CFE_df <- as.data.frame(t(GF_input[single_mut, ]))

################################################
### investigate survival based on signatures ###
################################################

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
BEM_CS_df <- BEM_CS_df[match(rownames(df_clinical),rownames(BEM_CS_df)),]
BEM_CFE_df <- BEM_CFE_df[match(rownames(df_clinical),rownames(BEM_CFE_df)),]


df_clinical_subset <- df_clinical %>% 
  select(OS_STATUS, OS_MONTHS, DFS_STATUS, DFS_MONTHS, AGE, SEX) %>%
  mutate(SAMPLE_ID = rownames(df_clinical), .before = 'OS_STATUS') %>%
  mutate(SEX = as.numeric(SEX == 'Female')) %>%
  mutate(OS_STATUS_LOG = OS_STATUS == '1:DECEASED', .after = 'OS_STATUS') %>%
  mutate(DFS_STATUS_LOG = DFS_STATUS == '1:Recurred/Progressed', .after = 'DFS_STATUS')
df_clinical_subset$DFS_STATUS_LOG[df_clinical$DFS_STATUS == ''] <- NA
# save 
write.table('data/paad_tcga_clinical.csv', x = df_clinical_subset, quote = T, row.names = F, 
            col.names = T, sep = ',')

# cox regression
cox_OS_CS <- cox_regression(BEM_df = BEM_CS_df, df_clinical = df_clinical_subset, 
               time_var = 'OS_MONTHS', status_var = 'OS_STATUS_LOG')

cox_DFS_CS <- cox_regression(BEM_df = BEM_CS_df, df_clinical = df_clinical_subset, 
                            time_var = 'DFS_MONTHS', status_var = 'DFS_STATUS_LOG')

cox_OS_CFE <- cox_regression(BEM_df = BEM_CFE_df, df_clinical = df_clinical_subset, 
                            time_var = 'OS_MONTHS', status_var = 'OS_STATUS_LOG')
kable(cox_OS_CFE)

cox_DFS_CFE <- cox_regression(BEM_df = BEM_CFE_df, df_clinical = df_clinical_subset, 
                             time_var = 'DFS_MONTHS', status_var = 'DFS_STATUS_LOG')
kable(cox_OS_CS)

# plot summary regression combining signatures and CS
cox_OS_CS <- cox_OS_CS %>% mutate(type = 'CELLector signatures')
cox_OS_CFE <- cox_OS_CFE %>% mutate(type = 'Single mutations')
cox_reg_OS <- rbind(cox_OS_CS, cox_OS_CFE)
# save
write.table(file = 'output/cox_regression_survival_summary.tsv', x = cox_reg_OS, quote = T, 
            row.names = F, sep = '\t', col.names = T)
plot_cox_summary(cox_reg = cox_reg_OS,  title_pl = 'Overal Survival', 
                 save_file = 'output/plot/clinical/cox_regression_survival_summary.png')

cox_DFS_CS <- cox_DFS_CS %>% mutate(type = 'CELLector signatures')
cox_DFS_CFE <- cox_DFS_CFE %>% mutate(type = 'Single mutations')
cox_reg_DFS <- rbind(cox_DFS_CFE, cox_DFS_CS)
# save
write.table(file = 'output/cox_regression_recurrence_summary.tsv', x = cox_reg_DFS, quote = T, 
            row.names = F, sep = '\t', col.names = T)
plot_cox_summary(cox_reg = cox_reg_DFS, title_pl = 'Recurrence', 
                 save_file = 'output/plot/clinical/cox_regression_recurrence_summary.png')


# K-M plot for selected results
# Surival
plot_km_curve(BEM_df = BEM_CS_df, df_clinical = df_clinical_subset, 
              time_var =  'OS_MONTHS', status_var = 'OS_STATUS_LOG', 
              sign_var = 'KRAS, TP53', cox_reg_df = cox_OS_CS, 
              save_file = 'output/plot/clinical/km_curve_survival_KRASTP53.png')

plot_km_curve(BEM_df = BEM_CFE_df, df_clinical = df_clinical_subset, 
              time_var =  'OS_MONTHS', status_var = 'OS_STATUS_LOG', 
              sign_var = 'KRAS', cox_reg_df = cox_OS_CFE, 
              save_file = 'output/plot/clinical/km_curve_survival_KRAS.png')

plot_km_curve(BEM_df = BEM_CFE_df, df_clinical = df_clinical_subset, 
              time_var =  'OS_MONTHS', status_var = 'OS_STATUS_LOG', 
              sign_var = 'TP53',cox_reg_df = cox_OS_CFE, 
              save_file = 'output/plot/clinical/km_curve_survival_TP53.png')

# Disease free
plot_km_curve(BEM_df = BEM_CS_df, df_clinical = df_clinical_subset, 
              time_var =  'DFS_MONTHS', status_var = 'DFS_STATUS_LOG', 
              sign_var = 'KRAS, TP53', cox_reg_df = cox_DFS_CS, 
              save_file = 'output/plot/clinical/km_curve_recurrence_KRASTP53.png')

plot_km_curve(BEM_df = BEM_CFE_df, df_clinical = df_clinical_subset, 
              time_var =  'DFS_MONTHS', status_var = 'DFS_STATUS_LOG', 
              sign_var = 'KRAS', cox_reg_df = cox_DFS_CFE, 
              save_file = 'output/plot/clinical/km_curve_recurrence_KRAS.png')

plot_km_curve(BEM_df = BEM_CFE_df, df_clinical = df_clinical_subset, 
              time_var =  'DFS_MONTHS', status_var = 'DFS_STATUS_LOG', 
              sign_var = 'TP53', cox_reg_df = cox_DFS_CFE, 
              save_file = 'output/plot/clinical/km_curve_recurrence_TP53.png')







