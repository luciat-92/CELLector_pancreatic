library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

options(header = T)
options(stringsAsFactors = F)

###################
#### load data ####
###################

drug_sens_CS <- read_csv('output/gdsctools_CS_cmap_cl_output.csv')
drug_sens_CSp <- read_csv('output/gdsctools_CSpartition_cmap_cl_output.csv')
drug_sens_CFE <- read_csv('output/gdsctools_CFE_cmap_cl_output.csv')
ic50_drug <- read_csv('data/GDSC2_IC50.csv')
drug_ann <- read_csv('data/GDSC2_drug_annotation.csv')
GF_CS <- read_csv('data/GDSC2_GF_CS_cmap_cl.csv')
GF_CSp <- read_csv('data/GDSC2_GF_CSpartition_cmap_cl.csv')
GF_CFE <- read_csv('data/GDSC2_GF_CFE_cmap_cl.csv')
####################################################

volcano_plot <- function(df, title_pl, save_file = NULL){
  
  tmp <- df %>% 
    select(FEATURE,DRUG_NAME,
           FEATURE_IC50_T_pval, 
           FEATURE_IC50_effect_size, 
           FEATURE_delta_MEAN_IC50, 
           ANOVA_FEATURE_FDR) %>%
    mutate(FEATURE_IC50_effect_size_signed = sign(FEATURE_delta_MEAN_IC50)*FEATURE_IC50_effect_size) %>%
    mutate(FEATURE_IC50_T_pval_log = -log10(FEATURE_IC50_T_pval))
  
  mycolors <- rev(brewer.pal(6, "Reds"))
  
  pl <- ggplot(tmp, aes(x = FEATURE_IC50_effect_size_signed, 
                        y = FEATURE_IC50_T_pval_log, color = ANOVA_FEATURE_FDR)) + 
    geom_point() + 
    theme_bw() + 
    theme(plot.title = element_text(hjust =0.5), axis.text = element_text(size = 12))+
    ylab('-log10(pvalue)') + xlab('Signed effect size') +
    ggtitle(title_pl) + 
    scale_color_stepsn(colours = mycolors, 
                       name = 'FDR %',n.breaks=6, breaks = seq(10, 100, 15), 
                       guide = guide_coloursteps(even.steps = FALSE,
                                                 show.limits = TRUE))
  
  if(!is.null(save_file)){
    ggsave(pl, filename = save_file, width = 5, height = 5, dpi = 200)
  }
  
  return(pl)
  
}
drug_sens_boxplot <- function(drug_id, feature_id, gf_mat, df_anova, ic50_drug, save_file = NULL){
  
  df_anova_spec <- df_anova %>% 
    filter(DRUG_ID == drug_id, FEATURE == feature_id)
  
  ylab_name <- paste0(df_anova_spec$DRUG_NAME, ' [', df_anova_spec$DRUG_TARGET,']', '\nlog(IC50)')
  xlab_name <- paste0('p-value ', round(df_anova_spec$FEATURE_IC50_T_pval, 5), 
                      '\nFDR ', round(df_anova_spec$ANOVA_FEATURE_FDR,1), '%')
  
  drug_id_complete <- paste0('Drug_', drug_id,'_IC50')
  ic50_drug_spec <- ic50_drug %>% select(COSMIC_ID, (!!drug_id_complete))
  gf_mat_spec <- gf_mat %>% select(COSMIC_ID, (!!feature_id))
  
  df <- inner_join(ic50_drug_spec,gf_mat_spec) %>%
    rename(drug = (!!drug_id_complete))
  df$feature <- 'FALSE'
  df$feature[as.vector(df[,feature_id]>0)] <- 'TRUE'
  
  
  pl <- ggplot(df, aes(x = feature, 
                       y = drug)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position=position_jitter(0.2)) + 
    stat_summary(fun="mean",geom="crossbar", color = 'red',
                 mapping=aes(ymin=..y.., ymax=..y..), 
                 width=1, position=position_dodge(),show.legend = FALSE) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5), axis.text = element_text(size = 12))+
    ylab(ylab_name) + xlab(xlab_name) +
    ggtitle(feature_id)
  
  # if(!is.null(save_file)){
  # ggsave(pl, filename = save_file, width = 3, height = 5, dpi = 200)
  # }
  
  return(pl)
  
}

####################################################
volcano_plot(drug_sens_CSp, "CELLector signatures (partition)")
volcano_plot(drug_sens_CS, "CELLector signatures")

# plot old-previous pair results
drug_sens_CS$ID_comp <- paste0(drug_sens_CS$FEATURE, '_', drug_sens_CS$DRUG_ID)
drug_sens_CS$type <- 'hierarchical'
drug_sens_CSp$ID_comp <- paste0(drug_sens_CSp$FEATURE, '_', drug_sens_CSp$DRUG_ID)
drug_sens_CSp$type <- 'partition'

common_id <- intersect(drug_sens_CS$ID_comp, drug_sens_CSp$ID_comp)
drug_sens_CS <- drug_sens_CS[match(common_id, drug_sens_CS$ID_comp), ]
drug_sens_CSp <- drug_sens_CSp[match(common_id, drug_sens_CSp$ID_comp), ]

df <- data.frame(hierarchical = -log10(drug_sens_CS$FEATURE_IC50_T_pval), partition = -log10(drug_sens_CSp$FEATURE_IC50_T_pval))
pl <- ggplot(df, aes(x = hierarchical, y = partition)) + 
  geom_point(size = 1.5) + 
  geom_abline(slope=1, intercept = 0, linetype = 'dashed', color = 'red')+
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5), axis.text = element_text(size = 12))

# boxplot specific
drug_sens_boxplot('1862', 'KRASmut, TP53mut', GF_CSp, drug_sens_CSp, ic50_drug)
drug_sens_boxplot('1862', 'KRASmut, TP53mut', GF_CS, drug_sens_CS, ic50_drug)
