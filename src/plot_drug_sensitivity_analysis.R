# plot drug sensitivity
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
drug_sens_CFE <- read_csv('output/gdsctools_CFE_cmap_cl_output.csv')
ic50_drug <- read_csv('data/GDSC2_IC50.csv')
drug_ann <- read_csv('data/GDSC2_drug_annotation.csv')
GF_CS <- read_csv('data/GDSC2_GF_CS_cmap_cl.csv')
GF_CFE <- read_csv('data/GDSC2_GF_CFE_cmap_cl.csv')

###################
#### functions ####
###################

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

# volcano plot
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


#####################
volcano_plot(drug_sens_CS, "CELLector signatures", "output/plot/drug_sensitivity/gdsctools_CS_volcano_plot.png")
volcano_plot(drug_sens_CFE, "Single mutations", "output/plot/drug_sensitivity/gdsctools_CFE_volcano_plot.png")

# boxplot selected CS
drug_sens_CS_sign <- drug_sens_CS %>% 
  filter(FEATURE_IC50_T_pval < 0.005)

for(i in 1:nrow(drug_sens_CS_sign)){
  
  feature_id = drug_sens_CS_sign$FEATURE[i]
  drug_id = drug_sens_CS_sign$DRUG_ID[i]
  pl <- drug_sens_boxplot(drug_id, feature_id, GF_CS, drug_sens_CS, ic50_drug)
  pl
  ggsave(pl, filename = sprintf('output/plot/drug_sensitivity/boxplot_IC50_feature_drug_CS_row%i.png', i), width = 3, height = 5, dpi = 200)
  
}

# corresponding CFE
drug_sens_CFE_selected <- drug_sens_CFE %>% 
  filter(DRUG_ID %in% drug_sens_CS_sign$DRUG_ID, 
         FEATURE %in% c('TP53', 'CDKN2A', 'SMAD4'))

for(i in 1:nrow(drug_sens_CFE_selected)){
  
  feature_id = drug_sens_CFE_selected$FEATURE[i]
  drug_id = drug_sens_CFE_selected$DRUG_ID[i]
  pl <- drug_sens_boxplot(drug_id, feature_id, GF_CFE, drug_sens_CFE, ic50_drug)
  pl
  ggsave(pl, filename = sprintf('output/plot/drug_sensitivity/boxplot_IC50_feature_drug_CFE_row%i.png', i), width = 3, height = 5, dpi = 200)
  
}

# boxplot selected CFE
drug_sens_CFE_sign <- drug_sens_CFE %>% 
  filter(FEATURE_IC50_T_pval < 0.005)

for(i in 1:nrow(drug_sens_CS_sign)){
  
  feature_id = drug_sens_CFE_sign$FEATURE[i]
  drug_id = drug_sens_CFE_sign$DRUG_ID[i]
  pl <- drug_sens_boxplot(drug_id, feature_id, GF_CFE, drug_sens_CFE, ic50_drug)
  pl
  # ggsave(pl, filename = sprintf('output/plot/boxplot_IC50_feature_drug_CFE_sign_row%i.png', i), width = 3, height = 5, dpi = 200)
  
}


