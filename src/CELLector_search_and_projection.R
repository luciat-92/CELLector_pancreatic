### load packages ###
library(CELLector)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridisLite)

options(header = T)
options(stringsAsFactors = F)

#### folder and input ####
setwd('/Volumes/home/workspace/CELLector_pancreatic/')
cmap_mutation <- '/Volumes/iorio/CellModelPassports/genes/mutations_current.csv'
cmap_model <- '/Volumes/iorio/CellModelPassports/models/model_list_20210611.csv'

#################
### functions ###
#################

plot_frequency <- function(mat, lines = F, count_thr = 1, title, save_file = NULL){
  
  # count mutations on organoid lines
  df_freq <- data.frame(mut = rownames(mat), freq = 100*rowSums(mat)/ncol(mat), 
                        count = rowSums(mat))
  df_freq$mut <- factor(df_freq$mut, levels = df_freq$mut[order(rowSums(mat),decreasing=TRUE)])
  ylab_name <- ifelse(lines,  'n. lines', 'n. patients')
  
  pl <- ggplot(subset(df_freq, count > count_thr), aes(x = mut, y = count)) + 
    geom_bar(stat = 'identity') + 
    theme_classic() + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust =0.5))+
    ylab(ylab_name) + 
    ggtitle(title)+
    scale_y_continuous(ylab_name, sec.axis = sec_axis(~ . * 100/ncol(mat), name = "%"))
  
  if(!is.null(save_file)){
    ggsave(plot = pl, width = 4, height = 6, dpi = 200,
           filename = save_file)
  }
  
  return(pl)
  
}

plot_frequency_line <- function(mat, lines = F, count_thr = 1, title, save_file = NULL){
  
  # count mutations on organoid lines
  df_freq <- data.frame(mut = rownames(mat), freq = 100*rowSums(mat)/ncol(mat), 
                        count = rowSums(mat))
  df_freq$mut <- factor(df_freq$mut, levels = df_freq$mut[order(rowSums(mat),decreasing=TRUE)])
  ylab_name <- 'n. mutations'
  
  pl <- ggplot(subset(df_freq, count > count_thr), aes(x = mut, y = count)) + 
    geom_bar(stat = 'identity') + 
    theme_classic() + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust =0.5))+
    ylab(ylab_name) + 
    ggtitle(title)+
    scale_y_continuous(ylab_name, sec.axis = sec_axis(~ . * 100/ncol(mat), name = "%"))
  
  if(!is.null(save_file)){
    ggsave(plot = pl, width = 4, height = 6, dpi = 200,
           filename = save_file)
  }
  
  return(pl)
  
}

convert_input_to_matrix <- function(df){
  
  rownames(df) <- df[,1]
  mat <- as.matrix(df[, -1])

  return(mat)
} 

summary_CSS <- function(Signatures, ModelMat, n_patients){
  
  # Signature: output of CELLector.createAllSignatures
  # ModelMat: output of CELLector.buildModelMatrix
  # n_patients : number of patients for search space
  
  df <- data.frame(Subtype = names(Signatures$ES), 
                   Signatures = Signatures$ES, 
                   Signatures_complete = Signatures$S,
                   N.patients = round(n_patients*Signatures$STS/100), 
                   P.patients = Signatures$STS, 
                   n_CL = NA, repr_CL = NA)
  
  id_mapped <- rownames(ModelMat)
  id_lines <- colnames(ModelMat)
  
  for(id_row in 1:nrow(df)){
    
    if(!df$Subtype[id_row] %in% id_mapped){
      df$repr_CL[id_row] <- 'lack of in vitro models'
      df$n_CL[id_row] <- 0
    }else{
      id_row_model <- which(id_mapped == df$Subtype[id_row])
      df$n_CL[id_row] <- sum(ModelMat[id_row_model, ])
      df$repr_CL[id_row] <- paste0(id_lines[ModelMat[id_row_model, ] == 1], collapse = ',')  
    }
  }
  return(df)
}

coverage_CSS <- function(navTable, ModelMat = NULL){
  
  if(!is.null(ModelMat)){
    navTable <- navTable[navTable$Idx %in% rownames(ModelMat),]
  }
  n_samples_cov <- navTable$AbsSupport[1]
  right_child <- navTable$Right.Child.Index[1]
  stop_cond <- right_child == 0
  
  while(!stop_cond){
    n_samples_cov <- n_samples_cov + navTable$AbsSupport[right_child]
    right_child <- navTable$Right.Child.Index[right_child]
    stop_cond <- right_child == 0
  }
  
  return(n_samples_cov)
  
}

coverage_CSS_max_granularity <- function(navTable, subset_proj){
  
  n_samples_projected_signatures <- data.frame(Subtype = subset_proj$Subtype, 
                                               Signatures = subset_proj$Signatures, 
                                               n_samples = NA)
  
  for(id_row in 1:nrow(subset_proj)){
    
    n <- subset_proj$N.patients[id_row]
    minus_n <- 0
    id_info <- which(subset_proj$Subtype[id_row] == navTable$Idx)
    
    # check if there is a left child
    if(navTable$Left.Child.Index[id_info]!=0){
      
      new_id <- navTable$Left.Child.Index[id_info]
      minus_n <- navTable$AbsSupport[navTable$Idx == new_id]
      
      # add all the right nodes
      exists_right <- navTable$Right.Child.Index[navTable$Idx == new_id] != 0
      
      while(exists_right){
        new_id <- navTable$Right.Child.Index[navTable$Idx == new_id]
        minus_n <- minus_n + navTable$AbsSupport[navTable$Idx == new_id]
        exists_right <- navTable$Right.Child.Index[navTable$Idx == new_id] != 0
      }
    }
    
    n_samples_projected_signatures$n_samples[id_row] <- n-minus_n
  }
  
  
  return(n_samples_projected_signatures)
  
}

line_signature_heatmap <- function(df, title_pl){
  
  signature_cell_lines <- unique(df$Signature)
  
  matrix_assign <- matrix(ncol = length(signature_cell_lines), nrow = nrow(df))
  colnames(matrix_assign) <- signature_cell_lines
  rownames(matrix_assign) <- df$CellLines
  
  for(id_col in 1:length(signature_cell_lines)){
    names_line <- df$CellLines[df$Signature == signature_cell_lines[id_col]]
    matrix_assign[names_line,id_col] <- df$CELLectorScores[df$Signature == signature_cell_lines[id_col]]
  }
  
  heatmap_score <- pheatmap(matrix_assign, cluster_rows = F, cluster_cols = F, 
                            main = title_pl, col = viridis(10), 
                            breaks = seq(0, 1, by = 0.1))
  
  return(heatmap_score)
  
}

CELLline_buildBEM_annotationProvided <- function(annCat = NULL,  
                                                varCat = NULL, Tissue, Cancer_Type, 
                                                Cancer_Type_details = NULL,sample_site = NULL, 
                                                excludeOrganoids = FALSE, humanonly = TRUE, 
                                                msi_status_select = NULL, gender_select = NULL, 
                                                mutational_burden_th = NULL, age_at_sampling = NULL, 
                                                ploidy_th = NULL, ethnicity_to_exclude = NULL, 
                                                GenesToConsider = NULL, VariantsToConsider = NULL){
  
  if (length(annCat) == 0){
    clAnnotation <- CELLector.CMPs_getModelAnnotation()
  }else{
    clAnnotation <- annCat
  }
  
  clAnnotation$cancer_type_detail <- str_sub(clAnnotation$cancer_type_detail, 
                                            3, end = str_length(clAnnotation$cancer_type_detail) - 3)
  
  if (length(varCat) == 0) {
    varCat <- CELLector.CMPs_getVariants()
  }else{  
    if (!excludeOrganoids) {
      id <- which(clAnnotation$tissue == Tissue & is.element(clAnnotation$cancer_type, 
                                                             Cancer_Type))
    }
    else {
      id <- which(clAnnotation$tissue == Tissue & is.element(clAnnotation$cancer_type, 
                                                             Cancer_Type) & clAnnotation$model_type != "Organoid")
    }
    cls <- clAnnotation$model_id[id]
    varCat <- varCat[which(is.element(varCat$model_id, cls)), 
    ]
    clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                  cls)), ]
    if (length(Cancer_Type_details) > 0) {
      id <- which(is.element(clAnnotation$cancer_type_detail, 
                             Cancer_Type_details))
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(sample_site) > 0) {
      id <- which(is.element(clAnnotation$sample_site, 
                             sample_site))
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(humanonly) > 0) {
      id <- which(clAnnotation$species == "Homo Sapiens")
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(msi_status_select) > 0) {
      id <- which(!is.na(clAnnotation$msi_status) & (clAnnotation$msi_status == 
                                                       msi_status_select | (msi_status_select == "MSI-L/H" & 
                                                                              (clAnnotation$msi_status == "MSI-L" | clAnnotation$msi_status == 
                                                                                 "MSI-H"))))
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(gender_select) > 0) {
      id <- which(is.element(clAnnotation$gender, gender_select))
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(ethnicity_to_exclude) > 0) {
      id <- which(!is.element(clAnnotation$ethnicity, ethnicity_to_exclude))
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(mutational_burden_th) > 0) {
      id <- which(round(clAnnotation$mutational_burden) >= 
                    mutational_burden_th[1] & round(clAnnotation$mutational_burden) <= 
                    mutational_burden_th[2])
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(ploidy_th) > 0) {
      id <- which(round(clAnnotation$ploidy) >= ploidy_th[1] & 
                    round(clAnnotation$ploidy) <= ploidy_th[2])
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
    if (length(age_at_sampling) > 0) {
      id <- which(round(clAnnotation$age_at_sampling) >= 
                    age_at_sampling[1] & round(clAnnotation$age_at_sampling) <= 
                    age_at_sampling[2])
      cls <- clAnnotation$model_id[id]
      varCat <- varCat[which(is.element(varCat$model_id, 
                                        cls)), ]
      clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id, 
                                                    cls)), ]
    }
  }
  if (length(GenesToConsider) > 0) {
    varCat <- varCat[which(is.element(varCat$gene_symbol, 
                                      GenesToConsider)), ]
  }
  if (length(VariantsToConsider) > 0) {
    sigs <- paste(varCat$gene_symbol, varCat$cdna_mutation, 
                  paste("p.", varCat$aa_mutation, sep = ""))
    varCat <- varCat[which(is.element(sigs, VariantsToConsider)), 
    ]
  }
  allModels <- sort(unique(varCat$model_id))
  allModel_ids <- varCat$model_id[match(allModels, varCat$model_id)]
  allGenes <- sort(unique(varCat$gene_symbol))
  BEM <- do.call(what = cbind, lapply(allModels, function(x) {
    is.element(allGenes, varCat$gene_symbol[varCat$model_id == 
                                              x]) + 0
  }))
  rownames(BEM) <- allGenes
  cls <- clAnnotation$model_name[match(allModel_ids, clAnnotation$model_id)]
  BEM <- data.frame(CMP_identifier = allModel_ids, CellLine = cls, 
                    t(BEM), stringsAsFactors = FALSE)
  return(BEM)
}

project_CSS <- function(CSS, mat_to_project, n_samples_CSS){
  
  Signatures <- CELLector.createAllSignatures(CSS$navTable)
  ModelMat <- CELLector.buildModelMatrix(Signatures$ES, mat_to_project, CSS$navTable)
  # score projected lines
  CSscores <- CELLector.Score(NavTab = CSS$navTable, CELLlineData = mat_to_project)
  CSscores_tab <- CSscores[!duplicated(CSscores$CellLines),]
  # get summary space:
  summary_res <- summary_CSS(Signatures = Signatures, 
                             ModelMat = ModelMat, 
                             n_patients = n_samples_CSS)
  
  # get percentage coverage
  subset_proj <- summary_res %>% filter(Signatures_complete %in% unique(CSscores_tab$Signature))
  n_samples_mapped_on_projected_maxg <- coverage_CSS_max_granularity(navTable = CSS$navTable, 
                                                                  subset_proj = subset_proj)
  n_samples_mapped_on_projected <- coverage_CSS(navTable = CSS$navTable, 
                                                ModelMat = ModelMat)
  
  return(list(ModelMat = ModelMat, 
              CS_score = CSscores_tab, 
              summary = summary_res, 
              n_samples_persign_max_gran = n_samples_mapped_on_projected_maxg, 
              n_samples_persign = n_samples_mapped_on_projected))
  
}

pie_chart_max_granularity <- function(n_samples, n_samples_mapped, 
                                      n_samples_mapped_with_lines, 
                                      save_file = NULL){
  
  # pie chart mapping
  missing <- n_samples - n_samples_mapped
  missing_in_lines <- n_samples_mapped - sum(n_samples_mapped_with_lines$n_samples)
  
  df_coverage <- data.frame(n = c(n_samples_mapped_with_lines$n_samples, 
                                  missing_in_lines, missing), 
                            type = c(n_samples_mapped_with_lines$Signatures, 
                                     'without organoid model', 
                                     'not in CELLector search space'))
  df_coverage <- df_coverage %>% 
    mutate(prop = n/n_samples*100, 
           prop_lab =  paste0(round(n/n_samples*100, digits = 1), '%'))
  df_coverage$type <- factor(df_coverage$type, levels = rev(df_coverage$type))

  nb.cols <- nrow(df_coverage) - 2
  mycolors <- c('grey90','grey50',rev(colorRampPalette(brewer.pal(9, "Oranges"))(nb.cols)))
  
  pie_chart <- ggplot(df_coverage, aes(x = "", y = prop, fill = type)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0)+
    geom_label_repel(aes(label = prop_lab), 
                     position = position_stack(vjust=0.5), 
                     min.segment.length = 0, 
                     force = 5) +
    scale_fill_manual(values = mycolors) +
    theme_void() + 
    theme(legend.title = element_blank())

  if(!is.null(save_file)){
    ggsave(pie_chart, filename = save_file, width = 5, height = 5, dpi = 200)
  }
  
  return(pie_chart)
  
}


pie_chart_mapped <- function(n_samples, n_samples_mapped, 
                                      n_samples_persign,
                                      save_file = NULL){
  
  # pie chart mapping
  missing <- n_samples - n_samples_mapped
  missing_in_lines <- n_samples_mapped - n_samples_persign
  
  df_coverage <- data.frame(n = c(n_samples_mapped, missing_in_lines, missing), 
                            type = c('in CELLector search space', 
                                     'without organoid model', 
                                     'not in CELLector search space'))
  
  df_coverage <- df_coverage %>% 
    mutate(prop = n/n_samples*100, 
           prop_lab =  paste0(round(n/n_samples*100, digits = 1), '%'))
  df_coverage$type <- factor(df_coverage$type, levels = rev(df_coverage$type))
  mycols <- c("grey90","grey50","orange")

  pie_chart <- ggplot(df_coverage, aes(x = "", y = prop, fill = type)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0)+
    geom_label_repel(aes(label = prop_lab), 
                     position = position_stack(vjust=0.5), 
                     min.segment.length = 0, 
                     force = 5) +
    scale_fill_manual(values = mycols) +
    theme_void()+
    theme(legend.title = element_blank())
  
  if(!is.null(save_file)){
    ggsave(pie_chart, filename = save_file, width = 5, height = 5, dpi = 200)
  }
  
  return(pie_chart)
  
}

##################
### load input ###
##################

TGCA_mat <- read.csv('data/tcga_snv_matrix_drivers.csv')
icgc_mat <- read.csv('data/icgc_snv_matrix_drivers.csv')
bailey_mat <- read.csv('data/bailey_snv_matrix_drivers.csv')
inhouse_mat <- read.csv('data/inhouse_snv_matrix_drivers.csv')

TGCA_mat <- convert_input_to_matrix(TGCA_mat)
icgc_mat <- convert_input_to_matrix(icgc_mat)
bailey_mat <- convert_input_to_matrix(bailey_mat)
inhouse_mat <- convert_input_to_matrix(inhouse_mat)
inhouse_mat <- inhouse_mat[rowSums(inhouse_mat)>0,]
  
### combine all model tables ###
common_mut <- intersect(intersect(rownames(TGCA_mat), rownames(icgc_mat)), rownames(bailey_mat))
tot_mat <- cbind(TGCA_mat[common_mut,], icgc_mat[common_mut,], bailey_mat[common_mut,])
tot_mat <- tot_mat[rowSums(tot_mat)>0,]
# save .RData format
res <- tot_mat
save(res, file = 'data/combined_snv_matrix_drivers.RData')

## load cell lines data and save matrices
modelAnnotation <- read_csv(cmap_model)
modelMutations <- read_csv(cmap_mutation)

# Note: impossible to retrieve organoid data, only 1 line for pancreatic cancer
PA_cl_BEM <- CELLline_buildBEM_annotationProvided(
  annCat = modelAnnotation,
  varCat = modelMutations,
  Tissue='Pancreas',
  Cancer_Type = 'Pancreatic Carcinoma',
  GenesToConsider = common_mut,
  excludeOrganoids = T)
PA_cl_BEM_mat <- t(convert_input_to_matrix(PA_cl_BEM[, -1]))
# save
save(PA_cl_BEM, file = 'data/cmap_PA_cl_snv_matrix_drivers.RData')

##################################
#### plot frequency mutations ####
##################################

# count mutations combined model table
model_mat <- list(tgca = TGCA_mat, icgc = icgc_mat, bailey = bailey_mat)
df_freq <- data.frame(mut = rep(common_mut, 3),
                      count =  unlist(lapply(model_mat, function(x) rowSums(x[common_mut,]))), 
                      set = unlist(lapply(names(model_mat), function(x) rep(x, length(common_mut)))))
df_freq$mut <- factor(df_freq$mut, levels = common_mut[order(rowSums(tot_mat),decreasing=TRUE)])

subset_mut <- df_freq %>% group_by(mut) %>% 
  summarise(n = sum(count)) %>% filter(n > 10)

pl <- ggplot(subset(df_freq, mut %in% subset_mut$mut), aes(x = mut, y = count, fill = set)) + 
  geom_bar(stat = 'identity') + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust =0.5))+
  ylab('count patients') + 
  scale_y_continuous('count patients', sec.axis = sec_axis(~ . * 100/ncol(tot_mat), name = "percentage"))+
  ggtitle('combined')
pl
ggsave(plot = pl, width = 7, height = 6, dpi = 200,
       filename = 'output/plot/combined_primary_frequency_mutations.png')

# plot distribution number of mutation per sample
df_nmut_per_sample <- data.frame(name = colnames(tot_mat),label = colnames(tot_mat), 
                                 n = colSums(tot_mat), lab = 'tot')
df_nmut_per_sample$label[df_nmut_per_sample$n<10] <- ""
pl <- ggplot(df_nmut_per_sample, aes(x=lab, y=n, label = label)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.4, alpha = 0.5) +
  geom_text(size = 3)+
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        plot.title = element_text(hjust =0.5))+
  ylab('n. mutations') + 
  ggtitle('combined')
pl
ggsave(plot = pl, width = 2, height = 6, dpi = 200,
       filename = 'output/plot/combined_primary_frequency_samples.png')


# organoid data to project space on
plot_frequency(inhouse_mat,  count_thr = 1, lines = T, title = 'inhouse', 
               save_file = 'output/plot/inhouse_frequency_mutations.png')

plot_frequency_line(t(inhouse_mat),  count_thr = 1, lines = T, title = 'inhouse', 
                    save_file = 'output/plot/inhouse_frequency_lines.png')

# cell lines to project space on
plot_frequency(PA_cl_BEM_mat,  count_thr = 1, lines = T, title = 'cmap cell lines', 
               save_file = 'output/plot/cmap_cl_frequency_mutations.png')

plot_frequency_line(t(PA_cl_BEM_mat),  count_thr = 1, lines = T, title = 'cmap cell lines', 
                    save_file = 'output/plot/cmap_cl_frequency_lines.png')

################################
#### create CELLector space ####
################################

inhouse_mat <- t(inhouse_mat)
inhouse_mat <- cbind(id = 1:nrow(inhouse_mat), CellLine=rownames(inhouse_mat), 
                     as.data.frame(inhouse_mat))
inhouse_mat$CellLine <- factor(inhouse_mat$CellLine)
res = inhouse_mat
save(res, file = 'data/inhouse_snv_matrix_drivers.RData')

tot_mat <- CELLector.unicizeSamples(tot_mat)

CSS <- CELLector.Build_Search_Space(ctumours = t(tot_mat),
                                    verbose = F,
                                    minGlobSupp = 0.02,
                                    cancerType = 'Pancreatic',
                                    mutOnly=T,  
                                    UD_genomics = T) 
# save space
save(CSS, file = 'output/CSS_output.RData')

Signatures <- CELLector.createAllSignatures(CSS$navTable)
Signatures_tab <- data.frame(S = Signatures$S, Perc = Signatures$STS)
n_samples_mapped <- coverage_CSS(navTable = CSS$navTable)

### plot sunburst
CELLector.visualiseSearchingSpace_sunBurst(CSS)

############################
### project on organoids ###
############################

project_CSS_inhouse <- project_CSS(CSS, inhouse_mat, n_samples_CSS = ncol(tot_mat))
write.table(file = 'output/summary_CSS_projected_inhouse.tsv', project_CSS_inhouse$summary, sep = '\t', 
            quote = F, col.names = T, row.names = F)
# plot tree
CELLector.visualiseSearchingSpace(searchSpace = CSS, CLdata = inhouse_mat)
# heatmap line scores
hm_plot <- line_signature_heatmap(project_CSS_inhouse$CS_score, title_pl = 'inhouse')

save_file = sprintf('output/plot/organoids_score_projected.png')
png(filename = save_file, width = 4, height = 6, units = 'in', res = 300)
hm_plot
dev.off()

# save output
project_CSS_obj <- project_CSS_inhouse
save(project_CSS_obj, file = 'output/project_CSS_inhouse.RData')

# create BEM for signatures
id <- match(rownames(project_CSS_inhouse$ModelMat), project_CSS_inhouse$summary$Subtype)
project_CSS_BEM <- cbind(data.frame(Signature = project_CSS_inhouse$summary$Signatures_complete[id]), 
                         as.data.frame(project_CSS_inhouse$ModelMat))
write.table(file = 'output/project_CSS_inhouse_BEM.tsv', x = project_CSS_BEM, 
            quote = F, col.names = T, row.names = F, sep = '\t')

# pie chart plots
pie_chart_max_granularity(n_samples = ncol(tot_mat), 
                          n_samples_mapped = n_samples_mapped, 
                          n_samples_mapped_with_lines = project_CSS_inhouse$n_samples_persign_max_gran, 
                          save_file = 'output/plot/piechart_coverage_inhouse_max_granularity.png')

pie_chart_mapped(n_samples = ncol(tot_mat), n_samples_mapped = n_samples_mapped, 
                 n_samples_persign = project_CSS_inhouse$n_samples_persign, 
                 save_file = 'output/plot/piechart_coverage_inhouse.png')


##########################
### project on cmap cl ###
##########################

project_CSS_cmap_cl <- project_CSS(CSS, PA_cl_BEM, n_samples_CSS = ncol(tot_mat))
# save summary
write.table(file = 'output/summary_CSS_projected_cmap_cl.tsv', project_CSS_cmap_cl$summary, sep = '\t', 
            quote = F, col.names = T, row.names = F)
# plot tree
CELLector.visualiseSearchingSpace(searchSpace = CSS, CLdata = PA_cl_BEM)
# heatmap line scores
hm_plot <- line_signature_heatmap(project_CSS_cmap_cl$CS_score, title_pl = 'cmap_cl')

save_file = sprintf('output/plot/cmap_cl_score_projected.png')
png(filename = save_file, width = 4, height = 8, units = 'in', res = 300)
hm_plot
dev.off()

# save output
project_CSS_obj <- project_CSS_cmap_cl
save(project_CSS_obj, file = 'output/project_CSS_cmap_cl.RData')

# create BEM for signatures
id <- match(rownames(project_CSS_cmap_cl$ModelMat), project_CSS_cmap_cl$summary$Subtype)
project_CSS_BEM <- cbind(data.frame(Signature = project_CSS_cmap_cl$summary$Signatures_complete[id]), 
                         as.data.frame(project_CSS_cmap_cl$ModelMat))
write.table(file = 'output/project_CSS_cmap_cl_BEM.tsv', x = project_CSS_BEM, 
            quote = F, col.names = T, row.names = F, sep = '\t')

# pie chart plots
pie_chart_max_granularity(n_samples = ncol(tot_mat), n_samples_mapped = n_samples_mapped, 
                          n_samples_mapped_with_lines = project_CSS_cmap_cl$n_samples_persign_max_gran, 
                          save_file = 'output/plot/piechart_coverage_cmap_cl_max_granularity.png')

pie_chart_mapped(n_samples = ncol(tot_mat), n_samples_mapped = n_samples_mapped, 
                 n_samples_persign = project_CSS_cmap_cl$n_samples_persign, 
                 save_file = 'output/plot/piechart_coverage_cmap_cl.png')



