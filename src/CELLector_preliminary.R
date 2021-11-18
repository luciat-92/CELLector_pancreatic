### load packages ###
library(devtools)
library(CELLector)
library(ggplot2)
library(pheatmap)
library(viridisLite)

options(header = T)
options(stringsAsFactors = F)

setwd('/Volumes/home/workspace/CELLector_pancreatic/')

## functions 
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


filter_mutation <- function(mat, mutation_name){
  
  index_keep <- mat[mutation_name,] == 1
  filtered_mat <- mat[, index_keep]
  filtered_mat <- filtered_mat[!rownames(filtered_mat) %in% mutation_name, ]
  
  return(filtered_mat)
}

convert_input_to_matrix <- function(df){
  
  rownames(df) <- df[,1]
  mat <- as.matrix(df[, -1])
  
  return(mat)
} 

search_space_and_project <- function(mat, mat_to_project, minGlobSupp = 0.02){
  
  # unicize the sample identifiers for the tumour data
  mat <- CELLector.unicizeSamples(mat)
  
  CSS <- CELLector.Build_Search_Space(ctumours = t(mat),
                                      verbose = F,
                                      minGlobSupp = minGlobSupp,
                                      cancerType = 'Pancreatic',
                                      mutOnly=T,  
                                      UD_genomics = T) 
  
  Signatures <- CELLector.createAllSignatures(CSS$navTable)
  ModelMat <- CELLector.buildModelMatrix(Signatures$ES, mat_to_project, CSS$navTable)
  Signatures_tab <- data.frame(S = Signatures$S, Perc = Signatures$STS)
  
  # visualize
  space_plot <- CELLector.visualiseSearchingSpace_sunBurst(CSS)
  CELLector.visualiseSearchingSpace(searchSpace = CSS, CLdata = mat_to_project)
  
  ### selecting all lines
  selectedCellLines <- CELLector.makeSelection(modelMat = ModelMat,
                                               n=ncol(ModelMat),
                                               searchSpace = CSS$navTable)
  #knitr::kable(selectedCellLines,align = 'l')
  
  ### Scoring cell lines based on the searching space assembled in the previous examples
  CSscores <- CELLector.Score(NavTab = CSS$navTable, CELLlineData = mat_to_project)
  #knitr::kable(CSscores,align = 'l')
  
  CSscores_tab <- CSscores[!duplicated(CSscores$CellLines),]
  
  search_output <- list(signatures = Signatures_tab, projection = CSscores_tab, plot = space_plot)
  
  return(search_output)
  
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

##################################
#### plot frequency mutations ####
##################################

plot_frequency(TGCA_mat,  count_thr = 5, title = 'tgca', 
               save_file = 'output/plot/tgca_frequency_mutations.png')
plot_frequency(icgc_mat,  count_thr = 1, title = 'icgc', 
               save_file = 'output/plot/icgc_frequency_mutations.png')
plot_frequency(bailey_mat,  count_thr = 10, title = 'bailey', 
               save_file = 'output/plot/bailey_frequency_mutations.png')

# organoid data to project space on
plot_frequency(inhouse_mat,  count_thr = 1, lines = T, title = 'inhouse', 
               save_file = 'output/plot/inhouse_frequency_mutations.png')
plot_frequency_line(t(inhouse_mat),  count_thr = 1, lines = T, title = 'inhouse', 
                    save_file = 'output/plot/inhouse_frequency_lines.png')

### convert data to project into proper format ###
inhouse_mat <- t(inhouse_mat)
inhouse_mat <- cbind(id = 1:nrow(inhouse_mat), CellLine=rownames(inhouse_mat), 
                     as.data.frame(inhouse_mat))
inhouse_mat$CellLine <- factor(inhouse_mat$CellLine)

### combine all model tables ###
common_mut <- intersect(intersect(rownames(TGCA_mat), rownames(icgc_mat)), rownames(bailey_mat))
tot_mat <- cbind(TGCA_mat[common_mut,], icgc_mat[common_mut,], bailey_mat[common_mut,])

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
df_nmut_per_sample <- data.frame(name = colnames(tot_mat),label = colnames(tot_mat), n = colSums(tot_mat), lab = 'tot')
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

##########################################
### build search space for each object ###
##########################################
model_mat <- list(tgca = TGCA_mat, icgc = icgc_mat, bailey = bailey_mat, combined = tot_mat)

search_output <- list(tgca = NULL, icgc = NULL, bailey = NULL, combined = NULL)
search_output_tp53 <- list(tgca = NULL, icgc = NULL, bailey = NULL, combined = NULL)

for(id in 1:length(model_mat)){
  
  mat <- model_mat[[id]]
  name_mat <- names(model_mat)[id]
  
  print(paste('########', name_mat, '########'))
  
  search_output[[id]] <- search_space_and_project(mat, inhouse_mat)
  #search_output[[id]]$plot
  heatmap_proj <- line_signature_heatmap(search_output[[id]]$projection, 
                                         title_pl = name_mat)
  
  save_file = sprintf('output/plot/lines_score_projected_from_%s.png', name_mat)
  png(filename = save_file, width = 4, height = 6, units = 'in', res = 300)
  heatmap_proj
  dev.off()
                      
  # subset of samples with TP53 mutation
  mat_filt <- filter_mutation(mat, 'TP53')
  search_output_tp53[[id]] <- search_space_and_project(mat_filt, inhouse_mat)
  #search_output_tp53[[id]]$plot
  heatmap_proj <- line_signature_heatmap(search_output_tp53[[id]]$projection, 
                                         title_pl = sprintf('%s (TP53)', name_mat))
  
  save_file = sprintf('output/plot/lines_score_projected_from_%s_TP53.png', name_mat)
  png(filename = save_file, width = 4, height = 6, units = 'in', res = 300)
  heatmap_proj
  dev.off()
  
}

################################################
### investigate survival based on signatures ###
################################################
#library(cgdsr) # Cancer Genomics Data Server
#cgds_db <- CGDS("http://www.cbioportal.org/")

#name_studies <- getCancerStudies(cgds_db)[,c(1,2)]
#id_paad <- which(grepl('paad',name_studies$cancer_study_id))
#name_paad <- name_studies$cancer_study_id[id_paad]

#for(id_paad in name_paad){
#  case_list <- getCaseLists(cgds_db,id_paad)[1,1]
#  clinical_data <- getClinicalData(cgds_db,case_list)  
#  print(head(clinical_data))
#}




