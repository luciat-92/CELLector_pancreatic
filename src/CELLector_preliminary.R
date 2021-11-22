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

summary_CSS <- function(Signatures, ModelMat, n_patients){
  
  # Signature: output of CELLector.createAllSignatures
  # ModelMat: output of CELLector.buildModelMatrix
  # n_patients : number of patients for search space
  
  df <- data.frame(Subtype = names(Signatures$ES), 
                   Signatures = Signatures$ES, 
                   Signatures_complete = Signatures$S,
                   N.patients = round(n_patients*Signatures$STS/100), 
                   P.patients = Signatures$STS, 
                   repr_CL = NA)
  
  id_mapped <- rownames(ModelMat)
  id_lines <- colnames(ModelMat)
  
  for(id_row in 1:nrow(df)){
    
    if(!df$Subtype[id_row] %in% id_mapped){
      df$repr_CL[id_row] <- 'lack of in vitro models'
    }else{
      id_row_model <- which(id_mapped == df$Subtype[id_row])
      df$repr_CL[id_row] <- paste0(id_lines[ModelMat[id_row_model, ] == 1], collapse = ',')  
    }
  }
  return(df)
}

coverage_CSS <- function(navTable){
  
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

coverage_CSS_mapped_lines <- function(navTable, subset_proj){
  
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
  #selectedCellLines <- CELLector.makeSelection(modelMat = ModelMat,
  #                                             n=ncol(ModelMat), 
  #                                             searchSpace = CSS$navTable)
  #knitr::kable(selectedCellLines,align = 'l')
  
  ### Scoring cell lines based on the searching space assembled in the previous examples
  CSscores <- CELLector.Score(NavTab = CSS$navTable, CELLlineData = mat_to_project)
  #knitr::kable(CSscores,align = 'l')
  
  CSscores_tab <- CSscores[!duplicated(CSscores$CellLines),]
  # get summary space:
  summary_res <- summary_CSS(Signatures = Signatures, 
                             ModelMat = ModelMat, 
                             n_patients = ncol(mat))
  
  # get percentage coverage
  subset_proj <- summary_res[summary_res$Signatures_complete %in% unique(CSscores_tab$Signature),]
  n_samples_mapped <- coverage_CSS(navTable = CSS$navTable)
  n_samples_mapped_with_lines <- coverage_CSS_mapped_lines(navTable = CSS$navTable, 
                                                           subset_proj = subset_proj)
  
  search_output <- list(signatures = Signatures_tab, projection = CSscores_tab, 
                        summary = summary_res, 
                        n_samples_mapped = n_samples_mapped, 
                        n_samples_mapped_with_lines = n_samples_mapped_with_lines, 
                        plot = space_plot)
  
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
res = inhouse_mat
save(res, file = 'data/inhouse_snv_matrix_drivers.RData')

### combine all model tables ###
common_mut <- intersect(intersect(rownames(TGCA_mat), rownames(icgc_mat)), rownames(bailey_mat))
tot_mat <- cbind(TGCA_mat[common_mut,], icgc_mat[common_mut,], bailey_mat[common_mut,])
# save .RData format
res <- tot_mat
save(res, file = 'data/combined_snv_matrix_drivers.RData')

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

# save summary table for combined 
summary_tab <- search_output[[4]]$summary
write.table(summary_tab, file = 'output/combined_summary_CSS.tsv', quote = F, 
            sep = '\t', row.names = F, col.names = T)

# pie chart mapping
n_samples <- ncol(tot_mat)
missing <- n_samples - search_output[[4]]$n_samples_mapped
missing_in_lines <- search_output[[4]]$n_samples_mapped - sum(search_output[[4]]$n_samples_mapped_with_lines$n_samples)

df_coverage <- data.frame(n = c(search_output[[4]]$n_samples_mapped_with_lines$n_samples, 
                                missing_in_lines, missing), 
                          type = c(search_output[[4]]$n_samples_mapped_with_lines$Signatures, 
                                   'without organoid model', 
                                   'not in CELLector search space'))
df_coverage <- df_coverage %>% 
               mutate(prop = n/n_samples*100, 
                      prop_lab =  paste0(round(n/n_samples*100, digits = 1), '%'))
df_coverage$type <- factor(df_coverage$type, levels = rev(df_coverage$type))
mycols <- c("grey90", "grey50", "lightsalmon4", "lightsalmon3", "lightsalmon1")

pie_chart <- ggplot(df_coverage, aes(x = "", y = prop, fill = type)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0)+
    geom_label_repel(aes(label = prop_lab), position = position_stack(vjust=0.5), min.segment.length = 0) +
    scale_fill_manual(values = mycols) +
    theme_void()

ggsave(pie_chart, filename = 'output/plot/combined_piechart_with_signatures.png',
       width = 4, height = 4, dpi = 200)

###
df_coverage <- data.frame(n = c(search_output[[4]]$n_samples_mapped, missing),
                          type = c( 'in CELLector search space', 
                                    'not in CELLector search space'))
df_coverage <- df_coverage %>% 
  mutate(prop = n/n_samples*100, 
         prop_lab =  paste0(round(n/n_samples*100, digits = 1), '%'))
df_coverage$type <- factor(df_coverage$type, levels = rev(df_coverage$type))
mycols <- c("grey90", "lightblue")

pie_chart <- ggplot(df_coverage, aes(x = "", y = prop, fill = type)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  geom_label_repel(aes(label = prop_lab), position = position_stack(vjust=0.5), min.segment.length = 0) +
  scale_fill_manual(values = mycols) +
  theme_void()

ggsave(pie_chart, filename = 'output/plot/combined_piechart.png',
       width = 4, height = 4, dpi = 200)



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




