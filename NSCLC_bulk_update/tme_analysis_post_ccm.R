# This script loads ccm fit
# Extract the TME labels
# Visualize CT-test results for specific TME vs other
# Visualize CT-test results for specific TME vs a different specific TME
# Visualize cell contribution results (box plot)
# NOTE: you have to be on the same branch in ccm.pipeline as the image that used to in the ccm run (in this case )

library(cytoreason.cc.client)
devtools::load_all("~/cytoreason.ccm.pipeline")
library('dplyr')
library(scales)
library(ggplot2)
out_path <- "~/"

# load latest NSCLC bulk model with TME (SO version)
ccm_wf <- "wf-c58e381360" #v41 with TME labels
ccm_wf <- "wf-26be3c3ad3" #v44 with TME on comparisons with tumor samples
ccm_fit <- as_ccm_fit(ccm_wf)

# cell type id to cell name mapping
cell_id_cell_name_map <- read.csv("~/cell_id_cell_type_map.csv")
cell_id_cell_name_map <- cell_id_cell_name_map %>%
  filter(!grepl("^CRCL_\\d+$", cell_type))
names(cell_id_cell_name_map)[names(cell_id_cell_name_map) == 'cell_type_id'] <- 'feature_id'
cell_id_cell_name_map
cell_id_cell_name_map <- cell_id_cell_name_map[!cell_id_cell_name_map$long_display_name == "CD16+ natural killer", ] #this cell appears twice



## ~ get TME labels for all samples per DS ~ ##
all_tme_labels <- data.frame()
for (dataset in names(ccm_fit$datasets)){
  curr_pdata <- pData(ccm_fit$datasets[[dataset]])
  curr_tme_data <- as.data.frame(curr_pdata[, 'data_tme_groups'])
  curr_tme_data$sample_id <- rownames(curr_pdata)
  colnames(curr_tme_data) <- c("data_tme_groups", "sample_id")
  curr_tme_data$dataset <- dataset
  all_tme_labels <- rbind(all_tme_labels, curr_tme_data)
}

all_tme_labels <- all_tme_labels[, c("dataset", "sample_id", "data_tme_groups")]
write.csv(all_tme_labels, paste0("~/TME_labels_", ccm_wf, ".csv"))


## ~ for BQ- get TME labels for all samples per DS ~ ##
all_tme_labels <- data.frame()
for (dataset in names(ccm_fit$datasets)){
  curr_data <- ccm_fit$datasets[[dataset]]
  curr_tme_data <- curr_data$tme_classifier$object$tme_predict
  all_tme_labels <- rbind(all_tme_labels, curr_tme_data)
}

check <- all_tme_labels %>%
  group_by(dataset) %>%
  summarise(unique_sample_count = n_distinct(sample_id)) %>%
  arrange(desc(unique_sample_count))

write.csv(all_tme_labels, paste0("~/TME_labels_for_BQ_", ccm_wf, ".csv"))


## ~ CT test results visualization ~ ##
comparisons <- c(
  "dom_lymph_myeloid_mild_stromal_vs_other",
  "dom_lymph_myeloid_vs_other",
  "epithelial_vs_other",
  "stromal_dom_myeloid_lymph_vs_other",
  "stromal_epithelial_vs_other"
)

levels <- c(
  "long_display_name"#,
  # "level_1",
  # "level_2"
)

not_subtypes <- c(
  "non malignant epithelial cell",
  "malignant epithelial cell",
  "stromal cell",
  "macrophage",
  "dendritic cell"
)

for (level in levels){
  for (comp in comparisons){
    comp_1 <- paste0("tumor_TME_", comp)
    comp_2 <- paste0("NSCLC:", comp)
    cttest_res_eset <- ccm_fit[["meta"]][[comp_1]][["ct_test"]][["fit"]][[comp_2]][["stats"]]
    print(dim(cttest_res_eset))
    cttest_res_merged <- merge(cttest_res_eset, cell_id_cell_name_map, by="feature_id") 
    cttest_res_merged$minus_log10_fdr_wiht_direction <- -log10(cttest_res_merged$fdr) * sign(cttest_res_merged$statistic)  
    cttest_res_merged <- cttest_res_merged[!cttest_res_merged$long_display_name != "CD16+ natural killer", ]
    cttest_res_merged <- cttest_res_merged[cttest_res_merged$long_display_name %in% not_subtypes, ]
    
    # plot
    x <- (ggplot(cttest_res_merged, aes_string(x = level, y = "minus_log10_fdr_wiht_direction", fill = "pvalue")) +
            geom_bar(stat = "identity") +
            xlab("Cell type") + ggtitle(paste0("CT-test:", comp))+
            scale_fill_gradientn(
              colors=c("red","white","blue"),
              values=rescale(c(0,0.05,1)), breaks=c(0.05, 0.25, 0.5,0.75),
              limits=c(0,1))+
            geom_hline(yintercept =1, linetype = "dashed",
                       colour = "darkred") +
            geom_hline(yintercept = -1, linetype = "dashed",
                       colour = "darkred") +
            coord_flip())
    
    # save plot
    # ggsave(x, file= paste0(out_path,"CT-test_", level,"_", comp, ".png"))
  }
}



## ~ cell contribution ~##

# TCGA
cell_contribution <- exprs(analysisResultExpressionSet(ccm_fit$datasets$TCGA, name = "cell_contribution"))
cell_contribution <- reshape2::melt(cell_contribution)
colnames(cell_contribution) <- c( "feature_id", "sample_id", "value")
cell_contribution_merged <- merge(cell_contribution, cell_id_cell_name_map, by="feature_id") 
cell_contribution_merged<- merge(cell_contribution_merged, all_tme_labels, by="sample_id")

# GSE174330
cell_contribution <- exprs(analysisResultExpressionSet(ccm_fit$datasets$GSE174330, name = "cell_contribution"))
cell_contribution <- reshape2::melt(cell_contribution)
colnames(cell_contribution) <- c( "feature_id", "sample_id", "value")
curr_pdata <- pData(ccm_fit$datasets$GSE174330) # take only tumor samples
tumor_samples <- rownames(curr_pdata[curr_pdata$group__TvsTadj=="B", ])
cell_contribution <- cell_contribution[cell_contribution$sample_id %in% tumor_samples, ]
cell_contribution_merged <- merge(cell_contribution, cell_id_cell_name_map, by="feature_id") 
cell_contribution_merged<- merge(cell_contribution_merged, all_tme_labels, by="sample_id")


tme_groups <- comparisons <- c(
  "dom_lymph_myeloid_mild_stromal",
  "dom_lymph_myeloid",
  "epithelial",
  "stromal_dom_myeloid_lymph",
  "stromal_epithelial"
)

cells <- unique(cell_contribution_merged$long_display_name)
  
for (cell in cells){
  curr_data <- cell_contribution_merged[cell_contribution_merged$long_display_name==cell, ]
  
  stats <- curr_data %>%
    group_by(data_tme_groups) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      std = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create a text label combining the statistics
  stats <- stats %>%
    mutate(label = sprintf("Median: %.2e\nMean: %.2e\nSTD: %.2e", median, mean, std))
  
  
  plot <- ggplot(curr_data, aes(x = data_tme_groups, y = value)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(width = 0.2, color = "blue", alpha = 0.2) +  # Plot each data point
    labs(x = "", y = "Cell contribution", title = cell) +
    theme_minimal() +
    geom_text(
      data = stats,
      aes(x = data_tme_groups, y = max(curr_data$value, na.rm = TRUE), label = label),
      inherit.aes = FALSE,
      vjust = .5,  # Position the text slightly above the boxplot
      size = 3       # Adjust text size as needed
    )
  ggsave(plot, file= paste0(out_path,"GSE174330_cell_contribution_", cell, ".png"), width = 10, bg='white')
  
}


## ~ counting samples per TME group - Tumor NSCLC ~##
tumor_nsclc_samples <- read.csv('~/tumor_NSCLC_samples.csv')
all_tme_labels_nsclc <- all_tme_labels[all_tme_labels$sample_id %in% tumor_nsclc_samples$sample_id & all_tme_labels$value == 1, ]
as.data.frame(table(all_tme_labels_nsclc$predict))




## ~ CT test - TME group1 vs TME group2 (in ccm fit we only get TME group1 vs other) ~##
run_ct_test_estimation <- function(k_to_explore, clustering_data, group_one, title=NaN){ # this function will work when there are only 2 clusters
  # define number of clusters
  sample_to_cluster <- as.data.frame(k_to_explore)
  #check
  colnames(clustering_data) %in% rownames(sample_to_cluster)
  ncol(clustering_data) == nrow(sample_to_cluster)
  # combine
  clustering_data <- rbind(clustering_data, t(sample_to_cluster))
  end <- as.numeric(dim(clustering_data)[1]) # number of rows: cells + cluster row
  rownames(clustering_data)[end] <- "kmeans_cluster"
  group <- clustering_data[end,]
  data <- clustering_data[1:(end-1),]
  
  group[group != group_one] <- "all"
  group[group == group_one] <- "one"
  print(paste0('results are for ', group_one, ' vs the other group'))
  res <- cytoreason.deconvolution::service_ct_test(data, group = group)
  
  plot_df <- res$pvalues %>%
    dplyr::mutate(feature_id = reorder(factor(rownames(res$pvalues)), estimate),
                  Direction = ifelse(estimate > 0, "Up", "Down"))
  
  # plot original
  print(ggplot(plot_df, aes(x = feature_id, y = estimate, fill = Direction)) +
          geom_bar(stat = "identity") +
          xlab("Cell type") + ggtitle(paste0("CT-test estimation result for cluster ", j," in K = ",max(unique_groups)))+
          scale_fill_brewer(palette = "Set1", direction = -1, aes(colour = Direction))+
          coord_flip())

  group <- clustering_data[end,]
  
  # plot I added
  plot_df$minus_log10_fdr_with_direction <- -log10(plot_df$FDR) * sign(plot_df$estimate)  
  plot_df <- plot_df[order(plot_df$minus_log10_fdr_with_direction), ]
  plot_df$feature_id <- factor(plot_df$feature_id, levels = plot_df$feature_id)
  print(ggplot(plot_df, aes(x = feature_id, y = minus_log10_fdr_with_direction, fill = p.value.adj)) +
          geom_bar(stat = "identity") +
          xlab("Cell type") + ggtitle(paste0("CT-test\n", title))+
          scale_fill_gradientn(
            colors=c("red","white","blue"),
            values=rescale(c(0,0.05,1)), breaks=c(0.05, 0.25, 0.5,0.75),
            limits=c(0,1))+
          geom_hline(yintercept =1, linetype = "dashed",
          colour = "darkred") +
          geom_hline(yintercept = -1, linetype = "dashed",
          colour = "darkred") +
          coord_flip())
  group <- clustering_data[end,]
 return(plot_df)
}


#### cell contribution per sample
all_cell_cont <- data.frame()
for (dataset in names(ccm_fit$datasets)){
  curr_data <- ccm_fit$datasets[[dataset]]
  curr_cell_cont_data <- as.data.frame(curr_data$cell_contribution$eset)
  print(dim(curr_cell_cont_data))
  cell_ids <- grep("^CRCL_[0-9]{7}$", colnames(curr_cell_cont_data), value = TRUE)
  curr_cell_cont_mat <- curr_cell_cont_data[, cell_ids]
  all_cell_cont <- rbind(all_cell_cont, curr_cell_cont_mat)
}
all_cell_cont <- t(all_cell_cont)
dim(all_cell_cont)
# take only tumor samples (from tumor_vs_tumor_adjacent comparison, label = NSCLC)
tumor_nsclc_samples <- read.csv('~/tumor_NSCLC_samples.csv')
all_cell_cont <- all_cell_cont[, colnames(all_cell_cont) %in% tumor_nsclc_samples$sample_id]
dim(all_cell_cont)
# convert cell id to cell name
rownames(all_cell_cont) <- cell_id_cell_name_map$long_display_name[match(rownames(all_cell_cont), cell_id_cell_name_map$feature_id)]
# remove parents
not_subtypes <- c(
  "non malignant epithelial cell",
  "malignant epithelial cell",
  "stromal cell",
  "macrophage",
  "dendritic cell"
)
all_cell_cont <- all_cell_cont[!rownames(all_cell_cont) %in% not_subtypes, ]


#### TME labels per sample
TME_labels <- all_tme_labels[all_tme_labels$value==1,]
TME_labels_input <- as.data.frame(TME_labels$predict)
rownames(TME_labels_input) <- TME_labels$sample_id
head(TME_labels_input)
# take only tumor samples (from tumor_vs_tumor_adjacent comparison, label = NSCLC)
TME_labels_input <- TME_labels_input[rownames(TME_labels_input) %in% tumor_nsclc_samples$sample_id, , drop = FALSE]
head(TME_labels_input)

# ordering the samples in the same or
TME_labels_input <- TME_labels_input[order(match(rownames(TME_labels_input), colnames(all_cell_cont))), , drop = FALSE]
colnames(TME_labels_input) <- c("TME")
dim(TME_labels_input)

TMEs = c(
  'stromal_dom_myeloid_lymph', #1
  'dom_lymph_myeloid', #2 
  'epithelial', #3
  'dom_lymph_myeloid_mild_stromal',#4
  'stromal_epithelial' #5
)

map_TME_to_num <- data.frame(
  TME = TMEs,
  number = 1:length(TMEs)
)
TME_labels_input_merged <- merge(TME_labels_input, map_TME_to_num, by = "TME", all.x = TRUE)
rownames(TME_labels_input_merged) <- rownames(TME_labels_input)

k_to_explore <- as.data.frame(TME_labels_input_merged['number'])
clustering_data <- all_cell_cont

# chosing the TME groups to compare
group1_cluster <- 2 #'dom_lymph_myeloid'
group2_cluster <- 5 #'stromal_epithelial'

k_to_explore <- k_to_explore[k_to_explore[,1] %in% c(group1_cluster, group2_cluster), ,  drop = FALSE]
clustering_data <- clustering_data[, colnames(clustering_data) %in% rownames(k_to_explore)]
title <- paste0(map_TME_to_num$TME[map_TME_to_num$number==group1_cluster], ' VS ', map_TME_to_num$TME[map_TME_to_num$number==group2_cluster])
print(title)

res <- run_ct_test_estimation(k_to_explore, clustering_data, group_one=group1_cluster, title=title)

write.csv(res, '~/plot_df.csv')
