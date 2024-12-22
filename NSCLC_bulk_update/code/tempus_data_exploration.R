###
# used this code to generate TEMPUS eset for NSCLC bulk model update
# used this section: "Generate full eset with Final annotations from data management - Excluding During, gap and Unknown"
# this code also read and save as csv the pdata of Pfizer datasets and TCGA
# generates boxplots of cell contribution of LUAD/LUSC (data queried from BQ)

library(ggplot2)
library(reshape2)

library(cytoreason.cc.client)
library(cytoreason.ccm.pipeline)
# library(cytoreason.shared.assets)
# library(cytoreason.curator.annotations)
library(cytoreason.datasource)



### Functions
cytocc_save <- function(obj) {
  return(obj)
}

convert_gene_symbol_to_entrez_id <- function(expr_df){
  temp <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(expr_df), "ENTREZID", "SYMBOL")
  n_nan <- sum(is.na(temp)) #1,336 Nas out of 19,377 (7%)
  print(paste0(n_nan, " gene symbols were not converted (", round(100*n_nan/dim(expr_df)[1]), "%)"))
  temp <- temp[!is.na(temp)]
  expr_df_sub <- expr_df[names(temp), ]
  rownames(expr_df_sub) <- temp
  return(expr_df_sub)
}


### load the data from Pfizer
setwd('~/p01-mnsclc-tempus-data/')
load('TEMPUS_Provided.log2TPM.corrected.RData')
load('expr_count_corrected_HGNC.RData')
patient_info <- read.csv('master_patient_info.v20240328.csv')
names(patient_info)[names(patient_info) == 'PATIENT_ID'] <- 'patient_id'
sample_info <- read.csv('master_sample_info.v20240328.csv')

expr_log2TPM <- expr_log2tpm_corrected_HGNC
expr_count <- expr_count_corrected_HGNC

rm(expr_log2tpm_corrected_HGNC, expr_count_corrected_HGNC)


### Explore metadata - #samples, #patients, #relevant samples (lung)
dim(patient_info)
length(unique(patient_info$patient_id))

dim(sample_info)
length(unique(sample_info$sample_id))

# samples in info but no data
length(sample_info$sample_id[!sample_info$sample_id %in% colnames(expr_count)])
# samples in data but no info
length(colnames(expr_count)[!colnames(expr_count) %in% sample_info$sample_id])
# common samples with data and info
samples_with_data_and_info <- intersect(sample_info$sample_id, colnames(expr_count))
patients_with_samples_with_data_and_info <- unique(sample_info$patient_id[sample_info$sample_id %in% samples_with_data_and_info])
length(samples_with_data_and_info)
length(patients_with_samples_with_data_and_info)

# take only samples with data and info
sample_info <- sample_info[sample_info$sample_id %in% samples_with_data_and_info, ]
patient_info <- patient_info[patient_info$patient_id %in% patients_with_samples_with_data_and_info, ]
expr_log2TPM <- expr_log2TPM[, samples_with_data_and_info]
expr_count <- expr_count[, samples_with_data_and_info]


patients_counts <- as.data.frame(table(sample_info$patient_id))
tissue_counts <- as.data.frame(table(sample_info$tissue_site_broad))

lung_samples <- sample_info[sample_info$tissue_site_broad == 'Lung',]
dim(lung_samples) # 291 samples
length(unique(lung_samples$patient_id)) # 283 patients
patients_counts_lung <- as.data.frame(table(lung_samples$patient_id))

write.csv(lung_samples, 'lung_samples.csv')


### Subset lung relevant samples from raw counts and log2TPM
# check that all samples have data
samples_no_data <- lung_samples$sample_id[!lung_samples$sample_id %in% colnames(expr_log2TPM)] # 6 samples have no data available
samples_with_data <- lung_samples$sample_id[lung_samples$sample_id %in% colnames(expr_log2TPM)] # 291 samples have data
expr_log2TPM_lung <- expr_log2TPM[, samples_with_data]
expr_count_lung <- expr_count[, samples_with_data]
sample_info_lung <- sample_info[sample_info$sample_id %in% samples_with_data, ]


### merge patients info into sample info table
pdata <- merge(sample_info_lung, patient_info, by="patient_id") 
rownames(pdata) <- pdata$sample_id

### Histogram of log2TPM values combined
df_melted <- melt(expr_log2TPM_lung)

combined_histogram <- ggplot(df_melted, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "green", color = "black") +
  theme_minimal() 
print(combined_histogram)


### convert gene symbols to entrez ids
expr_count_lung_sub <- convert_gene_symbol_to_entrez_id(expr_count_lung)
expr_log2TPM_lung_sub <- convert_gene_symbol_to_entrez_id(expr_log2TPM_lung)

# temp <- expr_count_lung_sub %>% rownames_to_column(var = "sample_id")
# write.table(temp, file.path("data/bulk_RNAseq/Tempus", "Tempus_Count_sub_entrezid.csv"), quote = FALSE, sep = ",", row.names = FALSE)


### Generate full eset: convert to entrez ids + upload to cyto-cc the full data (all tissues)
expr_count_sub <- convert_gene_symbol_to_entrez_id(expr_count)
expr_log2TPM_sub <- convert_gene_symbol_to_entrez_id(expr_log2TPM)

pdata_full <- merge(sample_info, patient_info, by="patient_id") 
rownames(pdata_full) <- pdata_full$sample_id

# non_na_columns <- pdata_full[, colSums(is.na(pdata_full)) != nrow(pdata_full)]

write.csv(pdata_full, 'tempus_pdata_original.csv')

#*****# 1. Add "tissue_site_broad_HISTOLOGY_TYPE2_AT_BIOPSY" column
col_to_add <- as.data.frame(paste(pdata_full$tissue_site_broad, pdata_full$HISTOLOGY_TYPE2_AT_BIOPSY, sep = "_"))
colnames(col_to_add) <- "tissue_site_broad_HISTOLOGY_TYPE2_AT_BIOPSY"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 2. Add "primary_metastasis_added" column
col_to_add <- as.data.frame(pdata_full$tissue_site_broad)
colnames(col_to_add) <- "primary_metastasis_added"
col_to_add$primary_metastasis_added[col_to_add$primary_metastasis_added == "Lung"] = 'primary'
col_to_add$primary_metastasis_added[!col_to_add$primary_metastasis_added == "primary"] = 'metastasis'
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 3. Add "tissue_site_broad_PROGRESSION_STATUS" column
col_to_add <- as.data.frame(paste(pdata_full$tissue_site_broad, pdata_full$PROGRESSION_STATUS, sep = "_"))
colnames(col_to_add) <- "tissue_site_broad_PROGRESSION_STATUS"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 4. Add "tissue_site_broad_race_concept_canonical_name" column
col_to_add <- as.data.frame(paste(pdata_full$tissue_site_broad, pdata_full$race_concept_canonical_name, sep = "_"))
colnames(col_to_add) <- "tissue_site_broad_race_concept_canonical_name"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 5. Add "tissue_site_broad_GENDER" column
col_to_add <- as.data.frame(paste(pdata_full$tissue_site_broad, pdata_full$GENDER, sep = "_"))
colnames(col_to_add) <- "tissue_site_broad_GENDER"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 6. Add "tissue_site_broad_SMOKING_AT_BIOPSY" column
col_to_add <- as.data.frame(paste(pdata_full$tissue_site_broad, pdata_full$SMOKING_AT_BIOPSY, sep = "_"))
colnames(col_to_add) <- "tissue_site_broad_SMOKING_AT_BIOPSY"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 7. Add "heavily_lightly_treated" column
col_to_add <- as.data.frame(pdata_full$ICI_LOT_COUNT)
colnames(col_to_add) <- "heavily_lightly_treated"
col_to_add[col_to_add==4] = NaN
col_to_add[col_to_add==1] = "lightly"
col_to_add[col_to_add==2] = "heavily"
col_to_add[col_to_add==3] = "heavily"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 8. Add "responder_non_responder_added" column
col_to_add <- as.data.frame(pdata_full$FIRST_RESPONSE)
colnames(col_to_add) <- "responder_non_responder_added"
col_to_add[col_to_add=="Partial Response"] = "responder"
col_to_add[col_to_add=="Progressive Disease"] = "non_responder"
col_to_add[col_to_add=="Stable Disease"] = "non_responder"
col_to_add[col_to_add=="Complete Response"] = "responder"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 9. Add "heavily_lightly_treated_responder_non_responder_added" column
col_to_add <- as.data.frame(paste(pdata_full$heavily_lightly_treated, pdata_full$responder_non_responder_added, sep = "_"))
colnames(col_to_add) <- "heavily_lightly_treated_responder_non_responder_added"
pdata_full <- cbind(pdata_full, col_to_add)


#*****# 10. Add "tissue_site_broad_heavily_lightly_treated_responder_non_responder_added" column
col_to_add <- as.data.frame(paste(pdata_full$tissue_site_broad, pdata_full$heavily_lightly_treated_responder_non_responder_added, sep = "_"))
colnames(col_to_add) <- "tissue_site_broad_heavily_lightly_treated_responder_non_responder_added"
pdata_full <- cbind(pdata_full, col_to_add)


write.csv(pdata_full, 'tempus_pdata_with_additional_columns.csv')


### Build the full eset
# pdata_anno_full <- new("AnnotatedDataFrame", data = pdata_full)
# eset_full <- Biobase::ExpressionSet(
#   assayData = as.matrix(expr_log2TPM_sub),
#   phenoData = pdata_anno_full,
#   annotation = "org.Hs.eg.db"
# )

ann_df <- annotateFeatures(expr_log2TPM_sub, 
                           annotation = "org.Hs.eg.db", 
                           columns = c("ENTREZID", "SYMBOL", "GENENAME", "CHR"))
eset_full <- ExpressionSet(assayData = as.matrix(expr_log2TPM_sub), 
                     phenoData = AnnotatedDataFrame(pdata_full),
                     featureData = AnnotatedDataFrame(ann_df),
                     annotation = "org.Hs.eg.db")
assayDataElement(eset_full, "counts") <- as.matrix(expr_count_sub)


### upload eset to cyto-cc
wf <- run_function_dist(cytocc_save, obj = eset_full)
# Cyto-CC workflow: wf-216de888a4 - with original pdata
# Cyto-CC workflow: wf-d2af40fa2b - with 1-4 columns added to pdata
# Cyto-CC workflow: wf-bd3ccb40ea - correcting featureData to include mapping of entrez id, gene symbol, gene name and CHR
# Cyto-CC workflow: wf-2e0defd1aa - checking if ccm can run with symbols instead of entrez ids (it fails)
# Cyto-CC workflow: wf-cec82ef68e - adding columns for comparisons of gender, smoking and heavily/lightly_treated_response




### Generate Lung eset: convert to entrez ids + upload to cyto-cc the primary tumor (lung) data
lung_samples <- pdata_full$sample_id[pdata_full$tissue_site_broad=='Lung']
print(length(lung_samples))
expr_log2TPM_sub_lung <- expr_log2TPM_sub[,lung_samples]
print(dim(expr_log2TPM_sub_lung))
ann_df_lung <- annotateFeatures(expr_log2TPM_sub_lung, 
                           annotation = "org.Hs.eg.db", 
                           columns = c("ENTREZID", "SYMBOL", "GENENAME", "CHR"))
pdata_lung <- pdata_full[pdata_full$sample_id %in% lung_samples,]
print(dim(pdata_lung))
expr_count_sub_lung <- expr_count_sub[, lung_samples]
print(dim(expr_count_sub_lung))
eset_lung <- ExpressionSet(assayData = as.matrix(expr_log2TPM_sub_lung), 
                           phenoData = AnnotatedDataFrame(pdata_lung),
                           featureData = AnnotatedDataFrame(ann_df_lung),
                           annotation = "org.Hs.eg.db")
assayDataElement(eset_lung, "counts") <- as.matrix(expr_count_sub_lung)

### upload eset to cyto-cc
wf <- run_function_dist(cytocc_save, obj = eset_lung)
# Cyto-CC workflow: wf-e01a273d27 - only lung samples 

asset_link <- paste0("ccw://", "wf-e01a273d27", ":0:output.rds")
curr_eset <- read_asset(asset_link)


asset_link <- paste0("ccw://", "wf-bd3ccb40ea", ":0:output.rds")
original_eset <- read_asset(asset_link)

asset_link <- paste0("ccw://", "wf-ee55ca0222", ":0:output.rds")
curr_eset <- read_asset(asset_link)
saveRDS(curr_eset, file = "TEMPUS_wf-ee55ca0222_eset.rds")



### Generate full eset with Final annotations from data management - Excluding During, gap and Unknown
expr_count_sub <- convert_gene_symbol_to_entrez_id(expr_count)
expr_log2TPM_sub <- convert_gene_symbol_to_entrez_id(expr_log2TPM)

final_annotations <- read.csv('Final annotations.csv')
keep_samples <- final_annotations$sample_id[final_annotations$pre_post %in% c("Pre","Post")]
length(keep_samples)
keep_samples_common <- keep_samples[keep_samples %in% colnames(expr_log2TPM_sub)]
expr_count_sub_filtered <- expr_count_sub[,keep_samples_common]
expr_log2TPM_sub_filtered <- expr_log2TPM_sub[,keep_samples_common]
dim(expr_count_sub_filtered)
dim(expr_log2TPM_sub_filtered)


# final_annotations <- final_annotations[samples_with_data_and_info,]
rownames(final_annotations) <- final_annotations$sample_id
final_annotations <- final_annotations[keep_samples_common,]
dim(final_annotations)


# setdiff(rownames(final_annotations), colnames(expr_log2TPM_sub))
# setdiff(colnames(expr_log2TPM_sub), rownames(final_annotations))

ann_df <- annotateFeatures(expr_log2TPM_sub_filtered, 
                           annotation = "org.Hs.eg.db", 
                           columns = c("ENTREZID", "SYMBOL", "GENENAME", "CHR"))
eset_full <- ExpressionSet(assayData = as.matrix(expr_log2TPM_sub_filtered), 
                           phenoData = AnnotatedDataFrame(final_annotations),
                           featureData = AnnotatedDataFrame(ann_df),
                           annotation = "org.Hs.eg.db")
assayDataElement(eset_full, "counts") <- as.matrix(expr_count_sub_filtered)


### upload eset to cyto-cc
wf <- run_function_dist(cytocc_save, obj = eset_full)

# Cyto-CC workflow: wf-c0239bd55f - 973 samples with final annotations (without counts because I had an error)
# Cyto-CC workflow: wf-921f391350 - 739 samples (without Gap/During/Unknown)
# Cyto-CC workflow: wf-ffe867c4d8 - 739 samples (without Gap/During/Unknown) - no spaces in column names
# Cyto-CC workflow: wf-e6baada5b9 - after the corrections weiging email- failed image error
# Cyto-CC workflow: wf-448893973d - trying again - failed because of tissue_comment_tmb_pd-l1_time column (changed to tissue_comment_tmb_pd_l1_time)
# Cyto-CC workflow: wf-987b9033da - trying again
# Cyto-CC workflow: wf-b9dcffe410 - after Adi fixed the column for TMB comparison 
# Cyto-CC workflow: wf-a01cffdf0d - changed in all column names from pd-l1 to pd_l1
# Cyto-CC workflow: wf-8606a5b125 - adding condition_subtype_race and tissue_comment_time_gender for the general race and gender comparisons
# Cyto-CC workflow: wf-65d4b339ac - adding columns for general comparisons that include metastasis column
# Cyto-CC workflow: wf-05381f03ed- adding columns for post_vs_pre paired analysis
# Cyto-CC workflow: wf-ee55ca0222 - making sure all columns have "_" and not spaces


eset <- read_asset("ccw://wf-ee55ca0222:0:output.rds")
pdata <- pData(eset)
unique(pdata$pre_post_year_organ_group_liver_included_paired)



eset_b@featureData
eset_full@featureData <- eset_b@featureData


check_empty_characters <- function(data) {
  # Apply a function to each element to check if it is an empty string
  empty_matrix <- apply(data, c(1, 2), function(x) x == "")
  return(empty_matrix)
}
a <- check_empty_characters(pdata_full)
# check the uploaded eset
eset_curr <- read_asset("ccw://wf-d2af40fa2b:0:output.rds")

eset_b <- read_asset("ccw://wf-84cd32937b:0:output.rds")


### generate eset
rownames(pdata) <- pdata$sample_id
pdata_anno <- new("AnnotatedDataFrame", data = pdata)
eset <- Biobase::ExpressionSet(
  assayData = as.matrix(expr_log2TPM_lung_sub),
  phenoData = pdata_anno
)
assayDataElement(eset, "counts") <- as.matrix(expr_count_lung_sub)

### upload eset to cyto-cc
wf <- run_function_dist(cytocc_save, obj = eset) 

# Cyto-CC workflow: wf-262e7830fc





### Read all other Pfizer eset and save their metadata
ds_wf_lst <- list("B9991004" = "wf-84cd32937b", 
                  "B8011001" = "wf-ac5285867f",
                  "B8011007" = "wf-233ef5d5dd",
                  "B9991027" = "wf-551747dbd3")

ds_wf_df <- data.frame(
  ds = names(ds_wf_lst),  # Extract the keys as the "ds" column
  wf = unlist(ds_wf_lst)  # Extract the values as the "wf" column
)

for (ds in ds_wf_df$ds){
  asset_link <- paste0("ccw://", ds_wf_df$wf[ds_wf_df$ds==ds], ":0:output.rds")
  pdata_filename <- paste0("pData_", ds, "_", ds_wf_df$wf[ds_wf_df$ds==ds], ".csv")
  curr_eset <- read_asset(asset_link)
  curr_p_data <- pData(curr_eset)
  write.csv(curr_p_data, pdata_filename)
}


### Read TCGA eset and it's metadata
asset_link <- paste0("ccw://", "wf-be092d2791", ":0:output.rds")
pdata_filename <- paste0("pData_", "TCGA", "_", "wf-be092d2791", ".csv")
curr_eset <- read_asset(asset_link)
curr_p_data <- pData(curr_eset)
write.csv(curr_p_data, pdata_filename)


#####
asset_link <- paste0("ccw://wf-551747dbd3:0:output.rds")
curr_eset <- read_asset(asset_link)
curr_expr <- exprs(curr_eset)
df_melted <- melt(curr_expr)

combined_histogram <- ggplot(df_melted, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "green", color = "black") 
print(combined_histogram)




# Function to generate histograms
generate_histograms <- function(df) {
  # Melt the data frame to long format for combined histogram
  df_melted <- melt(df)
  # # 1. Histograms for each column separately in subplots
  # separate_histograms <- ggplot(df_melted, aes(x = value)) +
  #   geom_histogram(binwidth = 30, fill = "blue", color = "black") +
  #   facet_wrap(~ variable, scales = "free") +
  #   theme_minimal() +
  #   labs(title = "Histograms for Each Column Separately")
  # 
  # print(separate_histograms)

  # 2. Histogram of all values combined
  combined_histogram <- ggplot(df_melted, aes(x = value)) +
    geom_histogram(binwidth = 0.5, fill = "green", color = "black") +
    theme_minimal() 
    # + labs(title = "Combined Histogram of All Values")
  print(combined_histogram)
}


result <- sapply(lung_samples$sample_id, function(x) x %in% colnames(expr_count_corrected_HGNC))
result[!result] # 6 samples out of 297 don't have gene expression data
sum(result)

result <- sapply(lung_samples$sample_id, function(x) x %in% colnames(expr_log2tpm_corrected_HGNC))
result[!result]

counts_lung <- expr_count_corrected_HGNC[, lung_samples$sample_id[result]]
log2tpm_lung <- expr_log2tpm_corrected_HGNC[, lung_samples$sample_id[result]]



generate_histograms(log2tpm_lung)




       
####### LUAD/LUSC cell contribution
data <- read.csv('~/bq-results-20241028-103649-1730111864169.csv')

# used this query from BQ to generate the data:
# With x as (
#   
#   WITH cell_data AS (
#     SELECT cells.*, cid2l.level_1_name AS level_1, cid2l.level_2_name AS level_2 FROM 
#     (SELECT cell_type_id AS cell_type,
#       cell_type_id as cell_type_id,
#       long_display_name
#       FROM `cytoreason.cr_annotations_prod.cell_types_versioned` where dbt_valid_to is null
#       union all
#       SELECT distinct 
#       long_display_name AS cell_type, 
#       cell_type_id as cell_type_id,
#       long_display_name
#       FROM `cytoreason.cr_annotations_prod.cell_types_versioned`) cells
#     JOIN
#     `cytoreason.cr_annotations_prod.current_ids_to_levels` cid2l ON cells.cell_type_id = cid2l.cell_type_id
#   )
#   SELECT 
#   cc.feature_id, cc.sample_id, cc.value, cc.version, cc.dataset_id, mgcs.group_label, cell_data.long_display_name, cell_data.level_1, cell_data.level_2
#   FROM 
#   (SELECT 
#     *, 
#     SUBSTR(sample_id, INSTR(sample_id, ' ') + 1) AS sample_id_new 
#     FROM 
#     `cytoreason.p01_io_nsclc.cell_contribution`
#     WHERE version=37
#   ) as cc 
#   JOIN 
#   `p01_io_nsclc.model_group_comparison_sample` as mgcs
#   ON 
#   cc.sample_id_new = mgcs.sample_id
#   JOIN 
#   cell_data 
#   ON 
#   cc.feature_id = cell_data.cell_type
#   where effect_id in ('LUAD_LUSC', 'soft_tissue_vs_primary', 'pleura_vs_primary', 'metastasis_vs_primary', 'lymph_node_vs_primary', 'liver_vs_primary', 'chest_vs_primary', 'brain_vs_primary', 'bone_vs_primary')
#   and mgcs.version=37
# )
# 
# select * from x 
# where long_display_name in ('malignant squamous epithelial cell of lung', 'malignant glandular epithelial cell of lung', 'type II alveolar cell', 'airway basal cell')
# and dataset_id = 'TEMPUS'

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)


########## all datasets with lines
ggplot(data[data$group_label=='LUAD',],aes(x = long_display_name, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers to avoid interference with lines
  geom_point(aes(color = as.factor(sample_id)), position = position_jitter(width = 0.2), alpha=0.5) +  # Add points with jitter
  geom_line(aes(group = sample_id, color = as.factor(sample_id)), position = position_jitter(width = 0.2), linewidth = 0.3) +
  labs(x = "cell name", y = "cell contribution", title = "Boxplot of B by A, Faceted by Dataset ID") +
  ggtitle("LUAD") +
  facet_wrap(~ dataset_id) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") 

ggplot(data[data$group_label=='LUSC',],aes(x = long_display_name, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers to avoid interference with lines
  geom_point(aes(color = 'b'), position = position_jitter(width = 0.2), alpha=0.5) +  # Add points with jitter
  # geom_line(aes(group = sample_id), position = position_jitter(width = 0.2), linewidth = 0.3) +
  labs(x = "cell name", y = "cell contribution", title = "Boxplot of B by A, Faceted by Dataset ID") +
  ggtitle("LUAD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") 


###### LUAD/LUSC cell abundance
# Step 1: Reshape the data to calculate differences for each sample
for (subtype in unique(data$group_label)){
  curr_data = data[data$group_label==subtype,]
  # curr_data = data[data$group_label==subtype & data$dataset_id=='TEMPUS',]
  
  # Reshape data for analysis - pivot
  data_wide <- curr_data %>%
    pivot_wider(
      id_cols = sample_id,  # Use sample_id as the unique identifier
      names_from = long_display_name,    # Take group values as column names
      values_from = value    # Use the value column for the data
    )
  
  # data_wide <- data %>%
  #   pivot_wider(names_from = long_display_name, values_from = value) %>%
  #   mutate(diff = `malignant squamous epithelial cell of lung` - `malignant glandular epithelial cell of lung`)  # Replace 'Group1' and 'Group2' with actual group names
  
  
  # plot distributions of cell abundance by cell type (check for normality)
  ggplot(data_wide, aes(x = `malignant squamous epithelial cell of lung`)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    labs(title = "Distribution of malignant squamous epithelial cell of lung", x = "cell contribution", y = "Density") +
    theme_minimal()
  
  ggplot(data_wide, aes(x = `malignant glandular epithelial cell of lung`)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    labs(title = "Distribution of malignant glandular epithelial cell of lung", x = "cell contribution", y = "Density") +
    theme_minimal()
  
  
  # Step 2: Check for normality in the differences
  shapiro_test_1 <- shapiro.test(data_wide$'malignant squamous epithelial cell of lung')
  shapiro_test_2 <- shapiro.test(data_wide$'malignant glandular epithelial cell of lung')
  
  
  # Step 3: Conduct the appropriate test based on normality
  if (min(shapiro_test_1$p.value, shapiro_test_2$p.value) > 0.05) {
    # Paired t-test if normally distributed
    test_result <- t.test(
      data_wide$`malignant squamous epithelial cell of lung`, 
      data_wide$`malignant glandular epithelial cell of lung`, 
      paired = TRUE)
    test_name <- "Paired t-test"
  } else {
    # Wilcoxon signed-rank test if not normally distributed
    test_result <- wilcox.test(data_wide$`malignant squamous epithelial cell of lung`, 
                               data_wide$`malignant glandular epithelial cell of lung`, 
                               paired = TRUE)
    test_name <- "Wilcoxon signed-rank test"
  }
  
  
  # Step 6: Create the boxplot with all dots and add text annotations
  curr_data <- curr_data %>%
    mutate(long_display_name = gsub("epithelial", "\nepithelial", long_display_name))
  
  curr_fig <- ggplot(curr_data, aes(x = long_display_name, y = value)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(width = 0.2, color = "blue", alpha = 0.2) +  # Plot each data point
    labs(x = "cell type", y = "cell contribution", title = subtype) +
    theme_minimal() +
    # Add text annotations for mean, median, standard deviation, and p-value
    annotate("text", x = 1.2, y = max(curr_data$value, na.rm = TRUE),
             label = paste0("P-value: ", sprintf("%.2e", test_result$p.value)),
             hjust = 0, vjust = 1, size = 4) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14))   # Rotate x-axis labels
    
  
  # Step 6: Violin plot
  # curr_data <- curr_data %>%
  #   mutate(long_display_name = gsub("epithelial", "\nepithelial", long_display_name))
  # 
  # medians <- curr_data %>%
  #   group_by(long_display_name) %>%
  #   summarise(median_value = median(value, na.rm = TRUE))
  # 
  # curr_fig <- ggplot(curr_data, aes(x = long_display_name, y = value)) +
  #  geom_violin(fill = "lightblue", color = "darkblue", alpha = 0.6) +  # Violin plot
  #  # geom_jitter(width = 0.2, color = "blue", alpha = 0.2) +  # Plot each data point
  #  # geom_hline(yintercept = median_diff, color = "red", linetype = "dashed", size = 1) +  # Add median line
  #  labs(x = "cell type", y = "cell contribution", title = subtype) +
  #  theme_minimal() +
  #  # Add text annotations for mean, median, standard deviation, and p-value
  #  annotate("text", x = 1.2, y = max(curr_data$value, na.rm = TRUE), 
  #           label = paste0("P-value: ", sprintf("%.2e", test_result$p.value)),
  #           hjust = 0, vjust = 1, size = 4) +
  #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  # Rotate x-axis labels
  #   geom_crossbar(data = medians, aes(x = long_display_name, y = median_value, 
  #                                     ymin = median_value, ymax = median_value),
  #                 color = "red", width = 0.3, inherit.aes = FALSE, linewidth = 0.2)
  
   ggsave(paste0("~/", subtype, "_cell_controbution_boxplot.png"), plot = curr_fig, width = 8, height = 6, bg='white',dpi = 300)
  # ggsave(paste0("~/", subtype, "_cell_controbution_TEMPUS_boxplot.png"), plot = curr_fig, width = 8, height = 6, bg='white',dpi = 300)
 
   
   # Step 7: calculate the paired diff between cell types grandular-squamous
   data_wide$diff <- data_wide$`malignant glandular epithelial cell of lung`-data_wide$`malignant squamous epithelial cell of lung`
   ggplot(data_wide, aes(x = diff)) +
     geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
     geom_density(color = "darkblue", size = 1) +
     labs(title = "Distribution of diff=(gradular-squamous)", x = "cell contribution", y = "Density") +
     theme_minimal()
               
}



# make a violin plot/ histogram to TEMPUS seperately with p-value (for LUAD and LUSC)
# make a violin plot/ histogram to TCGA seperately with p-value (for LUAD and LUSC)
# make graphs also for all samples with p-value (for LUAD and LUSC)
# look on the purity




######################################
# Boxplot per dataset with statistics
######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Inner function: Generate boxplot and calculate Wilcoxon signed-rank test
create_boxplot <- function(data_wide, curr_data, dataset) {
  # Calculate the paired Wilcoxon signed-rank test
  wilcox_test <- wilcox.test(data_wide[[2]], data_wide[[3]], paired = TRUE)
  p_value <- wilcox_test$p.value
  
  # Create the boxplot with p-value annotation
  curr_data <- curr_data %>%
    mutate(long_display_name = gsub("epithelial", "\nepithelial", long_display_name))
  
  plot <- ggplot(curr_data, aes(x = long_display_name, y = value)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(width = 0.2, color = "blue", alpha = 0.2) +  # Plot each data point
    labs(x = "cell type", y = "cell contribution", title = dataset) +
    theme_minimal() +
    annotate("text",  x = -Inf, y = Inf, 
             label = paste0("  P-value: ", sprintf("%.2e", p_value)),
             hjust = 0, vjust = 1, size = 4) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.title.x = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.title.y = element_text(angle = 90, size = 10))   
  
  return(plot)
}

generate_plots <- function(data) {
  
  for (group_label in unique(data$group_label)){
    plot_list <- list()
    
    for (dataset in unique(data$dataset_id)) {
      curr_data <- data[data$group_label==group_label & data$dataset_id==dataset,]
      # Pivot the table
      data_wide <- curr_data %>%
        pivot_wider(
          id_cols = sample_id,
          names_from = long_display_name,
          values_from = value
        )
      stats <- curr_data %>%
        group_by(long_display_name) %>%
        summarise(
          median = median(value, na.rm = TRUE),
          mean = mean(value, na.rm = TRUE),
          max = max(value, na.rm = TRUE),
          min = min(value, na.rm = TRUE)
        )
      stats$group_label <- group_label
      stats$dataset_id <- dataset
      print(stats)
      plot <- create_boxplot(data_wide, curr_data, dataset)
      plot_list[[dataset]] <- plot
    }
    
    # Arrange all subplots in a grid and save
    final_plot <- wrap_plots(plot_list, ncol = 5) +
      plot_annotation(title = paste(group_label), theme = theme(plot.title = element_text(size = 25)))
    
    # ggsave(paste0("~/", group_label, "_boxplot_per_dataset.png"), plot = final_plot, width = 22, height = 15)
  }
}

generate_plots(data)



######################################
# TEMPUS purity - cell abundance correlation
######################################

data_tempus <- read.csv('~/p01-mnsclc-tempus-data/bq-results-20241113-124928-1731502194741.csv')
data_tempus_filtered <- data_tempus[! duplicated(data_tempus[,c('feature_id', 'sample_id', 'value', 'long_display_name')]), ]


# load the annotation file
final_annotations <- read.csv('Final annotations.csv')

# merge the annotations of purity and
data_anno <- merge(data_tempus_filtered[, c('feature_id', 'sample_id', 'value', 'dataset_id', 'long_display_name')], 
                  final_annotations[, c('sample_id', 'organ', 'condition_subtype', 'purity', 'tnm_stage')], by="sample_id")
# rename to LUAD and LUSC
data_anno$condition_subtype[data_anno$condition_subtype=='lung adenocarcinoma'] = 'LUAD'
data_anno$condition_subtype[data_anno$condition_subtype=='lung squamous cell carcinoma'] = 'LUSC'

# remove samples that are not LUAD/LUSC
data_anno <- data_anno[data_anno$condition_subtype %in% c('LUAD', 'LUSC'),]


head(data_anno)




library(ggplot2)
library(ggpubr)  # For ggscatter and easy correlation calculation


scatter_plot_with_stats <- function(data, x, y, color = NULL, facet_wrap = NULL, title = NULL) {
  
  # Base plot setup
  plot <- ggplot(data, aes_string(x = x, y = y, color = color)) +
    geom_point(alpha=0.4) +
    geom_smooth(method = "lm", linetype = "dashed", color = "red") +  # Add fitted line
    theme_minimal() +
    labs(x = "Tumor purity", y = "cell contribution", title=title)
  
  # Calculate correlation and p-value for each facet separately, if facet_wrap is provided
  if (!is.null(facet_wrap)) {
    # Calculate correlation and p-value for each facet group
    data_stats <- data %>%
      group_by(.data[[facet_wrap]]) %>%
      summarize(
        cor_value = round(cor(.data[[x]], .data[[y]], use = "complete.obs"), 2),
        p_value = cor.test(.data[[x]], .data[[y]], use = "complete.obs")$p.value,
        .groups = "drop"
      )
    
    # Add text annotations to the plot for each facet
    plot <- plot +
      facet_wrap(as.formula(paste("~", facet_wrap))) #+
      # geom_text(data = data_stats,
      #           aes(label = sprintf("Correlation: %.2f\nP-value: %.2e", cor_value, p_value)),
      #           x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, inherit.aes = FALSE, color = "black", size = 3)
      
  } else {
    # If no facet_wrap, calculate overall correlation and p-value
    correlation <- cor.test(data[[x]], data[[y]], use = "complete.obs")
    cor_value <- round(correlation$estimate, 2)
    p_value <- correlation$p.value
    cor_text <- sprintf("Correlation: %.2f\nP-value: %.4f", cor_value, p_value)
  }
  # Add text annotation for overall correlation
  plot <- plot +
    annotate("text",  x = -Inf, y = Inf, 
             label = sprintf("Correlation: %.2f\nP-value: %.2e", cor_value, p_value),
             hjust = 0, vjust = 1, size = 4)
  
  # Display the plot
  return(plot)
}

generate_plots <- function(data) {
  plot_list <- list()
  for (group_label in unique(data$organ)){
    
    # curr_data <- data[data$organ==group_label,]
    # stats <- curr_data %>%
    #   group_by(long_display_name) %>%
    #   summarise(
    #     median = median(value, na.rm = TRUE),
    #     mean = mean(value, na.rm = TRUE),
    #     max = max(value, na.rm = TRUE),
    #     min = min(value, na.rm = TRUE)
    #   )
    # stats$group_label <- group_label
    print(stats)
    plot <- scatter_plot_with_stats(data=curr_data, x='purity', y='value', color='condition_subtype', facet_wrap='long_display_name')
    plot <- plot +  labs(title = group_label)
    plot_list[[group_label]] <- plot

    
    # Arrange all subplots in a grid and save
  final_plot <- wrap_plots(plot_list, ncol = 5) +
      plot_annotation(title = paste(group_label), theme = theme(plot.title = element_text(size = 25)))
    
    # ggsave(paste0("~/", group_label, "_boxplot_per_dataset.png"), plot = final_plot, width = 22, height = 15)
  }
}


###### scatter plot of purity vs cell contribution

# LUAD/LUSC in different colors
scatter_plot_with_stats(data=data_anno, x='purity', y='value', color='condition_subtype', facet_wrap='long_display_name')
scatter_plot_with_stats(data=data_anno[data_anno$organ=='lung',], x='purity', y='value', color='condition_subtype', facet_wrap='long_display_name')
generate_plots(data_anno)

# separate plots for LUAD and LUSC
scatter_plot_with_stats(data=data_anno[data_anno$condition_subtype=='LUAD',], x='purity', 
                        y='value', facet_wrap='long_display_name', title='LUAD\nAll samples (including metastasis)')
scatter_plot_with_stats(data=data_anno[data_anno$condition_subtype=='LUSC',], x='purity', 
                        y='value', facet_wrap='long_display_name', title='LUSC\nAll samples (including metastasis)')

scatter_plot_with_stats(data=data_anno[data_anno$condition_subtype=='LUAD' & data_anno$organ=='lung',], x='purity', 
                        y='value', facet_wrap='long_display_name', title='LUAD\nPrimary tumor samples')
scatter_plot_with_stats(data=data_anno[data_anno$condition_subtype=='LUSC' & data_anno$organ=='lung',], x='purity', 
                        y='value', facet_wrap='long_display_name', title='LUSC\nPrimary tumor samples')


scatter_plot_with_stats(data=data_anno[data_anno$condition_subtype=='LUSC' & data_anno$organ=='lung',], x='purity', 
                        y='value', facet_wrap='long_display_name', title='LUSC\nPrimary tumor samples')


##### sum the 4 cells to one cell type
data_sum <- data_anno %>% group_by(sample_id) %>% summarise(sum = sum(value))

data_sum_merged <- merge(data_sum, 
      final_annotations[, c('sample_id', 'organ', 'condition_subtype', 'purity', 'tnm_stage')], by="sample_id")

data_sum_merged$long_display_name <- 'cell'
data_sum_merged$condition_subtype[data_sum_merged$condition_subtype=='lung adenocarcinoma'] = 'LUAD'
data_sum_merged$condition_subtype[data_sum_merged$condition_subtype=='lung squamous cell carcinoma'] = 'LUSC'

scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUSC' & data_sum_merged$organ=='lung',], x='purity', 
                        y='sum', title='LUSC\nPrimary tumor samples')
scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUAD' & data_sum_merged$organ=='lung',], x='purity', 
                        y='sum', title='LUAD\nPrimary tumor samples')


scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUSC',], x='purity', 
                        y='sum', title='LUSC\nAll samples (including metastasis)')
scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUAD',], x='purity', 
                        y='sum', title='LUAD\nAll samples (including metastasis)')



##### sum the 2 malignant cells to one cell type
unique(data_anno$long_display_name)
data_anno <- data_anno[data_anno$long_display_name %in% c("malignant squamous epithelial cell of lung", "malignant glandular epithelial cell of lung"),]
unique(data_anno$long_display_name)
dim(data_anno)
data_sum <- data_anno %>% group_by(sample_id) %>% summarise(sum = sum(value))
head(data_sum)

data_sum_merged <- merge(data_sum, 
                         final_annotations[, c('sample_id', 'organ', 'condition_subtype', 'purity', 'tnm_stage')], by="sample_id")

data_sum_merged$condition_subtype[data_sum_merged$condition_subtype=='lung adenocarcinoma'] = 'LUAD'
data_sum_merged$condition_subtype[data_sum_merged$condition_subtype=='lung squamous cell carcinoma'] = 'LUSC'

# LUSC- primary+metastasis samples
plot <- scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUSC',], x='purity', 
                        y='sum', title='LUSC\nPrimary and metastasis samples')
plot <- plot + labs(x = "Tumor purity", y = "Malignant cell contribution")
ggsave(paste0("~/", "TEMPUS_LUSC_purity_cell_controbution_corr.png"), plot = plot, width = 8, height = 6, bg='white',dpi = 300)

# LUAD- primary+metastasis samples
plot <- scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUAD',], x='purity', 
                                y='sum', title='LUAD\nPrimary and metastasis samples')
plot <- plot + labs(x = "Tumor purity", y = "Malignant cell contribution")
ggsave(paste0("~/", "TEMPUS_LUAD_purity_cell_controbution_corr.png"), plot = plot, width = 8, height = 6, bg='white',dpi = 300)

# LUSC- primary samples
plot <- scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUSC' & data_sum_merged$organ=='lung',], x='purity', 
                                y='sum', title='LUSC\nPrimary tumor samples')
plot <- plot + labs(x = "Tumor purity", y = "Malignant cell contribution")
ggsave(paste0("~/", "TEMPUS_LUSC_purity_cell_controbution_corr_primary.png"), plot = plot, width = 8, height = 6, bg='white',dpi = 300)

# LUAD- primary samples
plot <- scatter_plot_with_stats(data=data_sum_merged[data_sum_merged$condition_subtype=='LUAD' & data_sum_merged$organ=='lung',], x='purity', 
                                y='sum', title='LUAD\nPrimary tumor samples')
plot <- plot + labs(x = "Tumor purity", y = "Malignant cell contribution")
ggsave(paste0("~/", "TEMPUS_LUAD_purity_cell_controbution_corr_primary.png"), plot = plot, width = 8, height = 6, bg='white',dpi = 300)






###### boxplot of purity distribution
create_boxplot <- function(curr_data, title=NULL) {
   # Create the boxplot with p-value annotation
  # curr_data <- curr_data %>%
  #   mutate(long_display_name = gsub("epithelial", "\nepithelial", long_display_name))
  
  plot <- ggplot(curr_data, aes(x = condition_subtype, y = purity)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(width = 0.2, color = "blue", alpha = 0.2) +  # Plot each data point
    labs(x = "", y = "Tumor purity", title = title) +
    theme_minimal() +
    # annotate("text",  x = -Inf, y = Inf, 
    #          label = paste0("  P-value: ", sprintf("%.2e", p_value)),
    #          hjust = 0, vjust = 1, size = 4) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.title.x = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.title.y = element_text(angle = 90, size = 10))   
  
  return(plot)
}



data <- data_anno[! duplicated(data_anno[,c('sample_id', 'organ', 'condition_subtype', 'purity')]), ]
plot <- create_boxplot(data, title = 'TEMPUS - Primary and metastasis samples')
ggsave(paste0("~/", "TEMPUS_purity_LUAD_LUSC.png"), plot = plot, width = 8, height = 6, bg='white',dpi = 300)

plot <- create_boxplot(data_anno[data_anno$organ=='lung',], title = 'TEMPUS - Primary tumor samples')


###### hoistogram of purity distribution
ggplot(data_anno, aes(x = purity)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(title = "Distribution of percent purity", x = "Purity", y = "Density") +
  theme_minimal() +
  facet_wrap(as.formula(paste("~", 'condition_subtype')))


