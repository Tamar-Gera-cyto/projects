#update eset based on needed comparisions for datasets that are missing in the curater
#rds_temp <- readRDS(get_task_outputs("wf-efc0d328f9", "0"))
#exp_name <- "E-MTAB-7845"
library(cytoreason.cc.client)
library(cytoreason.ccm.pipeline)
library(cytoreason.shared.assets)
library(cytoreason.curator.annotations)


###################################
# TCGA metadata - MSI/MSS
###################################
eset_curr <- read_asset("ccw://wf-dd4a53405f:0:output.rds")
exp_name <- experimentID(eset_curr)
exp_name
p_data_curr <- pData(eset_curr)
dim(p_data_curr)
View(p_data_curr)
a = p_data_curr[, c('sample_id','Tissue_comment','MSI_Status')]

write.csv(p_data_curr, 'p_data.csv')

#when eset is large - removing columns that are not needed in p_data_curr
#p_data_curr <- p_data_curr[ , -which(names(p_data_curr) %in% c("group__CD_ileum_inflamed_vs_healthy","group__CD_left_colon_vs_HC","group__CD_right_colon_vs_HC","group__CD_transv_colon_vs_HC","group__CDi_INF_vs_nINF","group__CDi_vs_Healthy","group__UC_colon_INF_vs_HC","group__UC_colon_non_INF_vs_HC","group__UC_colon_vs_HC"))]
# sub_p_data_curr = p_data_curr[p_data_curr$`diagnosis:ch1`=="ulcerative colitis",]
# row_names_sub_p_data <- row.names(sub_p_data_curr)
# length(row_names_sub_p_data)
# row_names_sub_p_data <- as.list(row_names_sub_p_data)
# random_names_sub_p_data = sample(row_names_sub_p_data, 165)
# random_names_sub_p_data = as.character(random_names_sub_p_data)
# list_for_config = paste(random_names_sub_p_data, collapse=",")
# list_for_config
# annots = get_annotations(paste('experiment_id: ', "E-MTAB-2967", sep = ''))
# dim(annots)
# View(annots)
# all(annots$sample_id==row.names(p_data_curr))

# including MSI-H and MSI-L
add_column_annot_1 <- as.data.frame(paste(p_data_curr$Tissue_comment, p_data_curr$MSI_Status,sep = "_"))
print(unique(add_column_annot_1))
colnames(add_column_annot_1) <- "tissue_comment_msi_status"

p_data_curr_updated <- cbind(p_data_curr,add_column_annot_1)

# general MSI/MSS (replacing "primary tumor_MSI-H" and "primary tumor_MSI-L" with "primary tumor_MSI")
add_column_annot_2 <- as.data.frame(paste(p_data_curr$Tissue_comment, p_data_curr$MSI_Status,sep = "_"))
colnames(add_column_annot_2) <- "tissue_comment_msi_vs_mss"
add_column_annot_2$tissue_comment_msi_vs_mss[add_column_annot_2$tissue_comment_msi_vs_mss=="primary tumor_MSI-L"] = "primary tumor-MSI"
add_column_annot_2$tissue_comment_msi_vs_mss[add_column_annot_2$tissue_comment_msi_vs_mss=="primary tumor_MSI-H"] = "primary tumor-MSI"
print(unique(add_column_annot_2$tissue_comment_msi_vs_mss))

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_2)

b <- p_data_curr_updated[, c('sample_id','Tissue_comment','MSI_Status', 'tissue_comment_msi_status', 'tissue_comment_msi_vs_mss')]


pData(eset_curr) <- p_data_curr_updated


# experimentID(eset_curr)
#experimentID(eset_curr) <- exp_name
# platformID(eset_curr)
#platformID(eset_curr) <- "GPL13158" #"rnaseq"
# hist(eset_curr@assayData$exprs)
# affy::plotDensity(eset_curr@assayData$exprs)
create_an_eset_manualy <- function(eset, exp_name) {
  saveRDS(eset, file = paste0("./output/", exp_name,".RDS"))
}
tags=list(
  list(name = "service", value = "create_an_eset_manualy"),
  list(name = "experiment_id", value = exp_name),
  list(name = "client", value = "public"))
image<-'eu.gcr.io/cytoreason/ci-cytoreason.united-package:master_1.0.0'
create_an_eset_manualy_on_cytocc <- run_function_dist(create_an_eset_manualy,
                                                      eset=eset_curr,
                                                      exp_name = exp_name,
                                                      tags=tags,
                                                      memory_request = '30Gi', #This should be enough. if not, try more
                                                      image = image,
                                                      replace_image_tags = TRUE,
                                                      force_execution = FALSE,
                                                      data_access = 'public')
create_an_eset_manualy_on_cytocc
# Cyto-CC workflow: wf-53cb855970




###################################
# TCGA metadata - age bins
###################################
eset_curr <- read_asset("ccw://wf-53cb855970:0:3a21753-6265893-3b15dae.RDS")
exp_name <- experimentID(eset_curr)
exp_name
p_data_curr <- pData(eset_curr)
dim(p_data_curr)
View(p_data_curr)
a = p_data_curr[, c('sample_id','Tissue_comment','MSI_Status')]

write.csv(p_data_curr, 'p_data.csv')


# Concatenate tissue to tissue_comment_msi_status
add_column_annot_1 <- as.data.frame(paste(p_data_curr$Tissue, p_data_curr$tissue_comment_msi_status, sep = "_"))
add_column_annot_1[add_column_annot_1=="NA_NA_NA"] = NA
print(unique(add_column_annot_1))
colnames(add_column_annot_1) <- "tissue_tissue_comment_msi_status"

p_data_curr_updated <- cbind(p_data_curr,add_column_annot_1)

# Concatenate tissue to tissue_comment_msi_vs_mss 
add_column_annot_2 <- as.data.frame(paste(p_data_curr$Tissue, p_data_curr$tissue_comment_msi_vs_mss,sep = "_"))
add_column_annot_2[add_column_annot_2=="NA_NA_NA"] = NA
print(unique(add_column_annot_2))
colnames(add_column_annot_2) <- "tissue_tissue_comment_msi_vs_mss"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_2)


# Adding age_group column for under 50 and 50 and above
p_data_curr_updated$age_group <- ifelse(p_data_curr_updated$Age < 50, "young", "old")

# concatenate age_group to tissue_tissue_comment_msi_status
add_column_annot_3 <- as.data.frame(paste(p_data_curr_updated$age_group, p_data_curr_updated$tissue_tissue_comment_msi_status, sep = "_"))
add_column_annot_3[add_column_annot_3=="NA_NA_NA_NA"] = NA
print(unique(add_column_annot_3))
colnames(add_column_annot_3) <- "age_group_tissue_tissue_comment_msi_status"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_3)


# concatenate age_group to tissue_tissue_comment_msi_vs_mss
add_column_annot_3 <- as.data.frame(paste(p_data_curr_updated$age_group, p_data_curr_updated$tissue_tissue_comment_msi_vs_mss, sep = "_"))
add_column_annot_3[add_column_annot_3=="NA_NA_NA_NA"] = NA
print(unique(add_column_annot_3))
colnames(add_column_annot_3) <- "age_group_tissue_tissue_comment_msi_vs_mss"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_3)


# Adding stage_IV_vs_others column
add_column_annot_4 <- as.data.frame((p_data_curr_updated$TNM.stage.2))
add_column_annot_4[add_column_annot_4=='Stage I'] = 'other'
add_column_annot_4[add_column_annot_4=='Stage II'] = 'other'
add_column_annot_4[add_column_annot_4=='Stage III'] = 'other'
colnames(add_column_annot_4) <- "TNM.stage_IV_vs_other"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_4)

# concatenate age_group with tissue and "TNM.stage_IV_vs_other"
add_column_annot_5 <- as.data.frame(paste(p_data_curr_updated$age_group, p_data_curr_updated$Tissue, p_data_curr_updated$TNM.stage_IV_vs_other, sep = "_"))
add_column_annot_5[add_column_annot_5=="NA_NA_NA"] = NA
print(unique(add_column_annot_5))
colnames(add_column_annot_5) <- "age_group_tissue_TNM.stage_IV_vs_other"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_5)

# concatenate tissue and case.control (tumor/normal)
add_column_annot_6 <- as.data.frame(paste(p_data_curr_updated$Tissue, p_data_curr_updated$case.control, sep = "_"))
add_column_annot_6[add_column_annot_6=="NA_control"] = 'control'
print(unique(add_column_annot_6))
colnames(add_column_annot_6) <- "tissue_case.control"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_6)

# concatenate tissue, TNM.stage.2 and MSI_status
add_column_annot_7 <- as.data.frame(paste(p_data_curr_updated$Tissue, p_data_curr_updated$TNM.stage.2, p_data_curr_updated$MSI_Status, sep = "_"))
add_column_annot_7[add_column_annot_7=="NA_NA_NA"] = NA
print(unique(add_column_annot_7))
colnames(add_column_annot_7) <- "tissue_TNM.stage.2_MSI_Status"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_7)


# concatenate tissue with mmr.cat
add_column_annot_8 <- as.data.frame(paste(p_data_curr_updated$Tissue, p_data_curr_updated$mmr.cat, sep = "_"))
print(unique(add_column_annot_8))
colnames(add_column_annot_8) <- "tissue_mmr.cat"

p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_8)

p_data_curr_updated <- replace(p_data_curr_updated, p_data_curr_updated == "NA_NA", NA)
p_data_curr_updated <- replace(p_data_curr_updated, p_data_curr_updated == "NA_NA_NA", NA)

write.csv(p_data_curr_updated, '~/p_data_26062024.csv')


pData(eset_curr) <- p_data_curr_updated


create_an_eset_manualy <- function(eset, exp_name) {
  saveRDS(eset, file = paste0("./output/", exp_name,".RDS"))
}
tags=list(
  list(name = "service", value = "create_an_eset_manualy"),
  list(name = "experiment_id", value = exp_name),
  list(name = "client", value = "public"))
image<-'eu.gcr.io/cytoreason/ci-cytoreason.united-package:master_1.0.0'
create_an_eset_manualy_on_cytocc <- run_function_dist(create_an_eset_manualy,
                                                      eset=eset_curr,
                                                      exp_name = exp_name,
                                                      tags=tags,
                                                      memory_request = '30Gi', #This should be enough. if not, try more
                                                      image = image,
                                                      replace_image_tags = TRUE,
                                                      force_execution = FALSE,
                                                      data_access = 'public')
create_an_eset_manualy_on_cytocc


# Cyto-CC workflow: wf-f5d39eef49



# check the uploaded eset
eset_curr <- read_asset("ccw://wf-f5d39eef49:0:3a21753-6265893-3b15dae.RDS")
exp_name <- experimentID(eset_curr)
exp_name
p_data_curr <- pData(eset_curr)


###################################
# TCGA metadata -
# add rectum_stage_IV_vs_other and colon_stage_IV_vs_other
###################################
eset_curr <- read_asset("ccw://wf-f5d39eef49:0:3a21753-6265893-3b15dae.RDS")
exp_name <- experimentID(eset_curr)
exp_name
p_data_curr <- pData(eset_curr)
dim(p_data_curr)
View(p_data_curr)

# concatenate tissue and "TNM.stage_IV_vs_other"
add_column_annot_1 <- as.data.frame(paste(p_data_curr$Tissue, p_data_curr$TNM.stage_IV_vs_other, sep = "_"))
add_column_annot_1[add_column_annot_1=="NA_NA"] = NA
print(unique(add_column_annot_1))
colnames(add_column_annot_1) <- "tissue_TNM.stage_IV_vs_other"

p_data_curr_updated <- cbind(p_data_curr,add_column_annot_1)
dim(p_data_curr_updated)
write.csv(p_data_curr_updated, '~/p_data_02072024.csv')

pData(eset_curr) <- p_data_curr_updated

create_an_eset_manualy <- function(eset, exp_name) {
  saveRDS(eset, file = paste0("./output/", exp_name,".RDS"))
}
tags=list(
  list(name = "service", value = "create_an_eset_manualy"),
  list(name = "experiment_id", value = exp_name),
  list(name = "client", value = "public"))
image<-'eu.gcr.io/cytoreason/ci-cytoreason.united-package:master_1.0.0'
create_an_eset_manualy_on_cytocc <- run_function_dist(create_an_eset_manualy,
                                                      eset=eset_curr,
                                                      exp_name = exp_name,
                                                      tags=tags,
                                                      memory_request = '30Gi', #This should be enough. if not, try more
                                                      image = image,
                                                      replace_image_tags = TRUE,
                                                      force_execution = FALSE,
                                                      data_access = 'public')
create_an_eset_manualy_on_cytocc


# Cyto-CC workflow: wf-ed4ebb340c



# check the uploaded eset
eset_curr <- read_asset("ccw://wf-ed4ebb340c:0:3a21753-6265893-3b15dae.RDS")
exp_name <- experimentID(eset_curr)
exp_name
p_data_curr <- pData(eset_curr)
dim(p_data_curr)
colnames(p_data_curr)





















#white columns in cutaror
#annots = get_annotations(paste('experiment_id: ', "GSE193677", sep = ''))
endo <- p_data_curr$Endoscopic_active
endo[is.na(endo)] <- "Inactive"
add_column_annot <- as.data.frame(paste(p_data_curr$tissue,endo,p_data_curr$condition,sep = "_"))
colnames(add_column_annot) <- "tissue_endo_condition"
p_data_curr_updated <- cbind(p_data_curr,add_column_annot)
add_column_annot_second <- as.data.frame(paste(p_data_curr$tissue, p_data_curr$Condition_comment,p_data_curr$condition,sep = "_"))
colnames(add_column_annot_second) <- "tissue_condition_comment_condition"
p_data_curr_updated_second <- cbind(p_data_curr_updated,add_column_annot_second)
clinical <- p_data_curr$Condition_comment
# rep_clinical = c('inactive (clinical +endo)'='inactive',
#                  'inactive (clinically) moderate (endo)'='inactive',
#                  'inactive (clinically) severe (endo)'='inactive',
#                  'inactive (clinically) mild (endo)'='inactive',
#                  'inactive (clinically)'='inactive',
#
#                      'active (clinically) mild (endo)'='active',
#                      'active (clinically) moderate (endo)'='active',
#                      'active (clinically) severe  (endo)'='active',
#                      'active (clinically) inactive (endo)'='active' ,
#                      'active (clinically)'='active')
clinical[clinical=="inactive (clinical +endo)"] <- "Inactive"
clinical[clinical=="inactive (clinically) moderate (endo)"] <- "Inactive"
clinical[clinical=="inactive (clinically) severe (endo)"] <- "Inactive"
clinical[clinical=="inactive (clinically) mild (endo)"] <- "Inactive"
clinical[clinical=="inactive (clinically)"] <- "Inactive"
clinical[clinical=="active (clinically) mild (endo)"] <- "Active"
clinical[clinical=="active (clinically) moderate (endo)"] <- "Active"
clinical[clinical=="active (clinically) severe  (endo)"] <- "Active"
clinical[clinical=="active (clinically) inactive (endo)"] <- "Active"
clinical[clinical=="active (clinically)"] <- "Active"
add_column_annot_third <- as.data.frame(paste(p_data_curr$tissue, clinical,p_data_curr$condition,sep = "_"))
colnames(add_column_annot_third) <- "tissue_clinical_condition"
p_data_curr_updated_second <- cbind(p_data_curr_updated_second,add_column_annot_third)
severity <- p_data_curr$Condition_comment
severity[severity=="inactive (clinical +endo)"] <- "Inactive"
severity[severity=="inactive (clinically) moderate (endo)"] <- "Moderate"
severity[severity=="inactive (clinically) severe (endo)"] <- "Severe"
severity[severity=="inactive (clinically) mild (endo)"] <- "Mild"
severity[severity=="active (clinically) mild (endo)"] <- "Mild"
severity[severity=="active (clinically) moderate (endo)"] <- "Moderate"
severity[severity=="active (clinically) severe  (endo)"] <- "Severe"
severity[severity=="active (clinically) inactive (endo)"] <- "Inactive"
severity[severity=="moderate (endo)"] <- "Moderate"
severity[severity=="severe (endo)"] <- "Severe"
severity[severity=="mild (endo)"] <- "Mild"
severity[severity=="inactive (endo)"] <- "Inactive"
add_column_annot_fourth <- as.data.frame(paste(p_data_curr$tissue, severity,p_data_curr$condition,sep = "_"))
colnames(add_column_annot_fourth) <- "tissue_severity_condition"
p_data_curr_updated_second <- cbind(p_data_curr_updated_second,add_column_annot_fourth)
# add_column_annot_second <- as.data.frame(paste(p_data_curr$Tissu, p_data_curr$condition,sep = "_"))
# colnames(add_column_annot_second) <- "tissue_type_condition"
# p_data_curr_updated_second <- cbind(p_data_curr_updated,add_column_annot_second)
pData(eset_curr) <- p_data_curr_updated_second
unique(p_data_curr_updated_second$tissue_type_condition)
p_data_curr_updated_sliced = subset(p_data_curr_updated_second, Tissu_type=="Colon")
eset_curr_sliced = eset_curr[,rownames(p_data_curr_updated_sliced)]
eset_curr_sliced_expres = exprs(eset_curr_sliced)
##add_column_annot <- as.data.frame(paste(p_data_curr$`Factor Value[disease]`, p_data_curr$`Factor Value[clinical information]`,sep = "_"))
##colnames(add_column_annot) <- "disease_clinical_information"
##p_data_curr_updated <- cbind(p_data_curr,add_column_annot)
#add_column_annot <- as.data.frame(p_data_curr$`Factor Value[disease]`)
#colnames(add_column_annot) <- "disease"
#p_data_curr_updated <- cbind(p_data_curr,add_column_annot)
#add_column_annot_first <- as.data.frame(annots$condition)
#colnames(add_column_annot_first) <- "condition"
#p_data_curr_updated <- cbind(p_data_curr,add_column_annot_first)
add_column_annot_second <- as.data.frame(annots$drug_dose_response)
colnames(add_column_annot_second) <- "drug_dose_response"
p_data_curr_updated <- cbind(p_data_curr_updated,add_column_annot_second)
colnames(p_data_curr_updated)[colnames(p_data_curr_updated) == "Tissue"] <- "tissue_original"
dim(p_data_curr_updated)
View(p_data_curr_updated)
##p_data_curr_updated <- cbind(p_data_curr,annots)
pData(eset_curr) <- p_data_curr_updated_second
experimentID(eset_curr)
#experimentID(eset_curr) <- exp_name
platformID(eset_curr)
#platformID(eset_curr) <- "GPL13158" #"rnaseq"
hist(eset_curr@assayData$exprs)
affy::plotDensity(eset_curr@assayData$exprs)
create_an_eset_manualy <- function(eset, exp_name) {
  saveRDS(eset,file = paste0("./output/", exp_name,".RDS"))
}
tags=list(
  list(name = "service", value = "create_an_eset_manualy"),
  list(name = "experiment_id", value = exp_name),
  list(name = "client", value = "public"))
image<-'eu.gcr.io/cytoreason/ci-cytoreason.united-package:master_1.0.0'
create_an_eset_manualy_on_cytocc <- run_function_dist(create_an_eset_manualy,
                                                      eset=eset_curr,
                                                      exp_name = exp_name,
                                                      tags=tags,
                                                      memory_request = '30Gi', #This should be enough. if not, try more
                                                      image = image,
                                                      replace_image_tags = TRUE,
                                                      force_execution = FALSE,
                                                      data_access = 'public')
create_an_eset_manualy_on_cytocc