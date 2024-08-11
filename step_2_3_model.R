####### Step 2 - Generate text output #######
# install.packages('cytoreason.cc.client')
library(cytoreason.cc.client)
# install.packages('cytoreason.ccm.pipeline')
library(cytoreason.ccm.pipeline)

disease_name <- 'Colorectal cancer'
ccm_wf_id <- "wf-7bf34a47ca" # V57 in Tableau
  # "wf-4894f31af4" # V56 in Tableau
  # "wf-08a2575e53" #V54 in Tableau
  # "wf-261a9efac0" #- V47
  # "wf-f29bf67dc9"
  # "wf-19411b52c0"
  # "wf-9ea6cf3dae"
  # "wf-6873269f1e"
  # "wf-6d43a4628e" 
  # "wf-913ab7ab26"
  # "wf-182bbe7410"
  #"wf-526fa2c5fa" 
  #"wf-927e6e1284"
  #"wf-eea4eae0a5"
memory_request <- '32Gi'
memory_request_downstream <- '[20Gi]'
# classifier <- "V3"
# special_details <- "check subtypes TCGA - with IO flag"
# TAGS$comment <- "check subtypes TCGA - with IO flag"
wf <- ccm_api_save_data(AssetData(ccm_wf_id),
                        DirectoryTXT(), #bio_qc = F,
                        image = ccm_cyto_cc("eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest",
                                            save_data.api = memory_request_downstream),
                        # image = ccm_cyto_cc("eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:feature_ENGINE_2044_support_json_config_0.57.17",
                        #                     save_data.api = memory_request_downstream),
                        # 
                        memory_request=memory_request,
                        tags=list(list(name="disease", value=disease_name),
                                  list(name="ccm-wf-id", value=ccm_wf_id)
                                  # list(name="classifier", value=classifier),
                                  # list(name="comment", value=TAGS$comment),
                                  #  list(name="project", value=TAGS$project),
                                  # list(name="special_details", value=special_details)
                        )
)

# save_data -- Tue May 28 14:39:54 2024: wf-0f0b3e1556
# save_data -- Sun Jun  2 11:59:14 2024: wf-eeed0ceb44 - adding stage comparisons 
# save_data -- Sun Jun  2 20:11:16 2024: wf-88fff233f8 - adding mss/msi comparisons for TCGA
# save_data -- Mon Jun  3 18:51:42 2024: wf-4a25915156 - correcting MSI/MSH comparisons of TCGA to be on primary tumor
# save_data -- Thu Jun  6 08:09:26 2024: wf-2bf3c97568 - GSE183984: dividing treatment_cetuximab post vs pre into 2 separate comparisons: primary_tumor_post_vs_pre and metastasis_post_vs_pre, and changing from Paired to Group
# save_data -- Sun Jun  9 19:54:14 2024: wf-075f4ca607 - including GSE195985 that failed ccm-run because of multiple SRRS, now by changing asset_id to the manually uploaded workflow:wf-329abf986f after Ori's correction 
# save_data -- Mon Jun 17 04:36:45 2024: wf-723b52a279 - adding GSE183984 cetuximab response_post_EGFR primary and metastasis separately and removing one TCGA MMRd_vs_MMRp comparison (keeping the one with more excluded samples: with 468 samples excluded vs 454)
# save_data -- Wed Jun 19 04:42:01 2024: wf-a581a4d905 - changing microsatellite stability to instability, bringing back MMRp_vs_MMRd of TCGA, changing "metastatic_colonadeno" to stage_IV and excluding samples (based on meta PCA):SRR11296798 from GSE146889, SRR6191654 from GSE104836, GSM5857288 from GSR195985
# save_data -- Thu Jun 20 04:18:37 2024: wf-1f4b933112 - adding "_" in the microsatellite instability comparisons and changing
# save_data -- Fri Jun 28 12:25:03 2024: wf-d7113acd3c - V46 after KOL comments
# save_data -- Tue Jul  2 14:23:00 2024: wf-99cf8bdfa9 - V47 excluding some of the comparisons, adding colon_stage_IV_vs_other and colon_stage_IV_vs_other, excluding some samples
# save_data -- Fri Jul 26 08:42:57 2024: wf-6c96c91d00 - V54 in Tableau - ran again because of the duplication in runs SUP-4892, also removing comparisons that were excluded in export: 'rectum_msi_h_vs_mss','rectum_msi_h_vs_msi_l','rectum_msi_l_vs_mss','rectum_msi_vs_mss' and 'colon_stage_IV_msi_l_vs_mss' that was added to overcome caching
# save_data -- Sat Aug  3 03:00:34 2024: wf-eb25e6ba4f - V56 in Tableau - fixing after first QA (GSE146889 - MSI_H_tumor and MSS_tumor taking only colorectal samples, TCGA-unite rectum_tumor_vs_normal and colon_tumor_vs_normal to tumor_vs_normal)
# save_data -- Sat Aug  3 11:45:28 2024: wf-1ff8dab65d - running agian because of the error in export v56
# save_data -- Sun Aug  4 11:54:24 2024: wf-040b2becdc - V57 in Tableau - re-dividing tumor_vs_normal of TCGA to colon/rectum and changing the name effect id of treatment_cetuximab to contain primary/metastasis information (makes term md creation easier)

####### Step 3 - upload text output #######

library(cytoreason.cc.client)

#########
# image #
#########
task_image = 'eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest'
memory_request = '64Mi'

#####################
# ARGS (positional) #
#####################
service = 'ccm'
dataset = "p00_io_colorectal_cancer" ## !!!!!!!! Pay attention to the DB name !!!!!!!! #
wf_id = "wf-040b2becdc"


###################
# ARGS (optional) #
###################
verbose = '--verbose'
ti = sprintf('-ti=%s', task_image)
tm = sprintf('-tm=%s', memory_request)


###########
# command #
###########
command <- sprintf('python /app/cytobigquery/exec_service.py %s %s %s %s %s %s',
                   service, wf_id, dataset, verbose, ti, tm)
############
# workflow #
############
tags=list(list(name="de-process", value='bigquery-upload'),
          list(name="service", value=service),
          list(name="export_workflow", value=wf_id),
          list(name="target_dataset", value=dataset)#,
          # list(name="project", value=TAGS$project)
)
task_env_vars = list(list("name"="DE_PROCESS","value"="BigQuery_upload")
)
wf <- run_command_dist(command,
                       outdir = "./output/",
                       image = task_image ,
                       task_env_vars = task_env_vars,
                       tags=tags,
                       memory_request = memory_request)
wf
# Cyto-CC workflow: wf-9e9e93017c
# Cyto-CC workflow: wf-33178bb067 - adding stage comparisons
# Cyto-CC workflow: wf-8119d5c3fd - adding mss/msi comparisons for TCGA
# Cyto-CC workflow: wf-7bcb7e4144 - correcting MSI/MSH comparisons of TCGA to be on primary tumor
# Cyto-CC workflow: wf-918e1a2f42 - GSE183984: dividing treatment_cetuximab post vs pre into 2 separate comparisons: primary_tumor_post_vs_pre and metastasis_post_vs_pre, and changing from Paired to Group
# Cyto-CC workflow: wf-3b4866dcbe - including GSE195985 that failed ccm-run because of multiple SRRS, now by changing asset_id to the manually uploaded workflow:wf-329abf986f after Ori's correction 
# Cyto-CC workflow: wf-603796b966 - adding GSE183984 cetuximab response_post_EGFR primary and metastasis separately and removing one TCGA MMRd_vs_MMRp comparison (keeping the one with more excluded samples: with 468 samples excluded vs 454)
# Cyto-CC workflow: wf-049a86b375 - V43 changing microsatellite stability to instability, bringing back MMRp_vs_MMRd of TCGA, changing "metastatic_colonadeno" to stage_IV and excluding samples (based on meta PCA):SRR11296798 from GSE146889, SRR6191654 from GSE104836, GSM5857288 from GSR195985
# Cyto-CC workflow: wf-b0c4b8d863 - V44 adding "_" in the microsatellite instability comparisons and changing
# Cyto-CC workflow: wf-0df366784a - V46 afeter KOL comments
# Cyto-CC workflow: wf-dc1feb3e3f - V47 - V47 excluding some of the comparisons, adding colon_stage_IV_vs_other and colon_stage_IV_vs_other, excluding some samples
# Cyto-CC workflow: wf-83f74466ce - V54 in Tableau - ran again because of the duplication in runs SUP-4892, also removing comparisons that were excluded in export: 'rectum_msi_h_vs_mss','rectum_msi_h_vs_msi_l','rectum_msi_l_vs_mss','rectum_msi_vs_mss' and 'colon_stage_IV_msi_l_vs_mss' that was added to overcome caching
# Cyto-CC workflow: wf-80ecc7dc0b - run again because of error in ETL trigger
# Cyto-CC workflow: wf-45b799b6ab - again with wf-eb25e6ba4f - worked!
# Cyto-CC workflow: wf-f607557120 - V57 in Tableau - re-dividing tumor_vs_normal of TCGA to colon/rectum and changing the name effect id of treatment_cetuximab to contain primary/metastasis information (makes term md creation easier)





####### Step 4 - download config json and convert to csv #######
library(hash)
library(dplyr)
library("cytoreason.ccm.pipeline")

prepare_config <- function(config_file){
  # Edits config file so it can be an input to build_designModelTermMetadata
  config_file_edited <- config_file
  column_names_dict <- hash()
  column_names_dict <- setNames(
    c("group", "pairing", "covariates", "timepoint", "exclude_samples"), # old
    c("model_group", "model_pairing", "model_covariates", "model_timepoint", "model_exclude_samples") # new
  ) 
  column_order <- c("ccm_exclude", "ccm_meta", "ccm_meta_pca",
                    "ccm_annotate",	"dataset_id",	"experiment_id",
                    "platform_id", "asset_id", "asset_data_access",
                    "effect_id", "model_type", "comparison", 
                    "comparison_levels", "model_group", "model_pairing",
                    "model_covariates", "model_timepoint", "model_exclude_samples",
                    "comments", "sample_tissue", "tissue", "condition")
  # rename the columns
  config_file_edited <- config_file_edited %>%
    rename(!!!column_names_dict)
  # remove unnecessary columns
  config_file_edited <- config_file_edited[column_order]
  # edit model_group column to remove the outer curly brackets and double quotes 
  config_file_edited$model_group <- gsub('^\\{"(.*?):', '\\1:', config_file_edited$model_group) # removes outer curly brackets, double quotes
  config_file_edited$model_group <- gsub('\\}$', '', config_file_edited$model_group) # removes the curly bracket at the end of the expression
  config_file_edited$model_group <- gsub('"', '', config_file_edited$model_group) # removes all "
  
  return(config_file_edited)
}

library("cytoreason.ccm.pipeline")
ccm_fit <- as_ccm_fit(ccm_run_wf)
config_file <- modelMetadata(ccm_fit)
config_file <- prepare_config(config_file)
write.csv(config_file, paste0('~/config_', ccm_run_wf, '.csv'))



library("cytoreason.ccm.pipeline")
library(cytoreason.cc.client)
library(cytoreason.shared.assets)
library(cytoreason.curator.annotations)

eset <- read_asset("ccw://wf-4f0b666280:0:ccm_fit.qs")
eset$datasets$TCGA__rnaseq$cell_pca$center
a <- analysisResultCell(eset$datasets$TCGA__rnaseq)
analysisResultCell(ccm_fit$datasets$GSE7307__GPL570)



crc_eset <- as_ccm_fit('wf-2ce5bee91e')
ccm_fit <- crc_eset
crc_cell_contribution <- analysisResultCell(eset$datasets$TCGA__rnaseq)
crc_gene_set <- analysisResultExpressionSet(crc_eset$datasets$GSE44076__GPL13667, "gene_set_activity")
eset <- assayDataExpression(crc_eset$datasets$GSE44076__GPL13667)
datasets_T = grabd_data_sets(list("Tumor"="Tumor"))
apply_model_analysis(ccm_fit, .filter = function(x) varMetadataElement(x, "Tumor") %in% "Tumor", function(model, ccm_dataset){ ... ))
#split row names by column time
tumor_vs_t_adj_groups <- split(rownames(pData(eset)), pData(eset)$sample_classification)
tumor_vs_t_adj_groups <- tumor_vs_t_adj_groups[c(1,2)]
run_method_dist(method = "cytopro_service_gene_cell_expression", ns = "cytoreason.deconvolution",
                eset =eset,
                cells_eset = crc_cell_contribution,
                group_list = tumor_vs_t_adj_groups,
                image = "eu.gcr.io/cytoreason/ci-cytoreason.deconvolution-package:cytopro_v2_release_patch_latest")



library(cytoreason.deconvolution)
service_ct_contribution(eset$datasets$GSE101588__rnaseq)

eset <- as_ccm_fit('wf-2ce5bee91e')
crc_cell_contribution <- exprs(analysisResultCell(eset$datasets$TCGA__rnaseq))


x <- analysisResultElement(eset$datasets$GSE101588__rnaseq, "cell_contribution")
y <- basis(x, format="eset", normalize="cell")

x <- eset$datasets$GSE101588__rnaseq$cell_contribution



eset_old <- read_asset("ccw://wf-261a9efac0:0:ccm_fit.qs")
crc_cell_contribution_old <- exprs(analysisResultCell(eset_old$datasets$TCGA__rnaseq))

eset_old$datasets$TCGA__rnaseq$cell_pca$center

x_old <- analysisResultElement(eset_old$datasets$GSE101588__rnaseq, "cell_contribution")
y_old <- basis(x_old, format="eset", normalize="cell")


eset_old$datasets$GSE101588__rnaseq$metadata



test_model_str ="Model=p00_io_colorectal_cancer;Version=47;Service=ccm"

