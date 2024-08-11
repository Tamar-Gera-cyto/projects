library(cytoreason.cc.client)
load_all("/home/coder/service-single-cell/")
library(cytoreason.single.cell, lib.loc = "/opt/R/4.2.3/lib/R/library")

# read the config which was an input for wf-f33e74cf3e (V9) in BQ: SELECT * FROM `cytoreason.p00_io_colorectal_cancer.upload_admin` where service_name='single_cell' order by version desc
config_tmp <- readRDS(get_task_inputs('wf-f33e74cf3e',"0", files_names_grepl_pattern = "config"))
config_ns <- readRDS(get_task_inputs('wf-e3c410ef58',"0", files_names_grepl_pattern = "config"))

# convert the "sample_annotations" to csv - these are the annotations that were used instead of the curator
pelka_sample_annotations <- as.data.frame(config_tmp$sample_annotations)
write.csv(pelka_sample_annotations, 'pelka_sample_annotations.csv')

pelka_sample_annotations <- read.csv('pelka_sample_annotations.csv')

# add a column "subtype" with the values of "mmrstatus" (because of nan values in "contrast_effect", "contrast_effect_group_a", "contrast_effect_group_b" of MMRd_vs_MMRp
pelka_sample_annotations$subtype = pelka_sample_annotations$mmrstatus
print(pelka_sample_annotations$subtype)
write.csv(pelka_sample_annotations, '~/pelka_sample_annotations.csv')

config <- sc_api_configuration(
  # asset_id [mandatory]: workflow ID of an exported single-cell dataset as
  # output by `sc_api_save_data()`
  asset_id = "wf-803f4fdde3",  
  cell_annotation = "clontology_merged_new",
  # experiment_id [optional]: experiment identifier.
  # If not provided, the default is to use the asset_id.
  # It does not need to be a GSE ID.
  # Sample annotations are looked up in the Currator using this key, but the
  # pipeline will not fail if it does not find any
  experiment_id = "pelka_2021_colon",
  # platform_id [optional]: platform ID, if not provided it will be "sc-rnaseq"
  # Note however that the platform ID in the text files used by the platform API will
  # always be set to "sc-rnaseq".
  platform_id = "sc-rnaseq",
  # biosample_id [optional]: a string that indicates the name of the cell annotation
  # variable (in the cell metadata), which holds the bio-sample identifier,
  # i.e. the key that is used to aggregate the single-cell data into a
  # pseudobulk data
  biosample_id = "patienttypeid",
  sample_annotations = read_data("~/pelka_sample_annotations.csv"),
  
  term_metadata = list(
    tumor = list(
      tumor_vs_tumor_adjacent = list(
        contrast_type = "disease_vs_control:adjacent"
        # term_is_io = TRUE #?#
      )
    ),
    MMRd_vs_MMRp = list(
      MMRd_vs_MMRp = list(
        contrast_type = "disease_vs_disease"#"inflamed_vs_non_inflamed", #"inflamed_vs_non.inflamed" 
        # term_is_io = TRUE #?#
      )
    )
  ),
  
  analysis_model = list(
    tumor = list(
      group = list(
        # this will fit the same term as use-case 2, but label it as "DZ_vs_HC"
        sample_classification = c(tumor_adjacent = "Tumor adjacent (tumor surrounding)", tumor = "Tumor")
      )
      # pairing = NULL,
      # covariates = c("assay","layer")
      # covariates = list()
      # covariates = NULL
      # covariates = "layer"
    ),
    MMRd_vs_MMRp = list(
      group = list(
        subtype = c(MMRp = "mmrp", MMRd = "mmrd")
      )
      # pairing = NULL,
      # covariates = c("assay","layer")
      # covariates = list()
      # covariates = NULL
      # covariates = "layer"
    )
  ),
  parameters = list(
    model=list(
      services = c("gx_diff", "ct_test"),
      gx_diff=list(),  
      # min_sample_size=15,
      ct_test=list()
    )
  )
)



# check the configuration locally
sc_validate_dataset_configuration(config)
# check that a valid dataset can be loaded from the config
pb_object <- prepare_single_cell_dataset(config)
# check a particular analysis model (if failing)
obj <- prepare_single_cell_dataset(config, validate = FALSE)

str(build_analysis_model_metadata(config$analysis_model$speciman_comparison, obj))
#str(build_analysis_model_metadata(config$analysis_model$Treatment, obj))
str(build_analysis_model_metadata(config$analysis_model$treatment, obj))
#str(build_analysis_model_metadata(config$analysis_model$Response, obj))
str(build_analysis_model_metadata(config$analysis_model$response, obj))


## Available model keys
knitr::kable(ldply(get_disease_model_metadata(), .id = "model_key", as.data.frame))

# Running the pipeline
# launch piepline
sc_api_run_analysis("CRC-CO", 
                    config, 
                    image = "eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package:SUP_4864_condition_subtype_in_columns_list_latest")
# [Tue Jul 23 09:01:44 2024] run_analysis - wf-5432744b6b
# [Tue Jul 23 09:09:55 2024] run_analysis - wf-fdff501873
# [Tue Jul 23 10:10:31 2024] run_analysis - wf-2f4d8b90f0 - run again (I changed back to "inflamed_vs_non_inflamed")
# [Wed Jul 24 10:40:11 2024] run_analysis - wf-eaf72232fe - final fix for now, fixing tumor_vs_tumor_adjacent comparison
# [Tue Jul 30 17:15:33 2024] run_analysis - wf-bee0288bb8 - add term_is_io=TRUE (Yacov asked)
# [Wed Jul 31 10:32:00 2024] run_analysis - wf-a01ad0330f - MT-3822: changing the effect ids of the comparisons from specimen_comparison and MMRStatus_comparison to tumor and MMRd_vs_MMRp respectively.
# [Thu Aug  1 13:31:46 2024] run_analysis - wf-782fc6d4d0 - correcting contrast type of MMRp vs MMRd to disease_vs_disease and mapping the subtype column with mmrstatus values

model_md <- get_disease_model_metadata("CRC-CO")


# BQ upload ----
command <- paste("python /app/cytobigquery/exec_service.py single_cell", 
                 wf_id = "wf-782fc6d4d0", 
                 dataset_name = "p00_io_colorectal_cancer", #"p01_io_nsclc", 
                 "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE,
    data_access = 'public'
  )
res
# Cyto-CC workflow: wf-fc082bb548 - CRC SC after fixing tumor_Vs_tumor adjacent comparison (used to be tumor_vs_normal)
# Cyto-CC workflow: wf-f24c76e512 - CRC SC after adding term_is_io=TRUE (Yacov asked)
# Cyto-CC workflow: wf-2a224a663c - CRC SC MT-3822: changing the effect ids of the comparisons from specimen_comparison and MMRStatus_comparison to tumor and MMRd_vs_MMRp respectively.
# Cyto-CC workflow: wf-72b4433314 - correcting contrast type of MMRp vs MMRd to disease_vs_disease and mapping the subtype column with mmrstatus values


