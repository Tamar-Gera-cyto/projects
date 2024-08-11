########################################################
# API
########################################################

# install.packages("devtools")
# install.packages("memoise") # version 2.0.0 and up
library("devtools")
library("memoise")
devtools::install_url("https://apps.private.cytoreason.com/sphere/client-tools/r-client/cytoreason.platform.client_1.1.1.tar.gz")
library(cytoreason.platform.client)

# cr_init_save("pyy:suhRtG3vNdZhwjZiOMZ5xPWCgqZMaMdb:eu")
cr_init("pxx:sPJnT3eNTDd00uDWLSw60DchhO4Ds9KE")
# cr_init_save("p01:Q6ZXFbHpI8B6HWAGV3JW75wbgTdbjxkU")
# cr_init_save("p00:v2zb0LG2vEqJljrRpYjFjmZpExCeS2X2:us")

Sys.setenv(PLATFORM_SERVER_URL="https://apps.private.cytoreason.com")
Sys.setenv(AUTH0_DOMAIN="cytoreason-pxx.us.auth0.com")


# Sys.setenv(PLATFORM_SERVER_URL="https://apps.private.cytoreason.com") # staging
# Sys.setenv(PLATFORM_SERVER_URL="https://apps.cytoreason.com") # production - for a model to be on production Noam has to approve it.


# query the main disease model inventory

avail_inventory <- cr_get_inventory()
# selected_model <- subset(cr_get_inventory(), model_name == "non-small cell lung carcinoma_lung")
# selected_model <- subset(cr_get_inventory(), model_name == "colorectal cancer_colon")
selected_model <- subset(cr_get_inventory(), model_name == "colorectal cancer_colon (crb)")


# retrieve a single model
m_inventory <- cr_model_get_inventory(model_name = selected_model$model_name)
duplicates <- m_inventory[duplicated(m_inventory), ]

unique(m_inventory$dataset)

m_inventory_samples <- cr_model_get_inventory_samples_metadata(model_name =selected_model$model_name)

tumor_vs_adjacent <- m_inventory[grepl("Tumor", m_inventory$effect_group_b), c("dataset", "effect", "effect_variable", "effect_group_a", "effect_group_b", "term")]
sc_terms <- m_inventory[grepl("pelka_2021_colon", m_inventory$dataset), c("dataset", "effect", "effect_variable", "effect_group_a", "effect_group_b", "term")]

cc <- cr_model_get_cell_contribution(model_name = selected_model$model_name)
cell_types <- unique(cc$feature_id)
cell_types <- cell_types[order(cell_types)]



m_inventory <- cr_model_get_inventory(model_name = selected_model$model_name)
# compare feature_id to the google bucket

config <- read.csv('~/ccm-metadata_NSCLC_p01_wf-f282d3ae0c.csv')
config <- config[is.na(config$ccm_exclude),]
unique(config$effect_id)

model_metadata_16_7 = list(
  # "LUAD_tumor", 
  # "LUSC_tumor", 
  # "LUAD_tumor_normal",
  # "subtype", 
  "tumor", 
  "tumor_normal", 
  "LUAD_LUSC", 
  "relapse", 
  "response", 
  "treatment", 
  "gender", 
  "race", 
  "LUAD_race", 
  "LUAD_relapse", 
  "smoking", 
  "LUAD_smoking", 
  "LUSC_smoking", 
  "line_of_therapy", 
  "progression", 
  "dose", 
  "treatment_combination", 
  "BOR", 
  "stage",
  "LUAD_stage", 
  "LUSC_stage", 
  "prior_therapy", 
  "non_squamous_PDL1_status", 
  "squamous_PDL1_status",
  "PDL1_status", 
  "non_squamous_TMB", 
  "TMB", 
  "non_squamous_race",           
  "squamous_race"
)

common <- c()
excluded_from_export <- c()
for (effect_id in unique(config$effect_id)){
  if (effect_id %in%  model_metadata_16_7){
    common <- append(common, effect_id)
  }
  else {
    excluded_from_export <- append(excluded_from_export, effect_id)
  }
}
not_in_config <- c()
for (effect_id in model_metadata_16_7){
  if (!(effect_id %in%  unique(config$effect_id))){
    not_in_config <- append(not_in_config, effect_id)
  }
}
length(unique(config$effect_id))
length(model_metadata_16_7)
length(common)
length(excluded_from_export)
length(not_in_config)

length(unique(m_inventory$dataset))
length(unique(config$dataset_id))
dataset_not_in_api <- c()
for (dataset in unique(config$dataset_id)){
  if (!(dataset %in%  unique(m_inventory$dataset))){
    dataset_not_in_api <- append(dataset_not_in_api, dataset)
  }
}
dataset_not_in_config <- c()
for (dataset in unique(m_inventory$dataset)){
  if (!(dataset %in%  unique(config$dataset_id))){
    dataset_not_in_config <- append(dataset_not_in_config, dataset)
  }
}




library(SeuratObject)
library(Seurat)
library(Matrix)

file_path <- "SC_Raw_Data_Files_February_2024_CRC_seurat_obj.gz"
readRDS(file_path)
data <- readMM(gzfile(file_path))

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = data)


library(cytoreason.cc.client)
eset <- read_asset("ccw://wf-4f0b666280:0:output.rds")
ccm_fit <- 
  analysisResultCell(ccm_fit$datasets$GSE7307__GPL570)










