######################################
# Converting the config to use "model_subset" column
######################################
# this script checks that the samples in each group are defined the same for ecah effect id when the config is built with model_subset

library("jsonlite")
library("dplyr")
library("readr")
library(jsonlite)
library(stringr)

library(cytoreason.cc.client)
library(cytoreason.ccm.pipeline)

convert_json_config_to_csv <- function(json_path, datasets=TRUE, comparisons=TRUE, outpath="~/", prefix=""){
  json_text <- read_file(json_path)
  data <- fromJSON(json_text)
  
  # Extract and flatten the "datasets" and "comparisons" data
  if (datasets){
    df_datasets <- as.data.frame(data$data$datasets)
    dataset_filename = paste0(outpath, "/", prefix, "_", "datasets.csv")
    write.csv(df_datasets, dataset_filename, row.names = FALSE)
    print(paste0("file created: ", dataset_filename))
  }
  if (comparisons){
    df_comparisons <- as.data.frame(data$data$comparisons)
    comparisons_filename = paste0(outpath, "/", prefix, "_", "comparisons.csv")
    write.csv(df_comparisons, comparisons_filename, row.names = FALSE)
    print(paste0("file created: ", comparisons_filename))
  }
}

subset_by_col_values <- function(df, column_name, values) {
  # Convert the string column name to a symbol
  column_name <- sym(column_name)
  # Filter the data frame based on the values in the specified column
  df %>% filter((!!column_name) %in% values)
}

subset_by_parsed_data <- function(df, parsed_data) {
  # Function to subset dataframe based on parsed_data criteria
  subset_df <- df
  for (key in names(parsed_data)) {
    subset_df <- subset_by_col_values(subset_df, key, parsed_data[[key]])
  }
  return(subset_df)
}

parse_model_subset <- function(model_subset){
  # Step 1: Add double quotes around keys
  expression_json <- str_replace_all(model_subset, "([a-zA-Z0-9_]+):", '"\\1":')
  
  # Step 2: Add double quotes around each value inside square brackets, but not around the brackets themselves
  expression_json <- str_replace_all(expression_json, "\\[([^\\]]+)\\]", function(x) {
    # Remove square brackets temporarily
    contents <- gsub("\\[|\\]", "", x)
    
    # Split the contents by commas (in case there are multiple items)
    values <- str_split(contents, ",")[[1]]
    
    # Add quotes around each value individually
    quoted_values <- paste0('"', trimws(values), '"', collapse = ", ")
    
    # Reassemble with brackets around the quoted values
    paste0("[", quoted_values, "]")})
  parsed_data <- fromJSON(expression_json)
  return(parsed_data)
}


# convert config json ccm input to csv
# config with model subset:
outpath <- "/home/coder/p01-mnsclc-tempus-data"
json_path <- paste0(outpath, "/", "config_704_wf-6783e55150.json")
prefix="704_wf-6783e55150"
convert_json_config_to_csv(json_path, datasets=TRUE, comparisons=TRUE, outpath=outpath, prefix=prefix)

# config_model_subset <- read_csv(paste0(outpath, "/", prefix, "_", "comparisons.csv"))
config_model_subset <- read_csv(paste0(outpath, "/758_wf-a3f8f4a776_v35_comparisons_model_subset.csv"))
config_model_subset <- config_model_subset[config_model_subset$dataset_id=='TEMPUS',]

# config without model_subset:
outpath <- "/home/coder/p01-mnsclc-tempus-data"
json_path <- paste0(outpath, "/", "758_wf-a3f8f4a776_v35.json")
prefix="758_wf-a3f8f4a776_v35"
convert_json_config_to_csv(json_path, datasets=TRUE, comparisons=TRUE, outpath=outpath, prefix=prefix)

config <- read_csv(paste0(outpath, "/", prefix, "_", "comparisons.csv"))

# load the pdata used in the run
# datasets <- read_csv(paste0(outpath, "/", prefix, "_", "datasets.csv"))
# pdata_wf_id <- datasets$asset_id[datasets$dataset_id=='TEMPUS']
# eset <- read_asset(paste0("ccw://", pdata_wf_id, ":0:output.rds"))
# pdata <- pData(eset)
pdata <- read_csv(paste0(outpath, "/Final annotations.csv"))
names(pdata) <- gsub("\\.", "_", names(pdata))


# iterate over the effect_ids in the config of the most updated version ("config")
not_equal_effect_ids <- c()
equals_effect_ids <- c()
config_updated_with_model_subset <- config

for (i in 1:nrow(config)){
  row_i <- config[i,]
  if (row_i$dataset_id=='TEMPUS'){
    effect_id_i <- row_i$effect_id
    # print(effect_id_i)
    if (effect_id_i %in% unique(config_model_subset$effect_id)){
      # print(effect_id_i)
      # get sample list in each group - updated version (supports only 2 groups for now)
      col_name <- sub("^(.*?):.*", "\\1", row_i$model_group)
      groupA_value <- sub(".*?:(.*?)(,|\\}).*", "\\1", row_i$model_group)
      groupB_value <- sub(".*?:.*?:(.*?)(,|\\}).*", "\\1", row_i$model_group)
      groupA_samples_true <- pdata$sample_id[pdata[col_name]==groupA_value & !is.na(pdata[col_name])]
      groupB_samples_true <- pdata$sample_id[pdata[col_name]==groupB_value & !is.na(pdata[col_name])]
      
      # get sample list in each group - group defenition by model subset
      model_subset_row <- config_model_subset[config_model_subset$effect_id==effect_id_i,]
      model_subset_i <- model_subset_row$model_subset # continue here, you want to parse the model subset and then subset pdata accordingly and check the groups
      parsed_data <- parse_model_subset(model_subset_i)
      filtered_pdata <- subset_by_parsed_data(pdata, parsed_data)
      col_name <- sub("^(.*?):.*", "\\1", model_subset_row$model_group)
      groupA_value <- sub(".*?:(.*?)(,|\\}).*", "\\1", model_subset_row$model_group)
      groupB_value <- sub(".*?:.*?:(.*?)(,|\\}).*", "\\1", model_subset_row$model_group)
      groupA_samples <- filtered_pdata$sample_id[filtered_pdata[col_name]==groupA_value & !is.na(filtered_pdata[col_name])]
      groupB_samples <- filtered_pdata$sample_id[filtered_pdata[col_name]==groupB_value & !is.na(filtered_pdata[col_name])]
      
      if (!setequal(groupA_samples_true, groupA_samples) | !setequal(groupB_samples_true, groupB_samples)){
        print(paste0(effect_id_i, ": groupA- ", length(groupA_samples_true), " ", length(groupA_samples), " groupB- ", length(groupB_samples_true), " ", length(groupB_samples)) )
        not_equal_effect_ids <- append(not_equal_effect_ids, effect_id_i)
      } else {
        equals_effect_ids <- append(equals_effect_ids, effect_id_i)
        config_updated_with_model_subset[i,]$model_group <- model_subset_row$model_group
        config_updated_with_model_subset[i,]$model_subset <- model_subset_row$model_subset
      }
    }  
  }
}

config_updated_with_model_subset



