
# # load packages
# devtools::load_all("~/cytoreason.xdesign")
# # devtools::load_all("~/cytoreason.integration")
# library("cytoreason.integration")
# devtools::load_all("~/cytoreason.ccm.pipeline")

# in R
# install.packages(c("cytoreason.xdesign", "cytoreason.integration"))
library("cytoreason.integration")
# library("cytoreason.xdesign")
devtools::load_all("~/cytoreason.xdesign")
devtools::load_all("~/cytoreason.ccm.pipeline")



model_metadata_final = list(
  "tumor", tumor=list("term" = 'NSCLC_vs_tumor_adjacent',
                      "contrast_type" = "disease_vs_control:adjacent"),
  "LUAD_tumor", LUAD_tumor=list("term" = 'LUAD_vs_LUAD_adjacent',
                                "contrast_type" = "disease_vs_control:LUAD"),
  "LUSC_tumor", LUSC_tumor=list("term" = 'LUSC_vs_LUSC_adjacent',
                                "contrast_type" = "disease_vs_control:LUSC"),
  
  
  "tumor_normal", tumor_normal=list("term" = 'NSCLC_vs_normal',
                                    "contrast_type" = "disease_vs_control"),
  "LUAD_tumor_normal", LUAD_tumor_normal=list("term" = 'LUAD_vs_normal',
                                              "contrast_type" = "disease_vs_control:LUAD"),
  "LUSC_tumor_normal", LUSC_tumor_normal=list("term" = 'LUSC_vs_normal',
                                              "contrast_type" = "disease_vs_control:LUSC"),
  
  
  "LUAD_LUSC", LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                "contrast_type" = "disease_vs_disease"),
  "PDL1_negative_LUAD_LUSC", PDL1_negative_LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                                            "contrast_type" = "disease_vs_disease"),
  "PDL1_positive_LUAD_LUSC", PDL1_positive_LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                                            "contrast_type" = "disease_vs_disease"),

  "TMB_low_LUAD_LUSC", TMB_low_LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                                "contrast_type" = "disease_vs_disease"),
  "asian_LUAD_LUSC", asian_LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                            "contrast_type" = "disease_vs_disease"),
  "white_LUAD_LUSC", white_LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                            "contrast_type" = "disease_vs_disease"),
  
  
  "relapse", relapse=list("term" = 'relapse_vs_non_relapse',
                          "contrast_type" = "disease_vs_disease"),
  "LUAD_relapse", LUAD_relapse=list("term" = 'relapse_vs_non_relapse',
                                    "contrast_type" = "disease_vs_disease"),
  "response", response=list("term" = 'responder_vs_non_responder',
                            "contrast_type" = "responder_vs_non.responder"),
  "treatment", treatment=list("term" = 'post_vs_pre',
                              "contrast_type" = "post_vs_pre"),
  
  "gender", gender=list("term" = 'male_vs_female'),
  # "contrast_type" = "disease_vs_disease"),
  "LUAD_gender", LUAD_gender=list("term" = 'male_vs_female',
                                  "contrast_type" = "disease_vs_disease"),
  "LUSC_gender", LUSC_gender=list("term" = 'male_vs_female',
                                  "contrast_type" = "disease_vs_disease"),
  
  "race", race=list('term' = c('white_vs_asian','white_vs_african_american','asian_vs_african_american')),
  #"contrast_type" = "disease_vs_disease"),
  "LUAD_race", LUAD_race=list('term' = c('white_vs_asian','white_vs_african_american', 'asian_vs_african_american'),
                              "contrast_type" = "disease_vs_disease"),
  "LUSC_race", LUSC_race=list('term' = c('white_vs_asian','white_vs_african_american'),
                              "contrast_type" = "disease_vs_disease"),
  
  "smoking", smoking=list("term" = c('former_vs_non_smoker','smoker_vs_non_smoker','smoker_vs_former'),
                          "contrast_type" = "disease_vs_disease"),
  "LUAD_smoking", LUAD_smoking=list("term" = c('former_vs_non_smoker','smoker_vs_non_smoker','smoker_vs_former'),
                                    "contrast_type" = "disease_vs_disease"),
  "LUSC_smoking", LUSC_smoking=list("term" = c('smoker_vs_former', 'smoker_vs_non_smoker', 'former_vs_non_smoker'),
                                    "contrast_type" = "disease_vs_disease"),
  
  
  "line_of_therapy", line_of_therapy=list("term" = 'second_line_vs_first_line',
                                          "contrast_type" = "drug_vs_drug"),
  "progression", progression=list("term" = 'yes_vs_no',
                                  "contrast_type" = "recurrent_vs_non.recurrent"),
  "dose", dose=list("term" = 'dose_600mg_vs_dose_300mg',
                    "contrast_type" = "drug_vs_drug"),
  "treatment_combination", treatment_combination=list("term" = c('AveOX40_vs_AveUto','AveUtoOX40_vs_AveUto','AveUtoOX40_vs_AveOX40'),
                                                      "contrast_type" = "drug_vs_drug"),
  "prior_therapy", prior_therapy=list("term" = c('recent_C.R_vs_PDx','IL1_vs_PDx','IL1_vs_recent_C.R'),
                                      "contrast_type" = "drug_vs_drug"),
  
  "PDL1_status", PDL1_status=list('term' = 'pos_vs_neg',
                                  "contrast_type" = "disease_vs_disease"),
  "LUAD_PDL1_status", LUAD_PDL1_status=list('term' = 'pos_vs_neg',
                                            "contrast_type" = "disease_vs_disease"),
  "LUSC_PDL1_status", LUSC_PDL1_status=list('term' = 'pos_vs_neg',
                                            "contrast_type" = "disease_vs_disease"),
  
  
  "best_overall_response", 
  best_overall_response=list("term" ='stable_disease_vs_progressive_disease',
                             "contrast_type" = "partial.responder_vs_non.responder"),
  best_overall_response=list("term" = 'partial_response_vs_progressive_disease',
                             "contrast_type" = "responder_vs_non.responder"),
  best_overall_response=list("term" = 'partial_response_vs_stable_disease',
                             "contrast_type" = "responder_vs_partial.responder"),
  
  
  "stage", stage=list('term' = c('stage_II_vs_stage_I','stage_III_vs_stage_I','stage_III_vs_stage_II','stage_IV_vs_stage_I','stage_IV_vs_stage_II','stage_IV_vs_stage_III')),
  # "contrast_type" = "disease_vs_disease"),
  "LUAD_stage", LUAD_stage=list('term' = c('stage_II_vs_stage_I','stage_III_vs_stage_I','stage_III_vs_stage_II','stage_IV_vs_stage_I','stage_IV_vs_stage_II','stage_IV_vs_stage_III'),
                                "contrast_type" = "disease_vs_disease"),
  "LUSC_stage", LUSC_stage=list('term' = c('stage_II_vs_stage_I','stage_III_vs_stage_I','stage_III_vs_stage_II','stage_IV_vs_stage_I','stage_IV_vs_stage_II','stage_IV_vs_stage_III'),
                                "contrast_type" = "disease_vs_disease"),
  
  
  "TMB", TMB=list('term' = 'high_vs_low',
                  "contrast_type" = "disease_vs_disease"),
  "LUAD_TMB", LUAD_TMB=list('term' = 'high_vs_low',
                            "contrast_type" = "disease_vs_disease"),
  
  
  ############## TEMPUS
  
  ### metastasis_vs_primary
  'metastasis_vs_primary', metastasis_vs_primary = list("term" = "metastasis_vs_primary",
                                                        "contrast_type" = "metastasis_vs_primary"),
  'LUAD_metastasis_vs_primary', LUAD_metastasis_vs_primary = list("term" = "metastasis_vs_primary",
                                                                  "contrast_type" = "metastasis_vs_primary"),
  'LUSC_metastasis_vs_primary', LUSC_metastasis_vs_primary = list("term" = "metastasis_vs_primary",
                                                                  "contrast_type" = "metastasis_vs_primary"),
  'LUAD_adrenal_gland_vs_primary', LUAD_adrenal_gland_vs_primary = list("term" = "adrenal_gland_vs_primary",
                                                                        "contrast_type" = "metastasis_vs_primary"),
  'LUAD_bone_vs_primary',  LUAD_bone_vs_primary = list("term" = "bone_vs_primary",
                                                       "contrast_type" = "metastasis_vs_primary"),
  'LUAD_brain_vs_primary',  LUAD_brain_vs_primary = list("term" = "brain_vs_primary",
                                                         "contrast_type" = "metastasis_vs_primary"),
  'LUAD_chest_vs_primary',  LUAD_chest_vs_primary = list("term" = "chest_vs_primary",
                                                         "contrast_type" = "metastasis_vs_primary"),
  'LUAD_liver_vs_primary',  LUAD_liver_vs_primary = list("term" = "liver_vs_primary",
                                                         "contrast_type" = "metastasis_vs_primary"),
  'LUAD_lymph_node_vs_primary',  LUAD_lymph_node_vs_primary = list("term" = "lymph_node_vs_primary",
                                                                   "contrast_type" = "metastasis_vs_primary"),
  'LUAD_pleura_vs_primary', LUAD_pleura_vs_primary = list("term" = "pleura_vs_primary",
                                                          "contrast_type" = "metastasis_vs_primary"),
  'LUSC_bone_vs_primary',  LUSC_bone_vs_primary = list("term" = "bone_vs_primary",
                                                       "contrast_type" = "metastasis_vs_primary"),
  'LUSC_chest_vs_primary',  LUSC_chest_vs_primary = list("term" = "chest_vs_primary",
                                                         "contrast_type" = "metastasis_vs_primary"),
  'LUAD_soft_tissue_vs_primary',    LUAD_soft_tissue_vs_primary = list("term" = "soft_tissue_vs_primary",
                                                                       "contrast_type" = "metastasis_vs_primary"),
  'bone_vs_primary',  bone_vs_primary = list("term" = "bone_vs_primary",
                                             "contrast_type" = "metastasis_vs_primary"),
  'brain_vs_primary',  brain_vs_primary = list("term" = "brain_vs_primary",
                                               "contrast_type" = "metastasis_vs_primary"),
  'chest_vs_primary',  chest_vs_primary = list("term" = "chest_vs_primary",
                                               "contrast_type" = "metastasis_vs_primary"),
  'LUSC_liver_vs_primary',  LUSC_liver_vs_primary = list("term" = "liver_vs_primary",
                                                         "contrast_type" = "metastasis_vs_primary"),
  'LUSC_lymph_node_vs_primary', LUSC_lymph_node_vs_primary = list("term" = "lymph_node_vs_primary",
                                                                  "contrast_type" = "metastasis_vs_primary"), 
  'LUSC_pleura_vs_primary',  LUSC_pleura_vs_primary = list("term" = "pleura_vs_primary",
                                                           "contrast_type" = "metastasis_vs_primary"),
  'liver_vs_primary', liver_vs_primary = list("term" = "liver_vs_primary",
                                              "contrast_type" = "metastasis_vs_primary"),
  'lymph_node_vs_primary', lymph_node_vs_primary = list("term" = "lymph_node_vs_primary",
                                                        "contrast_type" = "metastasis_vs_primary"),
  'pleura_vs_primary', pleura_vs_primary = list("term" = "pleura_vs_primary",
                                                "contrast_type" = "metastasis_vs_primary"),
  'soft_tissue_vs_primary', soft_tissue_vs_primary = list("term" = "soft_tissue_vs_primary",
                                                          "contrast_type" = "metastasis_vs_primary"),
  
  
  
  
  ### LUAD_vs_LUSC
  #'LUAD_LUSC', 
  'LUAD_LUSC_including_metastasis_liver_excluded', LUAD_LUSC_including_metastasis_liver_excluded = list(term = "LUAD_vs_LUSC", 
                                                                                                        "contrast_type" = "disease_vs_disease"), 
  'LUAD_LUSC_including_metastasis_liver_included', LUAD_LUSC_including_metastasis_liver_included = list(term = "LUAD_vs_LUSC", 
                                                                                                        "contrast_type" = "disease_vs_disease"),
  'african_american_LUAD_LUSC', african_american_LUAD_LUSC = list(term = "LUAD_vs_LUSC", 
                                                                  "contrast_type" = "disease_vs_disease"),
  
  
  ### responder_vs_non_responder
  'LUAD_ICI_chemotherapy_pre_responder_vs_non_responder',
  'LUAD_PD1_or_PDL1_pre_responder_vs_non_responder', 
  'LUAD_PD1_pre_responder_vs_non_responder', 
  'LUAD_heavily_pre_responder_vs_non_responder',
  'LUAD_lightly_pre_responder_vs_non_responder', 
  'LUAD_only_ICI_pre_responder_vs_non_responder',
  'LUAD_post_responder_vs_non_responder',
  'LUAD_pre_responder_vs_non_responder', 
  'LUSC_PD1_pre_responder_vs_non_responder', 
  'LUSC_ICI_chemotherapy_pre_responder_vs_non_responder', 
  'LUSC_PD1_or_PDL1_pre_responder_vs_non_responder', 
  'LUSC_lightly_pre_responder_vs_non_responder', 
  'LUSC_post_responder_vs_non_responder', 
  'heavily_pre_responder_vs_non_responder', 
  'PD1_pre_responder_vs_non_responder', 
  'LUSC_pre_responder_vs_non_responder', 
  'lightly_pre_responder_vs_non_responder',
  'responder_vs_non_responder_pre_including_metastasis_liver_excluded', 
  'responder_vs_non_responder_pre_including_metastasis_liver_included', 
  'post_responder_vs_non_responder',
  'pre_responder_vs_non_responder', 
  
  
  
  ### heavily_vs_lightly
  'LUAD_non_responder_pre_heavily_vs_lightly', LUAD_non_responder_pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                                                "contrast_type" = "disease_vs_disease"),
  'LUAD_post_heavily_vs_lightly', LUAD_post_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                      "contrast_type" = "disease_vs_disease"),
  'LUAD_pre_heavily_vs_lightly', LUAD_pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                    "contrast_type" = "disease_vs_disease"),
  'LUAD_pre_responder_heavily_vs_lightly', LUAD_pre_responder_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                                        "contrast_type" = "disease_vs_disease"),
  'LUSC_pre_heavily_vs_lightly', LUSC_pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                    "contrast_type" = "disease_vs_disease"),
  'LUSC_responder_pre_heavily_vs_lightly', LUSC_responder_pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                                        "contrast_type" = "disease_vs_disease"),
  'heavily_vs_lightly_pre_including_metastasis_liver_excluded', heavily_vs_lightly_pre_including_metastasis_liver_excluded = list(term = "heavily_vs_lightly", 
                                                                                                                                  "contrast_type" = "disease_vs_disease"),
  'heavily_vs_lightly_pre_including_metastasis_liver_included', heavily_vs_lightly_pre_including_metastasis_liver_included = list(term = "heavily_vs_lightly", 
                                                                                                                                  "contrast_type" = "disease_vs_disease"),
  'non_responder_pre_heavily_vs_lightly', non_responder_pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                                      "contrast_type" = "disease_vs_disease"),
  'post_heavily_vs_lightly', post_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                            "contrast_type" = "disease_vs_disease"),
  'pre_heavily_vs_lightly', pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                          "contrast_type" = "disease_vs_disease"), 
  'responder_pre_heavily_vs_lightly', responder_pre_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                              "contrast_type" = "disease_vs_disease"),
  'LUAD_pre_non_responder_heavily_vs_lightly', LUAD_pre_non_responder_heavily_vs_lightly = list(term = "heavily_vs_lightly", 
                                                                                       "contrast_type" = "disease_vs_disease"),
  
  
  ### post_vs_pre
  'ICI_chemotherapy_non_responder_post_vs_pre',ICI_chemotherapy_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'ICI_chemotherapy_post_vs_pre', ICI_chemotherapy_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'ICI_chemotherapy_responder_post_vs_pre', ICI_chemotherapy_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_ICI_chemotherapy_non_responder_post_vs_pre', LUAD_ICI_chemotherapy_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_ICI_chemotherapy_post_vs_pre', LUAD_ICI_chemotherapy_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_ICI_chemotherapy_responder_post_vs_pre', LUAD_ICI_chemotherapy_responder_post_vs_pre=list("term" = 'post_vs_pre',"contrast_type" = "post_vs_pre"),
  'LUAD_PD1_non_responder_post_vs_pre', LUAD_PD1_non_responder_post_vs_pre=list("term" = 'post_vs_pre',"contrast_type" = "post_vs_pre"),
  'LUAD_PD1_or_PDL1_non_responder_post_vs_pre', LUAD_PD1_or_PDL1_non_responder_post_vs_pre=list("term" = 'post_vs_pre',"contrast_type" = "post_vs_pre"),
  'LUAD_PD1_or_PDL1_post_vs_pre', LUAD_PD1_or_PDL1_post_vs_pre=list("term" = 'post_vs_pre',"contrast_type" = "post_vs_pre"),
  'LUAD_PD1_or_PDL1_responder_post_vs_pre', LUAD_PD1_or_PDL1_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_PD1_post_vs_pre', LUAD_PD1_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_PD1_responder_post_vs_pre', LUAD_PD1_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_heavily_post_vs_pre', LUAD_heavily_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_heavily_responder_post_vs_pre', LUAD_heavily_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_lightly_non_responder_post_vs_pre', LUAD_lightly_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_lightly_post_vs_pre', LUAD_lightly_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_lightly_responder_post_vs_pre', LUAD_lightly_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_only_ICI_non_responder_post_vs_pre', LUAD_only_ICI_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUAD_only_ICI_post_vs_pre', LUAD_only_ICI_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_ICI_chemotherapy_non_responder_post_vs_pre', LUSC_ICI_chemotherapy_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_ICI_chemotherapy_post_vs_pre', LUSC_ICI_chemotherapy_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_PD1_non_responder_post_vs_pre', LUSC_PD1_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_PD1_or_PDL1_non_responder_post_vs_pre', LUSC_PD1_or_PDL1_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_PD1_or_PDL1_post_vs_pre', LUSC_PD1_or_PDL1_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_PD1_post_vs_pre', LUSC_PD1_post_vs_pre=list("term" = 'post_vs_pre',"contrast_type" = "post_vs_pre"),
  'LUSC_lightly_post_vs_pre', LUSC_lightly_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'PD1_non_responder_post_vs_pre', PD1_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'PD1_or_PDL1_non_responder_post_vs_pre', PD1_or_PDL1_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'PD1_or_PDL1_post_vs_pre', PD1_or_PDL1_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'PD1_or_PDL1_responder_post_vs_pre', PD1_or_PDL1_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'PD1_post_vs_pre', PD1_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'PD1_responder_post_vs_pre', PD1_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'TMB_low_PDL1_high_post_vs_pre', TMB_low_PDL1_high_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'TMB_low_PDL1_low_post_vs_pre', TMB_low_PDL1_low_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'heavily_non_responder_post_vs_pre', heavily_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'heavily_post_vs_pre', heavily_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'heavily_responder_post_vs_pre', heavily_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'lightly_non_responder_post_vs_pre', lightly_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'lightly_post_vs_pre', lightly_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'lightly_responder_post_vs_pre', lightly_responder_post_vs_pre=list("term" = 'post_vs_pre',"contrast_type" = "post_vs_pre"),
  'non_responder_post_vs_pre', non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'only_ICI_post_vs_pre', only_ICI_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'post_vs_pre', post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'post_vs_pre_including_metastasis_liver_excluded', post_vs_pre_including_metastasis_liver_excluded=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'post_vs_pre_including_metastasis_liver_included', post_vs_pre_including_metastasis_liver_excluded=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'post_vs_pre_paired_including_metastasis_liver_excluded', post_vs_pre_paired_including_metastasis_liver_excluded=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'post_vs_pre_paired_including_metastasis_liver_included', post_vs_pre_paired_including_metastasis_liver_included=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'responder_post_vs_pre', responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre"),
  'LUSC_lightly_non_responder_post_vs_pre', LUSC_lightly_non_responder_post_vs_pre=list("term" = 'post_vs_pre', "contrast_type" = "post_vs_pre")
  
  #'response', 
  #'smoking', 
  #'white_LUAD_LUSC'
  #'gender', 
  #'LUAD_gender', 
  #'LUSC_gender', 
  #'race',
  #'LUAD_race',
  #'LUSC_race', 
  #'TMB',   
)



ccm_fit <- as_ccm_fit("dm://p01_io_nsclc:ccm:41")
# ccm_fit <- as_ccm_fit("wf-c95eb44499")
term_md_test <- designModelTermMetadata(
  "dm://p01_io_nsclc:ccm:41", #ccm_fit, #"dm://p01_io_nsclc:ccm:19", 
  # target_terms = model_metadata_final,
  group_level_variables = 'gd://p01_nsclc_v41@1_DCJGBvsF4YqBWipNhegjHg-DfxW99sQ2OVT0qTGXxY',
  group_level_variables_column = "comparison_levels_test",
  validate = FALSE
  # subset = TRUE, # will keep terms that are in model_metadata_final
  # flag_invalid_terms = FALSE # will not filter out invalid terms (for investigation)
)

# write.csv(term_md, "~/term_md.csv")
# term_md <- read_data("~/term_md.csv")

###############################################################   
# keep only effect ids that are in model_metadata_final
############################################################### 
print(dim(term_md))
# term_md <- term_md[term_md$effect_id %in% effect_ids_model_metadata_final, ]
term_md <- term_md[!term_md$.term_exclusion %in% c('error', 'filtered'), ]
print(dim(term_md))

###############################################################   
# change time "W0" and "D0" in TEMPUS -> "Pre" and "post" -> "Post"
############################################################### 
print(unique(term_md$`term:time`))
term_md$`term:time`[(term_md$`term:time`=="D0" | term_md$`term:time`=="W0" ) & !is.na(term_md$`term:time`) & term_md$dataset_id=="TEMPUS"] <- 'Pre'
term_md$`term:time`[term_md$`term:time`=="post" & !is.na(term_md$`term:time`)] <- 'Post' # for all datasets
print(unique(term_md$`term:time`))


write.csv(term_md, "~/term_md.csv")


###
# for LUAD_tumor normal and LUSC_tumor normal: contrast_type = "disease_vs_control/subtype" and also group

###############################################################   
# change contrast_type of LUAD_tumor_normal (because I get an error that multiple meta-term have the same contrast type)
# from disease_vs_control -> disease_vs_control/subtype
# and removing the meta of LUSC_tumor_normal because there are only 2 DS in this meta (if I will not remove I will still get the error)
############################################################### 
# print(dim(term_md)) 
# term_md <- term_md[!(term_md$effect_id=="LUSC_tumor_normal" & term_md$dataset=="meta"), ]
# print(dim(term_md)) 



###############################################################   
# look at columns without NULL in a text format 
############################################################### 
write_filtered_dataframe <- function(df, output_file){
  file_conn <- file(output_file, open = "w")
  for (i in 1:nrow(df)) {
    non_na_cols <- names(df)[!is.na(df[i, ])]
    non_na_values <- df[i, non_na_cols, drop = FALSE]
    writeLines(paste0("***", non_na_values[["effect_id"]], ":"), file_conn)
    for (col in non_na_cols) {
      value <- non_na_values[[col]]
      writeLines(paste0(col, ": ", value), file_conn)
    }
    # Add a blank line between rows (optional)
    if (i < nrow(df)) writeLines("", file_conn)
  }
  # Close the file connection
  close(file_conn)
}

term_md_filtered <- term_md[term_md$term_category=='contrast',]
columns_to_look <-  c("dataset", 
                      "effect_id",
                      "term_category", 
                      "term", 
                      "term_error",
                      "term_namespace", 
                      "contrast_type",
                      "contrast_effect",
                      "contrast_effect_group_a",
                      "contrast_effect_group_b"
                      )
more_columns_to_look <- grep("^(term:|group_a:|group_b:)", colnames(term_md), value = TRUE)
columns_to_look <- append(columns_to_look, more_columns_to_look)
term_md_filtered <- term_md_filtered[, columns_to_look]
write_filtered_dataframe(term_md_filtered, "output.txt")



################################    
# 3. Validate term md
################################
dim(term_md) # 814 121
term_md <- subset(term_md, is.na(term_error))
dim(term_md) # 810 121

# Drop D of D that one of their parent comparison is dropped because of sample number or not in the model
drop_DofD_w_o_parant = function(term_md){
  model_dod_terms = subset(term_md,contrast_category=="delta_of_delta")$term
  if (dim(data.frame(model_dod_terms))[1] == 0){
    print("No D of D comparisons in this model")
    return(term_md)
  }
  split_by_dots = str_split(model_dod_terms, ":")
  parant_terms = list()
  for (i in 1:length(split_by_dots)){
    split_by_vs = str_split(split_by_dots[[i]][2], "_vs_")
    parant_term_1  = paste(split_by_dots[[i]][[1]], split_by_vs[[1]][[1]], sep=":")
    parant_term_2  = paste(split_by_dots[[i]][[1]], split_by_vs[[1]][[2]], sep=":")
    parant_terms[[model_dod_terms[i]]] = c(parant_term_1, parant_term_2)
  }
  # will_be_droped = subset(term_md,!(contrast_category=="delta_of_delta") & grepl("sample size too low", term_md$term_error))$term
  model_terms = subset(term_md, !(contrast_category == "delta_of_delta"))$term
  bolean_if_drop = unlist(lapply(names(parant_terms), function(d_of_d_term){
    to_drop = TRUE
    # for (j in 1:length(parant_terms[[d_of_d_term]])){
    if (parant_terms[[d_of_d_term]][1] %in% model_terms & parant_terms[[d_of_d_term]][2] %in% model_terms){
      to_drop = FALSE
      # }
    }
    to_drop
  }))
  dod_to_drop = names(parant_terms[bolean_if_drop])
  term_md_new = subset(term_md, !(term %in% dod_to_drop))
  term_md_new
}

canonical_contrast_types_test = function(df){
  sliced_df = subset(df, !(contrast_category=="delta_of_delta") & dataset %in% c('meta'))
  sliced_df = sliced_df[, c("meta_term_id", "term", "contrast_type", "contrast_effect_group_a")]
  saved_contrast_types = list("disease_vs_control", "disease_vs_control:adjacent",
                              "disease_vs_control:benign", "disease_vs_disease:metastases",
                              "disease_vs_disease:left_colon_vs_right_colon", "response:chemotherapy")
  test = lapply(saved_contrast_types, function(one_contrast){
    if (one_contrast %in% unique(sliced_df$contrast_type)){
      sliced_df = subset(sliced_df, contrast_type==one_contrast)
      unique_term = unique(sliced_df["meta_term_id"])
      if (dim(unique_term)[1]==1){
        print(paste0(one_contrast ," is OK"))
      }else stop(paste0(one_contrast ," have more than one meta-term, please fix"))
    }
  })
}


term_md_filtered = drop_DofD_w_o_parant(term_md)
dim(term_md_filtered) # 810 121
# filter out meta analysis of 2 data sets
# term_md_filtered <- filter_small_meta_term_id(term_md_filtered, threshold = 2L)
# dim(term_md_filtered)

validate_term_metadata(term_md_filtered)
canonical_contrast_types_test(term_md_filtered)


# removing LUAD_tumor_normal and LUSC_tumor_normal meta because of the error that there is more than 1 meta with disease_vs_control 
term_md_filtered <- term_md_filtered[!term_md_filtered$term_id %in% c('a3479c5','e203d77'),]
term_md_filtered$bq_version <- "dm://p01_io_nsclc:ccm:39"
canonical_contrast_types_test(term_md_filtered)


dim(term_md_filtered) # 810 121
View(term_md_filtered)



################################    
# 4. Submit
################################
ccm_api_save_data(term_md_filtered,
                  image = ccm_cyto_cc(save_data.api = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:SCAL_858_new_stage_term_meta_data_0.67.1.2"))
# save_data -- Sun Nov 24 11:56:59 2024: wf-6fcf1e16f3 []
# save_data -- Mon Dec  2 09:26:36 2024: wf-af996eacdf


ccm_api_save_data(term_md_filtered,
                  image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:SCAL_858_new_stage_term_meta_data_0.67.1.2")
# save_data -- Mon Dec  2 09:31:46 2024: wf-e501078032 []



# removing LUAD_tumor_normal and LUSC_tumor_normal meta because of the error that there is more than 1 meta with disease_vs_control 
term_md_filtered[!term_md_filtered$term_id %in% c('a3479c5','e203d77'),]


library(cytoreason.cc.client)
library(cytoreason.ccm.client)
devtools::load_all("~/analysis-p01-public-sphere/")

update_sphere_disease_model_entry(term_metadata = term_md_filtered,
                                  image = "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest",  #develop@0.57.14
                                  disease_model = "NSCLC-LU_p01", # take from R/ccm-disease-model.R
                                  tags = list(message = "export v39"),
                                  group_level_variables = '1_DCJGBvsF4YqBWipNhegjHg-DfxW99sQ2OVT0qTGXxY',
                                  sheet = 'p01_nsclc_v39')




