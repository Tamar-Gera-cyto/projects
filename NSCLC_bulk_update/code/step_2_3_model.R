####### Step 2 - Generate text output #######
# install.packages('cytoreason.cc.client')
library(cytoreason.cc.client)
# install.packages('cytoreason.ccm.pipeline')
library(cytoreason.ccm.pipeline)

disease_name <- 'non-small cell lung carcinoma'
ccm_wf_id <-   'wf-8d8aeeee6c' # v43 CRS with image SCAL_858_new_stage_term_meta_data_0.65.3.1
              # 'wf-c58e381360' # with TMW image from Inbal
              # 'wf-aad38df64b' # v41 - removing for B9991027: LUSC_TMB, TMB_high_LUAD_LUSC and adding back LUSC_lightly_non_responder_post_vs_pre, LUAD_pre_non_responder_heavily_vs_lightly
              #'wf-5ab3b37d57' # with TME
              #"wf-089505682f" # v39 - adding smoking and PD1_or_PDL1_post_vs_pre for TEMPUS 
              # "wf-e81aed1f43" # v38 - correcting according to Lital's comments
              #"wf-7fed6b65c8" # v37 - with model subset, removung "per_subtype" comparisons and adding instead 2 rows of comparisons 
              # "wf-d7230cfbef" # v36 (same as v35 but with CRS)
              #"wf-a3f8f4a776" # v35 in tableau (same as v34 with prior v12)
              #"wf-f2ecc87b56" # v34 in tableau (combined model after adding the paired post_vs_pre (including metastasis samples) comparisons)
              #"wf-948b7d7d8d" # v33 in tableau (combined model after adding the general LUAD+LUSC post_pre, R/NR, heavily/lightly and LUAD/LUSC)
              #"wf-ed70cf5444" # V32 in tableau (combined old model+TEMPUS with corrections) prior v11
              #"wf-e378cc4755" # V31 in tableau (combined old model+TEMPUS) prior v10
              #"wf-00bf4f535a" # V25 in tableau (old version+corrections prior v6)
              # "wf-e4e1111047" # V24 in tableau (TEMPUS with 3 LUAD_LUSC comparisons)
              # "wf-21b281b964" # V23 in Tableau
              # "wf-0662ff9fab" # V22 in Tableau
              # "wf-362c9d2aa8" # V21 in Tableau
              # "wf-9585c946a5" # V20 in Tableau
memory_request <- '32Gi'
memory_request_downstream <- '[20Gi]'
wf <- ccm_api_save_data(AssetData(ccm_wf_id),
                        DirectoryTXT(), #bio_qc = F,
                        image = ccm_cyto_cc("eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest",
                          # "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:MT_4002_run_model_with_TME_classifier_0.67.1.2", #image that generates
                          # "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:master_latest", # used to be develop, changed on 29.10.24
                          # "eu.gcr.io/cytoreason/ci-cytoreason.ccm.pipeline-package:develop_latest",
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

# wf-9585c946a5: save_data -- Thu Aug 29 12:57:27 2024: wf-f09e280a82 - V20: Naama's API version with Avital new lung prior test_lung_v3:new_io_nsclc_lung_CRB_sig (no TEMPUS yet)
# wf-362c9d2aa8: save_data -- Fri Aug 30 17:05:31 2024: wf-3ff2381f0b - V21: Adding first comparisons for TEMPUS, same prior as V20
# wf-0662ff9fab: save_data -- Tue Sep  3 11:58:16 2024: wf-bcccac243f - V22: running on TEMPUS lung samples, adding smoking, gender and treatment_intensity_response
# wf-21b281b964: save_data -- Thu Sep 19 17:57:06 2024: wf-942cd11663 - V23: after my corrections to the old version, without TEMPUS
# wf-e4e1111047: save_data -- Thu Sep 26 11:34:51 2024: wf-5596293b96 - V24: only TEMPUS (with 3 DS for LUAD/LUSC comparison for meta PCA)
# wf-00bf4f535a: save_data -- Tue Oct  1 05:57:48 2024: wf-85b36ca8de - V25: old version+corrections prior v6
# wf-e378cc4755: save_data -- Mon Oct 21 05:25:34 2024: wf-816af3cb64 - V31: combined old model+TEMPUS prior v10
# wf-ed70cf5444: save_data -- Wed Oct 30 07:07:54 2024: wf-b2844987c0 - v32: combined old model+TEMPUS with corrections prior v11
# wf-948b7d7d8d: save_data -- Thu Oct 31 16:49:46 2024: wf-ed8d3bccea - v33: combined model after adding the general LUAD+LUSC post_pre, R/NR, heavily/lightly and LUAD/LUSC 
# wf-f2ecc87b56: save_data -- Fri Nov  1 17:51:37 2024: wf-7997c80af9 - v34: combined model after adding the paired post_vs_pre (including metastasis samples) comparisons
# wf-a3f8f4a776: save_data -- Mon Nov  4 07:54:57 2024: wf-b56bf47703 - v35:same as v34 with prior v12
# wf-d7230cfbef: save_data -- Wed Nov  6 17:37:06 2024: wf-f6d3f63b49 - v36: (same as v35 but with CRS)
# wf-7fed6b65c8: save_data -- Fri Nov  8 17:33:52 2024: wf-a73ccd1968 - v37: with model subset, removung "per_subtype" comparisons and adding instead 2 rows of comparisons 
# wf-e81aed1f43: save_data -- Mon Nov 18 07:45:26 2024: wf-e1330af1a4 - v38: correcting according to Lital's comments
# wf-089505682f: save_data -- Mon Nov 18 17:07:57 2024: wf-cd51390a6e - v39: adding smoking and PD1_or_PDL1_post_vs_pre for TEMPUS 
# wf-5ab3b37d57: save_data -- Mon Nov 25 08:32:28 2024: wf-f104d737f7 - with TME - failed
# wf-5ab3b37d57: save_data -- Mon Nov 25 17:38:55 2024: wf-747f38a1fc - with TME - failed
# wf-aad38df64b: save_data -- Sun Dec  1 08:49:01 2024: wf-307bd353fb - v41 - removing for B9991027: LUSC_TMB, TMB_high_LUAD_LUSC and adding back LUSC_lightly_non_responder_post_vs_pre, LUAD_pre_non_responder_heavily_vs_lightly
# wf-c58e381360: save_data -- Tue Dec 10 08:37:34 2024: wf-fb7ddf18da - with TME image from Inbal
# wf-c58e381360: save_data -- Wed Dec 11 10:34:16 2024: wf-2fdb05eea3 - with TME changing the image in the uploader to the same as the pipeline (Inbal's image)
# wf-8d8aeeee6c: save_data -- Thu Dec 19 09:12:54 2024: wf-3578263d70 - v43 CRS (CCM ran with image SCAL_858_new_stage_term_meta_data_0.65.3.1), step 2 with master latest

####### Step 3 - upload text output #######

library(cytoreason.cc.client)

#########
# image #
#########
# task_image = 'eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest'
task_image = 'eu.gcr.io/cytoreason/cd-py-bigquery:master_latest' # used to be develop, changed on 29.10.24
memory_request = '64Mi'

#####################
# ARGS (positional) #
#####################
service = 'ccm'
dataset = "p01_io_nsclc" ## !!!!!!!! Pay attention to the DB name !!!!!!!! #
wf_id = "wf-3578263d70"


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
                       outdir = "/outdir/",
                       image = task_image ,
                       task_env_vars = task_env_vars,
                       tags=tags,
                       memory_request = memory_request)
wf
# Cyto-CC workflow: wf-37f2a9c40a
# Cyto-CC workflow: wf-77e6d07018
# Cyto-CC workflow: wf-aacdb5f5bf
# Cyto-CC workflow: wf-fc2ca57440
# Cyto-CC workflow: wf-c160049f14
# Cyto-CC workflow: wf-d8d11346cd
# Cyto-CC workflow: wf-bb54aeea7e
# Cyto-CC workflow: wf-0d466849bd
# Cyto-CC workflow: wf-ca698f73e8 - V33 failed
# Cyto-CC workflow: wf-c087f31244 - V34 
# Cyto-CC workflow: wf-339a4cd372
# Cyto-CC workflow: wf-339a4cd372 - v33 again (ccm wf: wf-948b7d7d8d, txt dump: wf-ed8d3bccea) 
# Cyto-CC workflow: wf-365f0cc220 - v34 again (ccm wf: wf-f2ecc87b56, txt dump: wf-7997c80af9) 
# Cyto-CC workflow: wf-8c625958fe - v35
# Cyto-CC workflow: wf-ede7e0067c - v36
# Cyto-CC workflow: wf-10206be456 - v37
# Cyto-CC workflow: wf-9f1e96daf1 - v38
# Cyto-CC workflow: wf-4b76fe8518 - v39
# Cyto-CC workflow: wf-9c96e75c1a - v41
# Cyto-CC workflow: wf-e3ea0185d1 - v43


