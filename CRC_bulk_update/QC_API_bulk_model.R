# install.packages("devtools") 
library(devtools)
# install.packages("memoise") # version 2.0.0 and up
library(memoise)
library(ggplot2)
library(magrittr)
library(checkmate) # for convenient assertions


#############################
# staging: 
# install.packages("devtools")
devtools::install_url("https://apps.private.cytoreason.com/sphere/client-tools/r-client/cytoreason.platform.client_1.1.1.tar.gz")
library(cytoreason.platform.client)
cr_init_save("pxx:sPJnT3eNTDd00uDWLSw60DchhO4Ds9KE:us")
Sys.setenv(PLATFORM_SERVER_URL="https://apps.private.cytoreason.com")
#############################



#############################
# production:
devtools::install_url("https://apps.cytoreason.com/sphere/client-tools/r-client/cytoreason.platform.client_1.1.2.tar.gz")
library(cytoreason.platform.client)
cr_init_save("p00:v2zb0LG2vEqJljrRpYjFjmZpExCeS2X2:us")
#############################

#############################
# pre-production:
devtools::install_url("https://apps.cytoreason.com/sphere/client-tools/r-client/cytoreason.platform.client_1.1.2.tar.gz")
cr_init_save("p00:v2zb0LG2vEqJljrRpYjFjmZpExCeS2X2:us")
Sys.setenv(PLATFORM_SERVER_URL="https://apps.pre.cytoreason.com")
#############################


avail_inventory <- cr_get_inventory()


# select the model of interest
selected_model <- subset(cr_get_inventory(), model_name == "colorectal cancer_colon (crb)")
# selected_model <- subset(cr_get_inventory(), model_name == "ulcerative colitis_colon")
# selected_model <- subset(cr_get_inventory(), model_name == "breast cancer_breast")

# cell abundance
cell_cont <- cr_model_get_cell_contribution(model_name=selected_model$model_name)

# cell gene expression
cell_gene_exp <- cr_model_get_cell_gene_expression(
  model_name = selected_model$model_name,
  dataset = "TCGA",
  genes = c('1','2','3','4','5','6'),
  input_validation = FALSE
)

cell_gene_exp <- cr_model_get_cell_gene_expression(
  model_name = selected_model$model_name,
  dataset = "meta",
  genes = c('1','2','3','4','5','6'),
  input_validation = FALSE
)


# gene cell correlation
inv <- cr_model_get_inventory(model_name = selected_model$model_name,
                              effect = "disease_vs_control:adjacent/tissue",
                              analysis_level = "dataset")
gene_cell_corr <- inv[1:3,] %>% cr_model_get_gene_cell_correlation(genes=c('1'))#,cells=c('dendritic cell', 'macrophage'))
str(gene_cell_corr)

eb577f9-308cc75


# meta gene cell correlation
inv <- cr_model_get_inventory(model_name = "Crohn's disease_ileum (crb)",
                              effect = "disease_vs_control",
                              analysis_level = "")
res <- inv %>% cr_model_get_gene_cell_correlation(genes=c('1','2'),
                                                  cells=c('dendritic cell', 'macrophage'))
str(res)





# retrieve the cell signatures
m_signature <- cr_bioentity_get_cell_signature(signature = "crc_v6")
# bioentity: cell annotations
View(cr_bioentity_cell_annotation())

 
tabl_inventory <- read.csv('~/sample_inventory.csv', check.names=FALSE)

# Term inventory: effect summary
# * list available effects
# * count number of meta-terms (groups of terms that are integrated together)
# * count number of dataset-level terms
# * count number of meta-analysis terms
# retrieve model inventory
m_inventory <- cr_model_get_inventory(model_name = selected_model$model_name) #####TODO: go over the number of samples in each group
m_inventory_to_csv <- m_inventory[, !colnames(m_inventory) %in% c("group_a_sample_list", "group_b_sample_list")] # list columns cannot be saved to csv
write.csv(m_inventory_to_csv, 'crc_62_inventory_API.csv')


# check if there are duplicated rows
colnames_for_dup_check <- c("dataset",
                            "effect", 
                            # "effect_variable", 
                            "effect_class", 
                            "effect_group_a", 
                            "effect_group_b", 
                            "term",
                            "term_condition",
                            "term_drug",
                            "term_drug_dosage")#,
                            # "group_a_samples",
                            # "group_b_samples")
m_inventory_dup_check <- m_inventory[, colnames_for_dup_check]
duplicate_rows <- m_inventory_dup_check[duplicated(m_inventory_dup_check), ]
m_inventory_subset <- m_inventory[,c("effect", "term", "effect_group_a", "effect_group_b", "dataset", "group_a_samples", "group_b_samples")]


# summarize the number of meta-terms, and count datasets
# and meta-analyses available in each effect
DT::datatable(m_inventory %>%
                dplyr::group_by(effect) %>%
                dplyr::summarise(
                  effect_term = length(unique(meta_term_id)),
                  dataset = length(unique(dataset[!is_meta])),
                  meta = sum(is_meta)
                )
              , options = list(
                scrollX = TRUE,
                scrollY = TRUE,
                scrollCollapse = TRUE
              ))

# meta gene cell correlation
inv <- cr_model_get_inventory(model_name = "Crohn's disease_ileum (crb)",
                              effect = "disease_vs_control",
                              analysis_level = "meta")
res <- inv %>% cr_model_get_gene_cell_correlation(genes=c('1','2'),
                                                  cells=c('dendritic cell 11111', 'macrophage'))#cells=c('dendritic cell', 'macrophage'))
str(res)


inv <- cr_model_get_inventory(model_name = selected_model$model_name,
                              effect = "disease_vs_control:adjacent",
                              analysis_level = "meta")
res <- inv %>% cr_model_get_gene_cell_correlation(genes=c('3501'))
str(res)


# Example 1: gene set differences
## Treatment effect
# retrieve gene set differences
# for specific treatment_this is a specific use case
term_df <- cr_model_get_inventory(model_name = "rheumatoid arthritis_synovium", 
                                  effect = c("post_vs_pre","disease_vs_control"),
                                  analysis_level = "all", #select meta or all
                                  effect %in% "disease_vs_control"|
                                    (term_drug %in% "Rituximab" & !is.na(group_a_response) & !is_meta)
)
#dim(term_df)
# select specific effects to see how many samples
term_df <- cr_model_get_inventory(model_name = selected_model$model_name, 
                                  effect = c("disease_vs_control:adjacent/tissue"),
                                  analysis_level = "all")

term_df <- cr_model_get_inventory(model_name = selected_model$model_name, 
                                  effect = c("disease_vs_control:adjacent"),
                                  analysis_level = "all")
# find table dimension
dim(term_df)
# annotate with term inventory metadata_pick submodel***
result_df <- cr_model_get_geneset_differences(term_df, 
                                              collection = "kegg", 
                                              submodel = "all") #%>%
  dplyr::inner_join(term_df) %>%
  dplyr::mutate(term_response_class = dplyr::recode(term_response,
                                                    "tumor" = "R",
                                                    "Non-responder" = "NR",
                                                    .missing = "DZ"))
# apply treatment effect logic_specific use case
res_logic <- result_df %>%
  dplyr::mutate("gene_set_stat" = ifelse(fdr <= 0.05, sign(statistic), 0)) %>%
  reshape2::dcast(feature_id + submodel ~ term_response_class, value.var = "gene_set_stat") %>%
  tidyr::unite("R_NR_DZ", "R", "NR", "DZ", sep = "_")
head(res_logic)
# get genes that are: up in disease, down in responders, not changed in Non-responders
GENESET <- subset(res_logic, R_NR_DZ %in% "-1_0_1" & submodel %in% "bulk")[, "feature_id"]
# plot
ggplot(result_df, aes(x = statistic, y = -log10(fdr))) +
  geom_point(colour = "#808080", size = 1) +
  # geom_point(data = subset(result_df, feature_id %in% GENESET),
  #            size = 1, aes(colour = as.character(sign(statistic)))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "darkred") +
  scale_color_brewer("Direction", palette = "Set1", direction = -1)+
  #facet_grid(submodel ~ effect + term_response) +
  facet_grid(dataset ~ term + submodel)+
  #facet_wrap(submodel ~ effect + term_response) +
  theme(
    legend.position = "top",
    text = element_text(size = 11)
  )


## Leading-edge symbols
# get all gene set information (not use here)
gs <- cr_bioentity_geneset_annotation()
# retrieve individual gene annotations
gene_annotations <- cr_bioentity_gene_annotation()
# add gene SYMBOL to gene set differences result
res_symb <- result_df %>%
  tidyr::unnest(leading_edge_list) %>%
  dplyr::left_join(gene_annotations[, c("feature_id", "SYMBOL")],
                   by = c(leading_edge_list = "feature_id")) %>%
  #collapse back to regular format
  dplyr::group_by(dplyr::across(c("resultset_id", "feature_id"))) %>%
  dplyr::mutate(SYMBOL = list(SYMBOL)) %>%
  dplyr::distinct(resultset_id, feature_id, .keep_all=TRUE)
res_symb %>% dplyr::select(feature_id,SYMBOL) %>% View()


# Example 2: cell composition differences
```{r ct_diff_plot, warning = FALSE, fig.align = "center", fig.height = 8, fig.width = 9}
# get cell composition differences
result_df <- cr_model_get_inventory(model_name = selected_model$model_name,
                                    analysis_level = "all",
                                    effect = 'disease_vs_control:adjacent/tissue') %>%
  cr_model_get_cell_differences()
# transform for plotting
plot_df <- dplyr::inner_join(result_df, cr_bioentity_cell_annotation()) %>%
  dplyr::mutate(feature_id = reorder(factor(feature_id), -log10(fdr)),
                Direction = ifelse(estimate > 0, "Up", "Down"))
# plot
ggplot(plot_df, aes(x = feature_id, y = -log10(fdr), fill = Direction)) +
  geom_bar(stat = "identity") +
  xlab("Cell type") +
  scale_fill_brewer(palette = "Set1", direction = -1, aes(colour = Direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "darkred") +
  coord_flip() +
  facet_grid(cell_category ~ ., switch = "y", scale = "free_y", space = "free") +
  theme(legend.position = "top")
```
# Example 3: gene differences
## Genes explained by cells
```{r gene diff}
# get gene bulk and adjusted differences
result_df <- cr_model_get_inventory(model_name = selected_model$model_name,
                                    analysis_level = "all",
                                    effect = "disease_vs_control:adjacent") %>%
  cr_model_get_gene_differences(submodel = "all")
# look for genes that are explained by cell composition: significant in bulk, not-significant in adjusted
GENESET <- intersect(
  subset(result_df, fdr <= 0.05 & submodel %in% "bulk")[, "feature_id"],
  subset(result_df, fdr > 0.05 & submodel %in% "adjusted")[, "feature_id"]
)
# plot
g=ggplot(result_df, aes(x = estimate, y = -log10(fdr))) +
  geom_point(colour = "#808080", size = 1) +
  geom_point(data = subset(result_df, feature_id %in% GENESET),
             size = 1, aes(colour = as.character(sign(estimate)))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "darkred") +
  scale_color_brewer("Direction", palette = "Set1", direction = -1)+
  facet_grid(submodel ~ term) +
  theme(
    legend.position = "top",
    text = element_text(size = 12)
  )
```
```{r gene diff}
result_df <- cr_model_get_inventory(model_name = "rheumatoid arthritis_synovium",
                                    analysis_level = "meta",
                                    effect = "disease_vs_control") %>%
  cr_model_get_gene_differences(submodel = "all")

#check specific genes of intrest ***
result_df <- cr_model_get_inventory(model_name = selected_model$model_name,
                                    effect = "disease_vs_control:adjacent/tissue") %>%
  # cr_model_get_gene_differences(submodel = "all")
cr_model_get_gene_cell_correlation(result_df[1:5,])
#same for treatments (my addition)
result_df <- cr_model_get_inventory(model_name = "rheumatoid arthritis_synovium",
                                    effect = "post_vs_pre" ,
                                    term %in% "W12_vs_W0:Partial_response") %>%
  cr_model_get_gene_differences(submodel = "all")
# draw a volcano plot
genes_of_intrest <- c() #no bio expectations
significance_cutoff <- 0.05
result_df$diffexpressed <- 'not significant'
result_df$diffexpressed[result_df$estimate > 0.5 & result_df$fdr <= 0.05] <- "Up regulated"
result_df$diffexpressed[result_df$estimate < -0.5 & result_df$fdr <= 0.05] <- "Down regulated"
g=ggplot(result_df, aes(x = estimate, y = -log10(fdr))) +
  xlab("Estimate") + ylab("-log10(FDR)") +
  ggtitle(paste0("Volcano Plot of ", result_df$term %>% unique(), " term, Bulk,", result_df$dataset %>% unique(), " \n \n Significantly down regulated: " ,
                 result_df %>% dplyr::filter(feature_label%in% genes_of_intrest & diffexpressed == 'Down regulated') %>% nrow(),'/',length(genes_of_intrest),
                 "\n Significantly up regulated: " ,
                 result_df %>% dplyr::filter(feature_label%in% genes_of_intrest & diffexpressed == 'Up regulated') %>% nrow(),'/',length(genes_of_intrest)
  )) +
  geom_point(aes(color=diffexpressed),alpha=0.5) +
  scale_color_manual(values=c("blue", "black", "red"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(aes(
    linetype = paste0("FDR = ", significance_cutoff),
    yintercept = -log10(significance_cutoff)
  ),
  color = "gray",
  size = 1) +
  scale_linetype_manual(values = c("dashed")) +
  guides(linetype = guide_legend(title = ""))
down=result_df %>% dplyr::filter(feature_label%in% genes_of_intrest & diffexpressed == 'Down regulated')
up=result_df %>% dplyr::filter(feature_label%in% genes_of_intrest & diffexpressed == 'Up regulated')
## Adjusted vs. bulk
```{r gene diff adj_vs_bulk}
# get gene bulk and adjusted differences
result_df <- cr_model_get_inventory(model_name = selected_model$model_name,
                                    analysis_level = "all", 
                                    effect = "disease_vs_control:adjacent") %>%
  cr_model_get_gene_differences(submodel = "all")
cutoff <- 0.05
gene_res_logic <- result_df %>%
  reshape2::dcast(feature_id ~ submodel, value.var = "fdr") %>%
  dplyr::mutate(
    significant = paste0(as.numeric(bulk <= cutoff), as.numeric(adjusted <= cutoff)),
    affected_genesets = dplyr::recode(significant,
                                      "00" = "Not significant",
                                      "10" = "Significant before adjustment",
                                      "01" = "Significant after adjustment",
                                      "11" = "Significant before and\nafter adjustment"
    )
  )
# Draw a plot showing gene set differences FDR, before and following adjustment for cell composition
ggplot(gene_res_logic,
       aes(x = -log10(bulk),
           y = -log10(adjusted),
           colour = affected_genesets)
) +
  geom_point() + xlab("-log10(FDR) Bulk") + ylab("-log10(FDR) Adjusted") +
  geom_hline(
    yintercept = -log10(cutoff),
    linetype = "dotted",
    col = "red"
  ) +
  geom_vline(
    xintercept = -log10(cutoff),
    linetype = "dotted",
    col = "red"
  ) + guides(colour = guide_legend(title = "Adjustment Effect"))
```