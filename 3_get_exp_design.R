###############################################################################
#
# Extract experimental design from Schraivogel et al (2020) data. 
#
###############################################################################

# minimal library imports to avoid namespace collisions
library(magrittr) # for pipes

# retrieve top-level data directory
schraivogel_dir <-.get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")

# names of experiments to extract experimental design for
experiment_names <- c("ground_truth_tapseq", "ground_truth_perturbseq")

# loop over experiments
for(exper_name in experiment_names){
  
  # create directory to store auxiliary data ###
  aux_dir <- sprintf("%sprocessed/%s/aux", schraivogel_dir, exper_name)
  if(!dir.exists(aux_dir)) dir.create(aux_dir, recursive = TRUE)
  
  # read in processed gRNA data
  processed_gRNA_dir <- sprintf("%sprocessed/%s/gRNA", 
                                schraivogel_dir, exper_name)
  gRNA_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
  gRNA_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
  gRNA_expr_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)
  
  # read in Supplementary Table 2
  supp_filename <- sprintf("%sraw/supp_tables/41592_2020_837_MOESM4_ESM.xlsx", 
                           schraivogel_dir)
  supp_table <- readxl::read_excel(supp_filename, sheet = "Chr 8 control library") 
  supp_table <- supp_table %>%
    dplyr::select(Gene, Target) %>%
    dplyr::rename(gene = Gene, target = Target) %>%
    dplyr::filter(target != "non-targeting") %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(gene = ifelse(target == "enhancer",
                                strsplit(gene, split = "-")[[1]][1],
                                gene)) %>%
    dplyr::ungroup() %>%
    unique()
  
  # merge information from gRNA ODM with the supplementary table
  feature_ids <- gRNA_expr_odm %>% ondisc::get_feature_ids()
  experimental_design <- tibble::tibble(gRNA_id = feature_ids) %>%
    dplyr::group_by(gRNA_id) %>%
    dplyr::mutate(gene = ifelse(grepl("non-targeting", gRNA_id), 
                                "non-targeting", 
                                strsplit(gRNA_id, split = "[-_]")[[1]][1])) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(supp_table, by = "gene") %>% 
    dplyr::mutate(target = ifelse(is.na(target), "non-targeting", target),
                  known_effect = gene,
                  gene = dplyr::case_when(
                    target == "promoter" ~ sprintf("%s-TSS", gene),
                    target == "enhancer" ~ sprintf("%s-enh", gene),
                    target == "non-targeting" ~ "non-targeting")) %>%
    dplyr::rename(target = gene, target_type = target)
  
  # save the result
  exp_design_filename = sprintf("%s/experimental_design.rds", aux_dir)
  saveRDS(experimental_design, exp_design_filename) 
}