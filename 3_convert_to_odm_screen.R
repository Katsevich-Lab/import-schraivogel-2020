###############################################################################
#
# Import Schraivogel et al (2020) data.
#
# Specific datasets imported:
#   - Screen data from chromosome 8
#   - Screen data from chromosome 11
#
# Notes:
#   - This script took roughly 8 minutes to run on a MacBook Pro.
###############################################################################

### retrieve top-level data directory ###
schraivogel_dir <-.get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")

### define the datasets to be imported ###
### - `name` is the name of the dataset as we will call it
### - `label` is the experiment name used to define filenames on GEO
experiments <- dplyr::tibble(name = c("enhancer_screen_chr8",
                                      "enhancer_screen_chr11"),
                             label = c("TASC_SCREEN_chr8",
                                       "TASC_SCREEN_chr11"))


for(experiment in 1:nrow(experiments)){
  
  ### retrieve experiment names/labels ###
  exper_name <- experiments$name[experiment]
  exper_label <- experiments$label[experiment]
  cat(sprintf("Importing the %s dataset...\n", exper_name))

  ### create directories to store processed gene and gRNA data ###
  processed_gene_dir <- sprintf("%sprocessed/%s/gene",
                                schraivogel_dir, exper_name)
  processed_gRNA_expression_dir <- sprintf("%sprocessed/%s/grna_expression",
                                           schraivogel_dir, exper_name)
  processed_gRNA_assignment_dir <- sprintf("%sprocessed/%s/grna_assignment",
                                           schraivogel_dir, exper_name)
  for(dir in c(processed_gene_dir, 
               processed_gRNA_expression_dir,
               processed_gRNA_assignment_dir)){
    if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  #######################################################
  # 1. Read in the expression matrices
  #######################################################
  
  cat("Reading in the expression matrices...\n")
  geo_files <- list.files(sprintf("%sraw/geo", schraivogel_dir), full.names = TRUE)
  geo_expr_files <- geo_files[grepl("counts.csv", geo_files) & grepl(exper_label, geo_files)]
  
  read_file <- function(filename){
    batch <- strsplit(filename, split = ".counts.csv|_")[[1]] |> dplyr::last()
    readr::read_csv(filename) |>
      dplyr::rename_with(function(col_names)(c("feature_name", 
                                               sprintf("%s_%s", col_names[-1], batch))))
  }
  expr_data_list <- lapply(geo_expr_files, read_file)

  # make sure all expression datasets have the same feature names
  feature_names <- expr_data_list[[1]]$feature_name
  same_feature_names <- lapply(expr_data_list, 
                               function(expr_data)(all(expr_data$feature_name == feature_names))) |> 
    unlist() |> 
    all()
  if(!same_feature_names) stop("All data frames should have the same feature names!")
  
  # join all the expression datasets by horizontal concatenation
  all_expr_data <- lapply(expr_data_list, 
                   function(expr_data)(expr_data |> 
                                         dplyr::select(-feature_name)
                                         )) |> 
    dplyr::bind_cols() |>
    as.matrix()
  
  ### extract gene and gRNA expression data from the combined expression matrix ###
  gRNA_prefix <- "CROPseq_dCas9_DS_"
  gRNA_rows <- grepl(gRNA_prefix, feature_names)
  gene_expr_data <- all_expr_data[!gRNA_rows,]
  gRNA_expr_data <- all_expr_data[gRNA_rows,]
  rm(all_expr_data)

  #######################################################
  # 2. Create ODM for gene expression matrix
  #######################################################
  
  cat("Creating ODM for gene expression matrix...\n")
  cell_barcodes <- colnames(gene_expr_data)
  gene_names <- dplyr::tibble(gene_name = feature_names[!gRNA_rows])
  batches <- gene_expr_data |> 
    colnames() |> 
    sapply(function(col_name)(strsplit(col_name, split = "_")[[1]][2])) |> 
    unname()
  odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
  metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
  ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gene_expr_data,
                                             barcodes = cell_barcodes,
                                             features_df = gene_names,
                                             odm_fp = odm_fp,
                                             metadata_fp = metadata_fp) |>
    # add batch information to cell covariates
    ondisc::mutate_cell_covariates(batch = batches) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp)
  rm(gene_expr_data)
  
  
  #######################################################
  # 3. process gRNA metadata
  #######################################################

  cat("Processing gRNA metadata...\n")
  
  # read in Supplementary Table 2, choosing the relevant sheet of the Excel file
  supp_filename <- sprintf("%sraw/supp_tables/41592_2020_837_MOESM4_ESM.xlsx",
                           schraivogel_dir)
  if(exper_name == "enhancer_screen_chr8"){
    supp_table_ctrl <- readxl::read_excel(supp_filename, sheet = "Chr 8 control library") 
  } else if(exper_name == "enhancer_screen_chr11"){
    supp_table_ctrl <- readxl::read_excel(supp_filename, sheet = "Chr 11 control library") 
  } else{
    stop("Wrong experiment name!")
  }
  
  # massage the supplementary table
  supp_table_ctrl <- supp_table_ctrl |>
    dplyr::select(Gene, Target) |>
    dplyr::rename(gene = Gene, target_type = Target) |>
    dplyr::filter(target_type != "non-targeting") |>
    dplyr::group_by(gene) |>
    dplyr::mutate(gene = ifelse(target_type == "enhancer",
                                strsplit(gene, split = "-")[[1]][1],
                                gene)) |>
    dplyr::ungroup() |>
    unique() |>
    dplyr::mutate(known_effect = gene,
                  target = dplyr::case_when(
                    target_type == "promoter" ~ sprintf("%s-TSS", gene),
                    target_type == "enhancer" ~ sprintf("%s-enh", gene))) |>
    dplyr::select(-gene)
  
  # remove the prefixes from the gRNA expression matrix
  remove_prefix <- function(rowname){
    strsplit(rowname, gRNA_prefix)[[1]][2]
  }
  gRNA_names <- unname(sapply(feature_names[gRNA_rows], remove_prefix))
  
  # join supplementary table with list of gRNAs
  experimental_design <- tibble::tibble(gRNA_id = gRNA_names) |>
    dplyr::rowwise() |>
    dplyr::mutate(known_effect = dplyr::case_when(
      grepl("non-targeting", gRNA_id) ~ "non-targeting",
      grepl("chr(8|11):", gRNA_id) ~ NA_character_,
      TRUE ~ strsplit(gRNA_id, split = "[-_]")[[1]][1])) |>
    dplyr::left_join(supp_table_ctrl, by = "known_effect") |>
    dplyr::mutate(target = dplyr::case_when(
      is.na(known_effect) ~ strsplit(gRNA_id, split = "_")[[1]][1],
      known_effect == "non-targeting" ~ "non-targeting",
      TRUE ~ target),
      target_type = dplyr::case_when(
        is.na(known_effect) ~ "enhancer",
        known_effect == "non-targeting" ~ "non-targeting",
        TRUE ~ target_type)
      ) |>
      dplyr::ungroup() |>
    dplyr::select(gRNA_id, target, target_type, known_effect)

  #########################################################################
  # 4. create ODM for gRNA expression matrix
  #########################################################################
  cat("Creating ODM for gRNA expression matrix...\n")
  cell_barcodes <- colnames(gRNA_expr_data)
  gRNA_names <- dplyr::tibble(gRNA_name = gRNA_names)
  odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_expression_dir)
  metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_expression_dir)
  grna_expression_odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_expr_data,
                                             barcodes = cell_barcodes,
                                             features_df = gRNA_names,
                                             odm_fp = odm_fp,
                                             metadata_fp = metadata_fp) |>
    # add batch information to cell covariates
    ondisc::mutate_cell_covariates(batch = batches) |>
    # add gRNA metadata to feature covariates
    ondisc::mutate_feature_covariates(target = experimental_design$target,
                                      target_type = experimental_design$target_type,
                                      known_effect = experimental_design$known_effect)
  # save to disk
  grna_expression_odm |>
    ondisc::save_odm(metadata_fp = metadata_fp)
  rm(gRNA_expr_data)
  
  
  #########################################################################
  # 5. create ODM for gRNA assignment matrix
  #########################################################################
  cat("Thresholding gRNA expression matrix...\n")
  grna_mat <- grna_expression_odm[[1:nrow(grna_expression_odm), 1:ncol(grna_expression_odm)]]
  grna_ids <- gRNA_names$gRNA_name
  gRNA_assignment_list <- apply(grna_mat, 2, function(col)(grna_ids[col >= 8])) |>
    lapply(function(perts)(ifelse(length(perts) == 0, "", perts)))
  
  cat("Creating ODM for gRNA assignment matrix...\n")
  odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_assignment_dir)
  metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_assignment_dir)
  ondisc::convert_assign_list_to_sparse_odm(
    cell_barcodes = cell_barcodes,
    gRNA_ids = grna_ids,
    gRNA_assignment_list = gRNA_assignment_list,
    odm_fp = odm_fp,
    metadata_fp = metadata_fp
  ) |>
    # add gRNA metadata to feature covariates
    ondisc::mutate_feature_covariates(
      target = experimental_design$target,
      target_type = experimental_design$target_type,
      known_effect = experimental_design$known_effect
    ) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp)
  
  cat("Done.\n")
}
