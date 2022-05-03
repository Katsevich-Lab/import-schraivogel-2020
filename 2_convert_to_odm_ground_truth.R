###############################################################################
#
# Import Schraivogel et al (2020) data.
#
# Specific datasets imported:
#   - Ground truth TAP-seq experiment
#   - Ground truth perturb-seq experiment
#
# Notes:
#   - There is overlapping expression data in both ftp and geo directories. The
#     ftp directory contains gRNA expressions (in files Whole.nods.RDS and
#     TAD.nods.RDS) but the geo directory does not. Conversely, the geo directory
#     contains batch information but the ftp directory does not. Therefore,
#     information from both directories will need to be used.
#   - There were roughly 50 cell barcodes that appeared in more than one batch
#     of the perturb-seq experiment. For simplicity, I have excluded these cell
#     barcodes from the processed dataset.
#   - This dataset did not provide any ENSEMBL gene IDs, so I instead used gene
#     names instead of gene IDs.
#   - This script took roughly 10-15 minutes to run on a MacBook Pro.
###############################################################################

### retrieve top-level data directory ###
schraivogel_dir <-.get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")

### define the datasets to be imported ###
### - `name` is the name of the dataset as we will call it
### - `label` is the experiment name used to define filenames on GEO
### - `filename` is the name of the combined gene + gRNA expression file on FTP
experiments <- dplyr::tibble(name = c("ground_truth_tapseq",
                                     "ground_truth_perturbseq"),
                            label = c("TASC_DIFFEX",
                                      "WTX_DIFFEX"),
                            filename = c("TAP.nods.RDS",
                                         "Whole.nods.RDS"))

### loop over datasets to be imported ###
for(experiment in 1:nrow(experiments)){

  ### retrieve experiment names/labels ###
  exper_name <- experiments$name[experiment]
  exper_label <- experiments$label[experiment]
  exper_filename <- experiments$filename[experiment]
  cat(sprintf("Importing the %s dataset...\n", exper_name))

  ### create directories to store processed gene and gRNA data ###
  processed_gene_dir <- sprintf("%sprocessed/%s/gene",
                                schraivogel_dir, exper_name)
  processed_gRNA_dir <- sprintf("%sprocessed/%s/gRNA",
                                schraivogel_dir, exper_name)
  for(dir in c(processed_gene_dir, processed_gRNA_dir)){
    if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }

  #######################################################
  # 1. process batch information from GEO files
  #######################################################

  ### extract batch information from GEO files ###

  # get a list of all files in the geo directory
  geo_files <- list.files(sprintf("%sraw/geo", schraivogel_dir), full.names = TRUE)
  # determine number of batches based on the number of files with the correct
  # label (e.g. "TASC_DIFFEX") that end in "counts.csv"
  num_batches <- sum(grepl(exper_label, geo_files) & grepl("counts.csv", geo_files))
  # define a list with one element for each batch
  cells_to_batch <- vector("list", num_batches)
  # loop over batches
  for(batch_id in 1:num_batches){
    # the file for the given batch (e.g. 2) is determined based on the fact that
    # its filename contains a given substring (e.g. TASC_DIFFEX_sample2.counts.csv)
    substring <- sprintf("%s_sample%d.counts.csv", exper_label, batch_id)
    expression_file <- geo_files[grepl(substring, geo_files)]
    # extract the first line of this file, which consists of all the cell barcodes
    # this batch contains; split by ",", and remove the first entry because it is
    # blank
    cell_barcodes <- strsplit(readLines(expression_file, 1), split = ",")[[1]][-1]
    # create a tibble with first column storing cell barcode and second column
    # storing the batch number
    cells_to_batch[[batch_id]] <- dplyr::tibble(cell_barcode = cell_barcodes,
                                               batch = batch_id)
  }
  # concatenate all the batches together into one tibble
  cells_to_batch <- do.call("rbind", cells_to_batch) |>
    # prepend "batch_" to the batch column to make it a string
    dplyr::mutate(batch = paste0("batch_", batch))

  ### find cells without barcode collisions across batches ###
  cells_to_keep <- cells_to_batch |>
    dplyr::group_by(cell_barcode) |>
    dplyr::summarise(ncells = dplyr::n()) |>
    dplyr::filter(ncells == 1) |>
    dplyr::pull(cell_barcode)
  cat(sprintf("Keeping %d out of %d cells with unique barcodes.\n",
              length(cells_to_keep), nrow(cells_to_batch)))

  ### subset the cell-batch dictionary to cells without barcode collisions
  cells_to_batch <- cells_to_batch |>
    dplyr::filter(cell_barcode %in% cells_to_keep) |>
    # also change "."to "-" in the cell barcodes for compatibility with FTP files
    dplyr::mutate(cell_barcode = gsub("[.]", "-", cell_barcode))

  #########################################################################
  # 2. extract gene and gRNA expression information from FTP directory
  #########################################################################

  ### read raw expression matrix into memory and subset cells ###
  raw_expr_filename <- sprintf("%sraw/ftp/%s", schraivogel_dir, exper_filename)
  # Note: the matrix raw_expr_data below contains both the gene and gRNA
  # expression data, with gRNA rownames containing the prefix "CROPseq_dCas9_DS_"
  raw_expr_data <- readRDS(raw_expr_filename)
  raw_expr_data <- raw_expr_data[,cells_to_batch$cell_barcode]

  ### convert raw expression matrix to sparse dgRMatrix format
  if(!("dgRMatrix" %in% class(raw_expr_data))){
    cat("Converting to sparse dgRMatrix format for memory efficiency and ")
    cat("compatibility with ondisc...\n")
    raw_expr_data <- as.matrix(raw_expr_data)
    raw_expr_data <- as(raw_expr_data, "dgRMatrix")
  }

  ### extract gRNA expression data from the combined expression matrix ###
  gRNA_prefix <- "CROPseq_dCas9_DS_"
  gRNA_rows <- grepl(gRNA_prefix, rownames(raw_expr_data))
  gRNA_expr_data <- raw_expr_data[gRNA_rows,]
  # remove the prefixes from the gRNA expression matrix
  remove_prefix <- function(rowname){
    strsplit(rowname, gRNA_prefix)[[1]][2]
  }
  rownames(gRNA_expr_data) <- unname(sapply(rownames(gRNA_expr_data),
                                            remove_prefix))

  ### extract gene expression data from the combined expression matrix ###
  gene_rows <- !gRNA_rows
  gene_expr_data <- raw_expr_data[gene_rows,]
  rm(raw_expr_data)

  #######################################################
  # 3. process gRNA metadata
  #######################################################

  # read in Supplementary Table 2
  supp_filename <- sprintf("%sraw/supp_tables/41592_2020_837_MOESM4_ESM.xlsx",
                           schraivogel_dir)
  supp_table <- readxl::read_excel(supp_filename, sheet = "Chr 8 control library")
  supp_table <- supp_table |>
    dplyr::select(Gene, Target) |>
    dplyr::rename(gene = Gene, target = Target) |>
    dplyr::filter(target != "non-targeting") |>
    dplyr::group_by(gene) |>
    dplyr::mutate(gene = ifelse(target == "enhancer",
                                strsplit(gene, split = "-")[[1]][1],
                                gene)) |>
    dplyr::ungroup() |>
    unique()

  # join supplementary table with list of gRNAs
  experimental_design <- tibble::tibble(gRNA_id = rownames(gRNA_expr_data)) |>
    dplyr::group_by(gRNA_id) |>
    dplyr::mutate(gene = ifelse(grepl("non-targeting", gRNA_id),
                                "non-targeting",
                                strsplit(gRNA_id, split = "[-_]")[[1]][1])) |>
    dplyr::ungroup() |>
    dplyr::left_join(supp_table, by = "gene") |>
    dplyr::mutate(target = ifelse(is.na(target), "non-targeting", target),
                  known_effect = gene,
                  gene = dplyr::case_when(
                    target == "promoter" ~ sprintf("%s-TSS", gene),
                    target == "enhancer" ~ sprintf("%s-enh", gene),
                    target == "non-targeting" ~ "non-targeting")) |>
    dplyr::rename(target = gene, target_type = target)

  #########################################################################
  # 3. create ODM for gene expression matrix
  #########################################################################
  cat("Creating ODM for gene expression matrix...\n")
  cell_barcodes <- colnames(gene_expr_data)
  gene_names <- dplyr::tibble(gene_name = rownames(gene_expr_data))
  odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
  metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
  ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gene_expr_data,
                                             barcodes = cell_barcodes,
                                             features_df = gene_names,
                                             odm_fp = odm_fp,
                                             metadata_fp = metadata_fp) |>
    # add batch information to cell covariates
    ondisc::mutate_cell_covariates(batch = cells_to_batch$batch) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp)

  #########################################################################
  # 4. create ODM for gRNA expression matrix
  #########################################################################
  cat("Creating ODM for gRNA expression matrix...\n")
  cell_barcodes <- colnames(gRNA_expr_data)
  gRNA_names <- dplyr::tibble(gRNA_name = rownames(gRNA_expr_data))
  odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
  metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
  ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_expr_data,
                                             barcodes = cell_barcodes,
                                             features_df = gRNA_names,
                                             odm_fp = odm_fp,
                                             metadata_fp = metadata_fp) |>
    # add batch information to cell covariates
    ondisc::mutate_cell_covariates(batch = cells_to_batch$batch) |>
    # add gRNA metadata to feature covariates
    ondisc::mutate_feature_covariates(target = experimental_design$target,
                                      target_type = experimental_design$target_type,
                                      known_effect = experimental_design$known_effect) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp)

  cat("Done.\n")
}
