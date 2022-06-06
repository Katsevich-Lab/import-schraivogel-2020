schraivogel_dir <-.get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")

experiments <- list.dirs(paste0(schraivogel_dir, "processed"),
                         full.names = FALSE,
                         recursive = FALSE)

for(experiment in experiments){
  cat(sprintf("Adding gRNA assignments to experiment %s...\n", experiment))
  
  # read gene odm
  gene_odm_fp <- sprintf("%sprocessed/%s/gene/expression_matrix.odm", schraivogel_dir, experiment)
  gene_metadata_fp <- sprintf("%sprocessed/%s/gene/metadata.rds", schraivogel_dir, experiment)
  gene_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)
  
  # read grna odm
  grna_odm_fp <- sprintf("%sprocessed/%s/grna/raw_ungrouped.odm", schraivogel_dir, experiment)
  grna_metadata_fp <- sprintf("%sprocessed/%s/grna/raw_ungrouped_metadata.rds", schraivogel_dir, experiment)
  grna_odm <- ondisc::read_odm(grna_odm_fp, grna_metadata_fp)

  # find perturbation assignments  
  grna_mat <- grna_odm[[1:nrow(grna_odm), 1:ncol(grna_odm)]]
  grna_ids <- grna_odm |> ondisc::get_feature_ids()
  perturbations <- apply(grna_mat, 2, function(col)(grna_ids[col >= 8]))

  # add to response odm and save
  gene_odm |> 
    ondisc::mutate_cell_covariates(perturbation = perturbations) |>
    ondisc::save_odm(metadata_fp = gene_metadata_fp)
}