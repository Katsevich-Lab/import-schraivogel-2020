Schraivogel (2020) Data Documentation
================
Gene Katsevich;
June 8, 2022

# Overview

This repository contains code to import and process the Schraivogel 2020
data. Schraivogel et al developed the targeted perturb-seq (TAP-seq)
protocol, a new assay using direct capture both for gRNAs and genes. The
point of sequencing only targeted genes is to reduce sequencing costs by
focusing attention only on relevant genes. Schraivogel et al carried out
proof-of-concept “ground truth” experiments based on TAP-seq and
perturb-seq, as well as an at-scale experiment to find enhancers for 147
total genes across chromosomes 8 and 11.

The `schraivogel-2020` directory structure is as follows:

    ├── processed
    │   ├── enhancer_screen_chr11
    │   │   ├── gene
    │   │   │   ├── expression_matrix.odm
    │   │   │   └── metadata.rds
    │   │   ├── grna_assignment
    │   │   │   ├── raw_ungrouped.odm
    │   │   │   └── raw_ungrouped_metadata.rds
    │   │   └── grna_expression
    │   │       ├── raw_ungrouped.odm
    │   │       └── raw_ungrouped_metadata.rds
    │   ├── enhancer_screen_chr8
    │   │   ├── gene
    │   │   │   ├── expression_matrix.odm
    │   │   │   └── metadata.rds
    │   │   ├── grna_assignment
    │   │   │   ├── raw_ungrouped.odm
    │   │   │   └── raw_ungrouped_metadata.rds
    │   │   └── grna_expression
    │   │       ├── raw_ungrouped.odm
    │   │       └── raw_ungrouped_metadata.rds
    │   ├── ground_truth_perturbseq
    │   │   ├── gene
    │   │   │   ├── expression_matrix.odm
    │   │   │   └── metadata.rds
    │   │   ├── grna_assignment
    │   │   │   ├── raw_ungrouped.odm
    │   │   │   └── raw_ungrouped_metadata.rds
    │   │   └── grna_expression
    │   │       ├── raw_ungrouped.odm
    │   │       └── raw_ungrouped_metadata.rds
    │   └── ground_truth_tapseq
    │       ├── gene
    │       │   ├── expression_matrix.odm
    │       │   └── metadata.rds
    │       ├── grna_assignment
    │       │   ├── raw_ungrouped.odm
    │       │   └── raw_ungrouped_metadata.rds
    │       └── grna_expression
    │           ├── raw_ungrouped.odm
    │           └── raw_ungrouped_metadata.rds
    └── raw
        ├── ...

The contents of the `raw` directory are suppressed, as they are
unimportant. The `processed` directory contains four subdirectories
(`ground_truth_tapseq`, `ground_truth_perturbseq`,
`enhancer_screen_chr8`, `enhancer_screen_chr11`), which correspond to
distinct datasets. The first two are for the ground truth experiment,
and the last two are for the at scale enhancer screen experiment. All
four have the same subdirectory structure, containing three modalities:
`gene`, `grna_expression`, and `grna_assignment`. The first two contain
gene and gRNA expression ODMs, while the last contains a logical ODM
corresponding to Schraivogel’s gRNA assignments to cells.

# Ground truth experiments

There are two separate ground truth experiments: one TAP-seq and one
perturb-seq. These experiments are meant as a proof-of-concept for
TAP-seq, so they contain only positive and negative control
perturbations (perturbations for which the ground truth is known). Both
experiments have the same experimental design; they differ only in that
TAP-seq targets a small number of genes while perturb-seq targets the
whole transcriptome.

## Experimental design

As mentioned above, the experimental design for these two ground truth
experiments is the same. The information about the gRNAs is contained in
the last three columns of the feature covariate matrix of both
`gRNA_expression` and `gRNA_assignment` `ondisc` matrices:

``` r
# load the gRNA expression data
processed_dir <- sprintf("%s/processed", 
                         .get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR"))
processed_gRNA_dir <- sprintf("%s/ground_truth_tapseq/gRNA_expression", processed_dir)
gRNA_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
gRNA_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_expr_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)

# extract the feature covariates containing the gRNA metadata
exper_design <- gRNA_expr_odm |> 
  ondisc::get_feature_covariates() |>
  dplyr::select(target, target_type, known_effect) |>
  tibble::rownames_to_column(var = "gRNA") 

# print the first five rows
exper_design |>
  head(5) |> 
  kableExtra::kable(format = "html",
                    booktabs = TRUE, 
                    col.names = c("gRNA", "Target", "Target Type", "Known Effect"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
gRNA
</th>
<th style="text-align:left;">
Target
</th>
<th style="text-align:left;">
Target Type
</th>
<th style="text-align:left;">
Known Effect
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
CCNE2\_+\_95907328.23-P1P2
</td>
<td style="text-align:left;">
CCNE2-TSS
</td>
<td style="text-align:left;">
promoter
</td>
<td style="text-align:left;">
CCNE2
</td>
</tr>
<tr>
<td style="text-align:left;">
CCNE2\_+\_95907382.23-P1P2
</td>
<td style="text-align:left;">
CCNE2-TSS
</td>
<td style="text-align:left;">
promoter
</td>
<td style="text-align:left;">
CCNE2
</td>
</tr>
<tr>
<td style="text-align:left;">
CCNE2\_+\_95907406.23-P1P2
</td>
<td style="text-align:left;">
CCNE2-TSS
</td>
<td style="text-align:left;">
promoter
</td>
<td style="text-align:left;">
CCNE2
</td>
</tr>
<tr>
<td style="text-align:left;">
CCNE2\_-\_95907017.23-P1P2
</td>
<td style="text-align:left;">
CCNE2-TSS
</td>
<td style="text-align:left;">
promoter
</td>
<td style="text-align:left;">
CCNE2
</td>
</tr>
<tr>
<td style="text-align:left;">
CPQ\_+\_97657557.23-P1P2
</td>
<td style="text-align:left;">
CPQ-TSS
</td>
<td style="text-align:left;">
promoter
</td>
<td style="text-align:left;">
CPQ
</td>
</tr>
</tbody>
</table>

The full table contains a total of 86 rows, so there are 86 gRNAs in
this experiment. Breaking these down by their target,

``` r
exper_design |> pull(target) |> table()
```

    ## 
    ##     CCNE2-TSS       CPQ-TSS     DSCC1-TSS    FAM83A-TSS     GATA1-enh 
    ##             4             4             4             4             4 
    ##       HS2-enh    LRRCC1-TSS       MYC-enh non-targeting      OXR1-TSS 
    ##             4             4             4            30             4 
    ##   PHF20L1-TSS     RIPK2-TSS      STK3-TSS      UBR5-TSS     ZFPM2-enh 
    ##             4             4             4             4             4

we see that we have 30 non-targeting gRNAs as well as four gRNAs each
for 14 targets
(![30 + 4\\times14 = 86](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;30%20%2B%204%5Ctimes14%20%3D%2086 "30 + 4\times14 = 86")).
Of these targets, 10 are gene TSSs and 4 are well-characterized
enhancers of four genes.

## TAP-seq data details

Let’s now take a look at the TAP-seq data:

``` r
# load the gRNA expression data
processed_gRNA_expression_dir <- sprintf("%s/ground_truth_tapseq/gRNA_expression", processed_dir)
gRNA_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_expression_dir)
gRNA_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_expr_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)
gRNA_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 86 features and 21977 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero, target, target_type, known_effect.

``` r
# load the gRNA assignment data
processed_gRNA_assignment_dir <- sprintf("%s/ground_truth_tapseq/gRNA_assignment", processed_dir)
gRNA_assignment_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_assignment_dir)
gRNA_assignment_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_assignment_dir)
gRNA_assign_odm <- ondisc::read_odm(gRNA_assignment_odm_fp, gRNA_assignment_metadata_fp)
gRNA_assign_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 86 features and 21977 cells.
    ##  A cell covariate matrix with columns n_nonzero.
    ##  A feature covariate matrix with columns n_nonzero, target, target_type, known_effect.

``` r
# load the gene expression data
processed_gene_dir <- sprintf("%s/ground_truth_tapseq/gene", processed_dir)
gene_odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
gene_metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
gene_expr_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)
gene_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 72 features and 21977 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

This experiment has 21977 cells across 2 batches. The gRNA data come in
the form of expressions and are not thresholded. There are a total of 86
gRNAs, as discussed above. A total of 72 genes are measured. Based on
the paper, there are supposed to be 74 genes measured: 14 that were
targeted and 60 presumably unrelated genes. Of the two missing genes,
one is HS2 (whose enhancer was targeted) and one is a presumably
unrelated gene. We can look further into why this is the case, but
perhaps it’s not urgent.

## Perturb-seq data details

Next, we turn to the perturb-seq data:

``` r
# load the gRNA expression data
processed_gRNA_dir <- sprintf("%s/ground_truth_perturbseq/gRNA_expression", processed_dir)
gRNA_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
gRNA_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_expr_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)
gRNA_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 85 features and 37918 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero, target, target_type, known_effect.

``` r
# load the gRNA assignment data
processed_gRNA_assignment_dir <- sprintf("%s/ground_truth_perturbseq/gRNA_assignment", processed_dir)
gRNA_assignment_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_assignment_dir)
gRNA_assignment_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_assignment_dir)
gRNA_assign_odm <- ondisc::read_odm(gRNA_assignment_odm_fp, gRNA_assignment_metadata_fp)
gRNA_assign_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 85 features and 37918 cells.
    ##  A cell covariate matrix with columns n_nonzero.
    ##  A feature covariate matrix with columns n_nonzero, target, target_type, known_effect.

``` r
# load the gene expression data
processed_gene_dir <- sprintf("%s/ground_truth_perturbseq/gene", processed_dir)
gene_odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
gene_metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
gene_expr_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)
gene_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 17107 features and 37918 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

This experiment has 37918 cells across 4 batches. The gRNA data come in
the form of expressions and are not thresholded. There are a total of 85
gRNAs, which is one fewer than the 86 in the experimental design.
Perhaps the missing one (STK3\_-\_99837866.23-P1P2) got removed during
QC by Schraivogel et al? I am not sure. Unlike the TAP-seq experiment,
we have measured the whole transcriptome (a total of 17107 genes).

# At-scale enhancer screen

As far as I can tell, lentiviral transfection was done with libraries
targeting enhancers on chromosomes 8 and 11, but library preparation,
sequencing, and association analysis were done separately for these two
chromosomes. For example, only genes on chromosome 8 (chromosome 11)
were targeted in the sequencing for chromosome 8 (chromosome 11).
Therefore, the data for these two chromosomes are processed as separate
experiments in separate `ondisc` matrices.

## Experimental design

The gRNAs used for the at scale screen are a superset of those used for
the ground truth experiment In particular, the same set of 30 negative
control gRNAs were used for chromosomes 8 and 11 as in the ground truth
experiment. Furthermore, the
![4 \\times 14 = 56](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4%20%5Ctimes%2014%20%3D%2056 "4 \times 14 = 56")
positive controls (targeting 10 TSSs and 4 enhancers) from the ground
truth experiment were also used in the chromosome 8 screen (all the
genes targeted in the ground truth experiment were on chromosome 8,
while the enhancers targeted in this experiment came from a mix of
chromosomes 8 and 11). There were also
![4 \\times 14 = 56](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4%20%5Ctimes%2014%20%3D%2056 "4 \times 14 = 56")
positive controls (targeting 10 TSSs and 4 enhancers) used in the
chromosome 11 screen. The TSSs are of genes from chromosome 11 (and thus
distinct from the chromosome 8 screen), while the enhancers are the same
as those used in the chromosome 8 screen. In addition to the positive
and negative controls, the chromosome 8 screen had 4035 gRNAs targeting
1019 enhancers and the chromosome 11 screen had 3019 gRNAs targeting 767
enhancers (most enhancers were targeted by four gRNAs, while a handful
were targeted by fewer than four).

The experimental design information is encoded in the feature covariates
of the gRNA ODM in the same format as described above for the ground
truth experiment.

## Details of chromosome 8 screen

``` r
# load the gRNA expression data
processed_gRNA_dir <- sprintf("%s/enhancer_screen_chr8/gRNA_expression", processed_dir)
gRNA_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
gRNA_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_expr_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)
gRNA_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 4119 features and 112260 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero, target, target_type, known_effect.

``` r
# load the gRNA assignment data
processed_gRNA_assignment_dir <- sprintf("%s/enhancer_screen_chr8/gRNA_assignment", processed_dir)
gRNA_assignment_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_assignment_dir)
gRNA_assignment_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_assignment_dir)
gRNA_assign_odm <- ondisc::read_odm(gRNA_assignment_odm_fp, gRNA_assignment_metadata_fp)
gRNA_assign_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 4119 features and 112260 cells.
    ##  A cell covariate matrix with columns n_nonzero.
    ##  A feature covariate matrix with columns n_nonzero, target, target_type, known_effect.

``` r
# load the gene expression data
processed_gene_dir <- sprintf("%s/enhancer_screen_chr8/gene", processed_dir)
gene_odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
gene_metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
gene_expr_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)
gene_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 72 features and 112260 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

This experiment has 112260 cells across 14 batches. The gRNA data come
in the form of expressions and are not thresholded. There are a total of
4119 gRNAs. All control gRNAs are present, whereas 2 enhancer-targeting
gRNAs are missing from the data. Note that only 72 chromosome 8 genes
were targeted.

## Details of chromosome 11 screen

``` r
# load the gRNA expression data
processed_gRNA_dir <- sprintf("%s/enhancer_screen_chr11/gRNA_expression", processed_dir)
gRNA_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
gRNA_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_expr_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)
gRNA_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 3103 features and 120310 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero, target, target_type, known_effect.

``` r
# load the gRNA assignment data
processed_gRNA_assignment_dir <- sprintf("%s/enhancer_screen_chr11/gRNA_assignment", processed_dir)
gRNA_assignment_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_assignment_dir)
gRNA_assignment_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_assignment_dir)
gRNA_assign_odm <- ondisc::read_odm(gRNA_assignment_odm_fp, gRNA_assignment_metadata_fp)
gRNA_assign_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 3103 features and 120310 cells.
    ##  A cell covariate matrix with columns n_nonzero.
    ##  A feature covariate matrix with columns n_nonzero, target, target_type, known_effect.

``` r
# load the gene expression data
processed_gene_dir <- sprintf("%s/enhancer_screen_chr11/gene", processed_dir)
gene_odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
gene_metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
gene_expr_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)
gene_expr_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 82 features and 120310 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

This experiment has 120310 cells across 13 batches. The gRNA data come
in the form of expressions and are not thresholded. There are a total of
3103 gRNAs. All negative control gRNAs are present, whereas 1 positive
control and 1 enhancer-targeting gRNAs are missing from the data. Note
that only 82 chromosome 11 genes were targeted.
