# Differential Expression Analysis


# LIBRARIES --------------------------------------------------------------------
library(yaml)
library(cli)
library(BiocParallel)    # required for parallelization
library(ggplot2)
library(DESeq2)
library(MASS)
library(ashr)
library(writexl)
library(dynamicTreeCut)
library(dplyr)
library(pheatmap)

# FUNCTION DEFINITION ----------------------------------------------------------
plot_volcano <- function(res, title, padj_cutoff = 0.05, lfc_cutoff = 1) {
  df <- as.data.frame(res)
  df$GeneID <- rownames(df)
  
  # Safety NA filtering
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange),]
  
  df$negLog10Padj <- -log10(df$padj)
  
  # Classification
  df$Status <- "Not significant"
  df$Status[df$padj <= padj_cutoff & df$log2FoldChange >= lfc_cutoff] <- "Significant_up"
  df$Status[df$padj <= padj_cutoff & df$log2FoldChange <= -lfc_cutoff] <- "Significant_down"
  
  p <- ggplot(df, aes(x = log2FoldChange, y = negLog10Padj)) +
    geom_point(aes(color = Status), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("grey70", "blue", "red")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black") +
    theme_bw() +
    labs(title = title,
         x = "log2 Fold Change",
         y = "-log10 adjusted p-value") +
    theme(legend.title = element_blank())
  
  return(p)
}

# Definition of function required to check outputs
check_output <- function(file_path) {
  if (file.exists(file_path)) {
    cli_alert_success(
      paste0("File created successfully -> ", file_path)
    )
  } else {
    cli_abort(
      paste0("Failed to create file -> ", file_path)
    )
  }
}

# Function to get shrunk results
get_shrunk_results <- function(dds_fit, group1, group2) {
  lfcShrink(
    dds_fit,
    contrast = c("CONDITION", group1, group2),
    type = "ashr"
  )
}

# ARGUMENTS AND CONFIG ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)    # Enable arguments

counts_file <- args[1]      # count matrix file
metadata_file <- args[2]    # metadata file
config_file <- args[3]

config <- yaml::read_yaml(config_file)

## PARALLELIZATION -------------------------------------------------------------
n_cores <- as.numeric(Sys.getenv("THREADS", unset = 1))    # set the number of cores used during analysis

if (.Platform$OS.type == "unix") {                         # Enable the use of multiple cores depending on the Operative System
  register(MulticoreParam(n_cores))
} else {
  register(SnowParam(n_cores))
}

# ---------------------------- SET PARAMETERS ----------------------------------

output_dir <- paste0("runs/", config$run_name, "/", config$output_dir, "/DIFF_EXPR/")
attribute <- config$featurecounts$attribute
geature_type <- config$featurecounts$feature_type
padj_cutoff <- config$de_analysis$padj_cutoff
lfc_cutoff <- config$de_analysis$lfc_cutoff
reference_condition <- config$reference_condition
cli_inform("DEG filtering thresholds:")
cli_inform(paste0("  padj <= ", padj_cutoff))
cli_inform(paste0("  |log2FC| >= ", lfc_cutoff))
cli_inform("Reference condition set to: ", reference_condition)

# INPUT DATA -------------------------------------------------------------------

raw_data <- read.delim(counts_file, skip = 1)           # import raw data
cli_alert_success("INFO: Data imported successfully")
metadata <- read.delim(metadata_file, header = TRUE)    # import metadata
cli_inform("INFO: metadata file imported:\n")
metadata
metadata$CONDITION <- factor(metadata$CONDITION)        # Set conditions in metadata file as a factor

# write diagnostic file
write.csv(
  data.frame(condition_levels = levels(metadata$CONDITION)),
  paste0(output_dir, "debug/factor_levels_used.csv"),
  row.names = FALSE)
# Check reference condition
if(is.null(reference_condition) || reference_condition == "") {
  cli_abort("Reference condition must be defined in config.yaml")
}

if (!reference_condition %in% levels(metadata$CONDITION)) {
  
  cli_alert_danger(paste(
    "ERROR: reference_condition = ", reference_condition,
    "not found in metadata$CONDITION.\n",
    "",
    paste0("Available levels: ", paste(levels(metadata$CONDITION), collapse = ", ")),
    "\n A diagnostic file has been generated:\n",
    paste0(output_dir, "debug(factor_levels_used.csv\n\n"),
    "Fix: update config.yaml at line 6 or check metadata$CONDITION spelling"
  ))
}

metadata$CONDITION <- relevel(metadata$CONDITION, ref = reference_condition)
cli_inform("Using reference condition: ", reference_condition)
cli_inform("Factor levels after relevel: ")
print(levels(metadata$CONDITION))

write.csv(
  data.frame(
    reference_used = reference_condition,
    final_levels = paste(levels(metadata$CONDITION), collapse = ",")
  ),
  paste0(output_dir, "debug/reference_after_relevel.csv"),
  row.names = FALSE
)

# PREFILTERING -----------------------------------------------------------------
## Setting identifier column name ----------------------------------------------
# 1. Take attribute from config if present
if (!is.null(config$featurecounts$attribute) &&
    config$featurecounts$attribute != "") {
    attribute <- config$featurecounts$attribute
    cli_inform(paste("PARAMETER: Using attribute from config:", attribute))
# 2. Fallback - feature_type will be used
} else if (!is.null(config$featurecounts$feature_type) &&
           config$featurecounts$feature_type != "") {
  attribute <- config$featurecounts$feature_type
  cli_alert_warning(paste("WARNING: attribute missing, using feature_type:", attribute))

# 3. Fallback - using column name
} else {
  attribute <- colnames(raw_data)[1]
  cli_alert_warning(paste("WARNING: using column name as attribute:", attribute))
}

## Clean columns and sample names ----------------------------------------------
counts <- raw_data[, -(2:6)]         # Remove surplus columns
rownames(counts) <- counts$Geneid    # Set GeneID as a rowname
counts <- counts[, -1]               # Remove GeneID column

# This step is required to be sure that sample names do not contain relative path
pattern_corretto <- paste0("^runs\\.", config$run_name, "\\.", config$output_dir, "\\.ALIGNMENT\\.(.*)\\.bam$")
colnames(counts) <- sub(pattern_corretto, "\\1", colnames(counts))

# Verify if counts and metadata files are correct
if (!all(colnames(counts) %in% metadata$SAMPLE) || !all(metadata$SAMPLE %in% colnames(counts))) {
  cli_abort("Mismatch between counts columns and metadata samples")
  }
# Verify the presence of duplicates
if (any(duplicated(rownames(counts)))) {
  cli_abort("Duplicate gene IDs found")
}
# Verify the presence of NA values
if (any(is.na(counts))) {
  cli_abort("NA values found in counts matrix")
}
# Verify the presence of non-numeric values
if(!all(sapply(counts, is.numeric))) {
  cli_abort("Counts matrix contains non-numeric values")
}

## Reorder and Check After Matching --------------------------------------------
metadata <- metadata[match(colnames(counts), metadata$SAMPLE), ]
if (!all(metadata$SAMPLE == colnames(counts))) {
  cli_abort("Sample order mismatch after matching")
}

# Convert decimal values in integers
counts[] <- lapply(counts, function(x) {
  as.integer(round(as.numeric(x)))
})

# CREATE DESeq2 DATASET --------------------------------------------------------
dds_raw <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design =~ CONDITION)

print(paste("The dimension of the dataset is:", paste(dim(dds_raw), collapse = " x ")))
dim(dds_raw[rowSums(counts(dds_raw)) >= 10,])    # Get dimension of dataset without rows whose sum is <10
cli_inform("Filtering low count rows")
keep <- rowSums(counts(dds_raw)) >= 10
print(paste("The dimension of the dataset excluding low count rows is:", paste(dim(dds_raw[keep,]), collapse = " x ")))
filtered_dds <- dds_raw[keep,]

# Differential Expression Analysis ---------------------------------------------
dds_fit <- DESeq(filtered_dds, parallel = TRUE)

saveRDS(dds_fit, paste0(output_dir, "rds_files/dds_fit.rds"))    # Save dds_fit
check_output(paste0(output_dir, "rds_files/dds_fit.rds"))

condition_levels <- levels(metadata$CONDITION)    # Obtaining conditions
comparisons <- combn(condition_levels, 2, simplify = FALSE)    # Obtaining combinations

## Results ----------------------------------------------------------------------
results_list <- list()

if (config$de_analysis$lfc_source == "shrunk") {
  cli_inform("INFO: shrinkage is used")
}

for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  comp_name <- paste(group1, "VS", group2, sep = "_")
  
  cli_inform("Running: {comp_name}")
  
  # Raw results
  
  res <- results(dds_fit, contrast = c("CONDITION", group1, group2))
  results_list[[comp_name]] <- res
}

# Save results file
saveRDS(results_list, paste0(output_dir, "rds_files/results_list.rds"))
check_output(paste0(output_dir, "rds_files/results_list.rds"))

## DE genes extraction ---------------------------------------------------------
sig_results_list <- list()
for (name in names(results_list)) {
  res <- as.data.frame(results_list[[name]])    # Convert in dataframe
  res$GeneID <- rownames(res)
  res <- res[!is.na(res$padj) &              # Remove NA values
             !is.na(res$log2FoldChange) &    # Remove NA values
             res$padj <= padj_cutoff &          # Filtering padj
             abs(res$log2FoldChange) >= lfc_cutoff, ]    # Filtering Log2FC
  sig_results_list[[name]] <- res
}

all_sig_genes <- unique(unlist(lapply(sig_results_list, function(x) x$GeneID)))

# Conversion of all_sig_genes in data.frame for saving
for (name in names(sig_results_list)) {
  write.csv(sig_results_list[[name]],
            paste0(output_dir, "tables/", name, ".csv"),
            row.names = FALSE)
  check_output(paste0(output_dir, "tables/", name, ".csv"))
}
write_xlsx(sig_results_list, path = paste0(output_dir, "tables/significative_genes.xlsx"))
check_output(paste0(output_dir, "tables/significative_genes.xlsx"))

# Create Log2FC matrix
log2fc_matrix <- lapply(names(results_list), function(name) {
  
  res <- as.data.frame(results_list[[name]])
  
  data.frame(GeneID = rownames(res), log2FC = res$log2FoldChange, check.names = FALSE)
})
names(log2fc_matrix) <- names(results_list)
log2fc_merged <- Reduce(function(x, y) {
  merge(x, y, by = "GeneID", all = TRUE)
}, log2fc_matrix)
colnames(log2fc_merged)[-1] <- names(results_list)
log2fc_merged <- log2fc_merged[order(log2fc_merged$GeneID), ]

write.csv(log2fc_merged, paste0(output_dir, "tables/log2FC_matrix.csv"), row.names = FALSE)
check_output(paste0(output_dir, "tables/log2FC_matrix.csv"))
write_xlsx(log2fc_merged, paste0(output_dir, "tables/log2FC_matrix.xlsx"))
check_output(paste0(output_dir, "tables/log2FC_matrix.xlsx"))

# VST (VARIANCE STABILIZING TRANSFORMATION) -------------------------------------
# Transformation of data and matrix creation

vsdata <- vst(dds_fit, blind = FALSE)    # VST transformation
vsdata_matrix <- assay(vsdata)                  # Create a matrix for clustering

# PCA Plot ---------------------------------------------------------------------
pcaData <- plotPCA(vsdata, intgroup = "CONDITION", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = CONDITION)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 0.6, linetype = "dashed") +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()

# Save PCA plot
pca_path <- paste0(output_dir, "plots/PCA.pdf")
ggsave(pca_path, plot = p, width = 8, height = 6)
print("PCA plot correctly saved")

# Save VST data 
saveRDS(vsdata, paste0(output_dir, "rds_files/vsdata.rds"))
check_output(paste0(output_dir, "rds_files/vsdata.rds"))
write.csv(vsdata_matrix, file = paste0(output_dir, "tables/vsdata_matrix.csv"), quote = FALSE)
check_output(paste0(output_dir, "tables/vsdata_matrix.csv"))
vs_df <- data.frame(GeneID = rownames(vsdata_matrix), vsdata_matrix, check.names = FALSE) # Set GeneID as 1st column
write_xlsx(vs_df, path = paste0(output_dir, "tables/vsdata_matrix.xlsx"))
check_output(paste0(output_dir, "tables/vsdata_matrix.xlsx"))

# VOLCANO PLOTS ----------------------------------------------------------------
volcano_plots <- list()
for (comp_name in names(results_list)) {
  cli_inform("Progessing: {comp_name}")
  comps <- strsplit(comp_name, "_VS_")[[1]]
  group1 <- comps[1]
  group2 <- comps[2]

  cli_inform("{group1} - {group2}")
  
  # Select LFC source for plotting 
  if (config$de_analysis$lfc_source == "shrunk") {
    res_plot <- results_list[[comp_name]]
  } else {
    res_plot <- get_shrunk_results(dds_fit, group1, group2)
  }
  
  # Create volcano plot
  volcano <- plot_volcano(res = res_plot,
                          title = comp_name,
                          padj_cutoff = padj_cutoff,
                          lfc_cutoff = lfc_cutoff)
  volcano_plots[[comp_name]] <- volcano
  
  ggsave(paste0(output_dir, "plots/volcano", comp_name, ".pdf"), width = 210, height = 297, units = "mm")
  check_output(paste0(output_dir, "plots/volcano", comp_name, ".pdf"))
}

# CLUSTERING -------------------------------------------------------------------

gene_var <- rowVars(vsdata_matrix)    # compute gene variance
threshold <- quantile(gene_var, config$clustering$threshold)    # Define threshold (quantile-based cutoff)

# filter genes
vsdata_filtered <- vsdata_matrix[gene_var > threshold,]

# Logs
cli_inform("Clustering genes retained: {nrow(vsdata_filtered)}")
cli_inform("Clustering genes removed: {sum(gene_var <= threshold)}")
cli_inform("Variance threshold used: {threshold}")

saveRDS(vsdata_filtered, paste0(output_dir, "rds_files/vsdata_filtered_clustering.rds"))
check_output(paste0(output_dir, "rds_files/vsdata_filtered_clustering.rds"))

write.csv(vsdata_filtered, paste0(output_dir, "tables/vsdata_filtered_clustering.csv"))
check_output(paste0(output_dir, "tables/vsdata_filtered_clustering.csv"))

# Distance computation
distance_condition <- dist(t(vsdata_filtered))
distance_gene <- as.dist(1- cor(t(vsdata_filtered), method = config$clustering$cor_method))

cli_inform("INFO: Performing gene clustering")
cluster_gene <- hclust(distance_gene, method = config$clustering$hclust_method)
saveRDS(cluster_gene, paste0(output_dir, "rds_files/hclust_genes.rds"))
pdf(paste0(output_dir, "plots/hclust_genes.pdf"))
plot(as.dendrogram(cluster_gene), horiz = TRUE)
dev.off()

cli_inform("INFO: Performing sample clustering")
cluster_samples <- hclust(distance_condition, method = config$clustering$hclust_method)
saveRDS(cluster_samples, paste0(output_dir, "rds_files/hclust_samples.rds"))
pdf(paste0(output_dir, "plots/hclust_samples.pdf"))
plot(cluster_samples)
dev.off()

# Default strategy: static cut is used for reproducibility and interpretability
cut_method <- config$clustering$cut_method

if (is.null(cut_method) || cut_method == "") {
  cut_method <- "static"
}

if (cut_method == "dynamic") {
  cli_inform("INFO: Using dynamic tree cut")
  
  deepSplit <- config$clustering$deepSplit
  if (is.null(deepSplit)) deepSplit <- 2
  
  clusters <- cutreeDynamic(dendro = as.dendrogram(cluster_gene),
                            distM = as.matrix(distance_gene),
                            deepSplit = deepSplit)
} else {
  cli_inform("INFO: Using static cut")
  k <- config$clustering$k
  if (is.null(k)) cli_abort("k must be defined when cut_method is not dynamic")
  clusters <- cutree(cluster_gene, k = k)
}
cli_inform("Clusters obtained:")
print(table(clusters))

cluster_assignment <- data.frame(geneIDs = rownames(vsdata_filtered),
                                 cluster = clusters)
write.csv(cluster_assignment, paste0(output_dir, "tables/cluster_assignment.csv"))
check_output(paste0(output_dir, "tables/cluster_assignment.csv"))
write_xlsx(cluster_assignment, paste0(output_dir, "tables/cluster_assignment.xlsx"))
check_output(paste0(output_dir, "tables/cluster_assignment.xlsx"))

## Ordering by cluster ----------------------------------------------------------
gene_order <- rownames(vsdata_filtered)[cluster_gene$order]
ordered_sig_genes <- gene_order[gene_order %in% all_sig_genes]

# HEATMAP ----------------------------------------------------------------------

## Subset expression matrix (filtered + ordered genes)

cluster_df <- cluster_assignment
colnames(cluster_df)[1] <- "geneIDs"

cluster_df <- cluster_df[match(rownames(vsdata_filtered), cluster_df$geneIDs), ]

gene_order <- cluster_df$geneIDs[order(cluster_df$cluster)]

heatmap_matrix <- vsdata_filtered[gene_order, ]

## Check
if (any(is.na(cluster_df$geneIDs))) {
  cli_alert_warning("Some genes in heatmap are missing cluster assignment")
}

## Pheatmap with samples--------------------------------------------------------

mat <- vsdata_filtered[gene_order, , drop = FALSE]

# Safety check
if(!all(gene_order %in% rownames(vsdata_filtered))) {
  cli_alert_warning("gene_order contains genes not in vsdata_filtered")
}

if (anyDuplicated(gene_order)) {
  cli_alert_warning("gene_order contains duplicated genes")
}

annotation_row <- data.frame(cluster = cluster_assignment$cluster)
rownames(annotation_row) <- cluster_assignment$geneIDs
annotation_row <- annotation_row[rownames(mat), , drop = FALSE]
# Annotation_colors
cluster_levels <- sort(unique(annotation_row$cluster))

# Setting color blinded-friendly colors
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

cluster_colors <- setNames(
  okabe_ito[seq_along(cluster_levels)],
  cluster_levels
)

# Setting color gradient from blue to red
annotation_colors <- list(cluster = cluster_colors)
my_colors <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

# Generating heatmap
cli_inform("INFO: Generating heatmap")

pdf(paste0(output_dir, "plots/heatmap_samples_vst.pdf"), width = 8.27, height = 11.69)
pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         color = my_colors,
         scale = "row",
         show_rownames = FALSE)
dev.off()
check_output(paste0(output_dir, "plots/heatmap_samples_vst.pdf"))

## Pheatmap with conditions ----------------------------------------------------
# Mean VST expression per condition
vsdata_condition <- sapply(condition_levels, function(cond) {
  rowMeans(vsdata_filtered[, metadata$CONDITION == cond, drop = FALSE])
})
vsdata_condition <- as.matrix(vsdata_condition)
vsdata_condition <- vsdata_condition[gene_order,, drop = FALSE]

## Heatmap
pdf(paste0(output_dir, "plots/heatmap_conditions_vst.pdf", width = 8.27, height = 11.69))

pheatmap(vsdata_condition, cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         color = my_colors,
         scale = "row",
         show_rownames = FALSE)
dev.off()

# Log2FoldChange Heatmap -------------------------------------------------------
## Create Log2FoldChange Matrix
log2fc_list <- lapply(comparisons, function(comp) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  comp_name <- paste(comp[1], "VS", comp[2], sep = "_")
  
  res_plot <- as.data.frame(get_shrunk_results(dds_fit, group1, group2))
  data.frame(GeneID = rownames(res_plot), log2FC = res_plot$log2FoldChange, check.names = FALSE)
})

# Assigning comparison names
comp_names <- sapply(comparisons, function(comp)
  paste(comp[1], "VS", comp[2], sep = "_"))

names(log2fc_list) <- comp_names

# Rename columns before merging
log2fc_list_named <- mapply(function(df, name) {
  colnames(df)[2] <- name
  df
}, log2fc_list, comp_names, SIMPLIFY = FALSE)

log2fc_matrix <- Reduce(function(x, y) {
  full_join(x, y, by = "GeneID")
}, log2fc_list_named)

log2fc_matrix <- as.data.frame(log2fc_matrix)
rownames(log2fc_matrix) <- log2fc_matrix$GeneID
log2fc_matrix$GeneID <- NULL

## DEG filter
deg_genes <- intersect(all_sig_genes, rownames(log2fc_matrix))

# Keep only DEGs present in matrix
log2fc_matrix <- log2fc_matrix[deg_genes, , drop = FALSE]

# Ordered genes
log2fc_matrix <- log2fc_matrix[match(gene_order, rownames(log2fc_matrix)),]
## Heatmap Log2FC
pdf(paste0(output_dir, "plots/heatmap_conditions_log2FC.pdf"), width = 8.27, height = 11.69)
pheatmap(log2fc_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         na_col = "grey90",
         show_rownames = FALSE)
dev.off()

check_output(paste0(output_dir, "plots/heatmap_conditions_log2FC.pdf"))
