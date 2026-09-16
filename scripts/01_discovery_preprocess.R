#!/usr/bin/env Rscript

# Discovery preprocessing and gene-level input construction.
# No cross-context analysis is performed here.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
phase_root <- normalizePath(analysis_root_env, mustWork = TRUE)

suppressPackageStartupMessages({
  library(affy)
})

processed_dir <- file.path(phase_root, "data", "processed_round2a")
qc_dir <- file.path(phase_root, "results", "round2a", "qc")
report_dir <- file.path(phase_root, "reports")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

write_matrix <- function(mat, path, id_name = "gene_symbol") {
  out <- data.frame(id = rownames(mat), mat, check.names = FALSE)
  names(out)[1] <- id_name
  write_tsv(out, path)
}

read_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  x <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- x[[1]]
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  mat
}

read_geo_annotation <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  header <- grep("^ID\\t", lines)[1]
  if (is.na(header)) stop("GEO annotation table header not found: ", path)
  ann <- read.delim(textConnection(lines[header:length(lines)]), check.names = FALSE,
                    quote = "", comment.char = "", stringsAsFactors = FALSE)
  required <- c("ID", "Gene symbol", "Gene ID")
  if (!all(required %in% names(ann))) stop("Required annotation fields missing: ", path)
  ann[, required]
}

clean_mapping <- function(annotation, probe_ids) {
  ann <- annotation[match(probe_ids, annotation$ID), , drop = FALSE]
  symbol <- trimws(ann$`Gene symbol`)
  gene_id <- trimws(ann$`Gene ID`)
  reliable <- !is.na(symbol) & nzchar(symbol) & symbol != "---" &
    !grepl("///", symbol, fixed = TRUE)
  symbol[!reliable] <- NA_character_
  gene_id[is.na(gene_id) | !nzchar(gene_id) | gene_id == "---" |
            grepl("///", gene_id, fixed = TRUE)] <- NA_character_
  data.frame(probe_id = probe_ids, gene_symbol = symbol, gene_id = gene_id,
             reliable_mapping = reliable, stringsAsFactors = FALSE)
}

summarize_gene_level <- function(probe_matrix, annotation, unit) {
  mapping <- clean_mapping(annotation, rownames(probe_matrix))
  keep <- mapping$reliable_mapping
  mapped_matrix <- probe_matrix[keep, , drop = FALSE]
  mapped <- mapping[keep, , drop = FALSE]
  variances <- apply(mapped_matrix, 1, var, na.rm = TRUE)
  variances[!is.finite(variances)] <- -Inf
  gene_order <- unique(mapped$gene_symbol)
  selected_index <- vapply(gene_order, function(g) {
    idx <- which(mapped$gene_symbol == g)
    idx[which.max(variances[idx])]
  }, integer(1))
  main <- mapped_matrix[selected_index, , drop = FALSE]
  rownames(main) <- mapped$gene_symbol[selected_index]
  selected <- mapped[selected_index, , drop = FALSE]
  selected$across_sample_variance <- variances[selected_index]
  selected$selection_rule <- "highest across-sample variance probe per gene"

  by_gene <- split(seq_len(nrow(mapped_matrix)), mapped$gene_symbol)
  median_mat <- t(vapply(by_gene, function(idx) {
    if (length(idx) == 1L) return(mapped_matrix[idx, ])
    apply(mapped_matrix[idx, , drop = FALSE], 2, median, na.rm = TRUE)
  }, numeric(ncol(mapped_matrix))))
  median_mat <- median_mat[rownames(main), , drop = FALSE]

  duplicate_genes <- sum(table(mapped$gene_symbol) > 1L)
  summary <- data.frame(
    unit = unit,
    original_probe_count = nrow(probe_matrix),
    mapped_probe_count = sum(keep),
    unmapped_or_ambiguous_probe_count = sum(!keep),
    unique_gene_count = length(unique(mapped$gene_symbol)),
    duplicate_mapped_gene_count = duplicate_genes,
    final_gene_count = nrow(main),
    primary_method = "highest across-sample variance probe",
    sensitivity_input = "median across probes per gene",
    stringsAsFactors = FALSE
  )
  list(main = main, median = median_mat, selected = selected,
       mapping = mapping, summary = summary)
}

open_png_pdf <- function(stem, code, width = 11, height = 8.5) {
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 180)
  code(); dev.off()
  pdf(paste0(stem, ".pdf"), width = width, height = height, useDingbats = FALSE)
  code(); dev.off()
}

manifest <- read.delim(file.path(phase_root, "data", "metadata", "sample_manifest.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)

# GSE180687: all 13 CEL files participate in RMA.
cel_dir <- file.path(phase_root, "data", "raw_round1_signoff", "GSE180687", "extracted")
cel_files <- sort(list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE))
if (length(cel_files) != 13L) stop("Expected 13 GSE180687 CEL files")
affy_raw <- ReadAffy(filenames = cel_files)
eset <- rma(affy_raw, verbose = FALSE)
g180_matrix <- exprs(eset)
g180_gsm <- sub("_.*", "", basename(sampleNames(eset)))
colnames(g180_matrix) <- g180_gsm
if (!identical(sort(g180_gsm), sort(manifest$GSM[manifest$dataset == "GSE180687"]))) {
  stop("GSE180687 RMA sample identities do not match manifest")
}
g180_matrix <- g180_matrix[, manifest$GSM[manifest$dataset == "GSE180687"], drop = FALSE]
write_matrix(g180_matrix, file.path(processed_dir, "GSE180687_RMA_expression.tsv"), "probe_id")

g180_meta <- manifest[match(colnames(g180_matrix), manifest$GSM), ]
g180_cols <- c("GSM", "sample_title", "treatment_state", "resistance_state", "biological_group",
               "biological_replicate", "include_exclude", "exclusion_reason")
write_tsv(g180_meta[, g180_cols], file.path(processed_dir, "GSE180687_endothelial_metadata.tsv"))

group_colors <- setNames(c("#4C78A8", "#59A14F", "#E15759"),
                         c("control endothelial", "sensitive endothelial", "resistant endothelial"))
sample_colors <- unname(group_colors[g180_meta$biological_group])
open_png_pdf(file.path(qc_dir, "GSE180687_RMA_expression_distribution"), function() {
  old <- par(mfrow = c(1, 2), mar = c(8, 4, 3, 1)); on.exit(par(old))
  boxplot(as.data.frame(g180_matrix), las = 2, cex.axis = 0.7, col = sample_colors,
          outline = FALSE, ylab = "RMA log2 expression", main = "GSE180687 RMA boxplot")
  ds <- lapply(seq_len(ncol(g180_matrix)), function(i) density(g180_matrix[, i]))
  plot(ds[[1]], col = sample_colors[1], xlab = "RMA log2 expression", ylab = "Density",
       main = "GSE180687 RMA density")
  for (i in 2:length(ds)) lines(ds[[i]], col = sample_colors[i])
  legend("topright", legend = names(group_colors), col = group_colors, lwd = 2, cex = 0.7, bty = "n")
})

g180_pca <- prcomp(t(g180_matrix), center = TRUE, scale. = FALSE)
open_png_pdf(file.path(qc_dir, "GSE180687_RMA_PCA"), function() {
  plot(g180_pca$x[, 1], g180_pca$x[, 2], col = sample_colors, pch = 16, cex = 1.3,
       xlab = paste0("PC1 (", round(100 * summary(g180_pca)$importance[2, 1], 1), "%)"),
       ylab = paste0("PC2 (", round(100 * summary(g180_pca)$importance[2, 2], 1), "%)"),
       main = "GSE180687 raw CEL → RMA PCA")
  text(g180_pca$x[, 1], g180_pca$x[, 2], labels = colnames(g180_matrix), pos = 3, cex = 0.65)
  legend("bottomleft", legend = names(group_colors), col = group_colors, pch = 16, cex = 0.7, bty = "n")
})

g180_cor <- cor(g180_matrix)
open_png_pdf(file.path(qc_dir, "GSE180687_RMA_sample_correlation"), function() {
  heatmap(g180_cor, Rowv = NA, Colv = NA, symm = TRUE, scale = "none",
          col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(101),
          margins = c(8, 8), cexRow = 0.7, cexCol = 0.7,
          main = "GSE180687 RMA sample correlation")
})
write_tsv(data.frame(GSM = rownames(g180_cor), g180_cor, check.names = FALSE),
          file.path(qc_dir, "GSE180687_RMA_sample_correlation.tsv"))

median_cor <- vapply(seq_len(nrow(g180_cor)), function(i) median(g180_cor[i, -i]), numeric(1))
pca_z <- sapply(seq_len(min(5, ncol(g180_pca$x))), function(j) {
  d <- mad(g180_pca$x[, j], constant = 1.4826)
  if (!is.finite(d) || d == 0) rep(0, nrow(g180_pca$x)) else abs((g180_pca$x[, j] - median(g180_pca$x[, j])) / d)
})
max_pca_z <- apply(pca_z, 1, max)
g180_qc <- data.frame(GSM = colnames(g180_matrix), biological_group = g180_meta$biological_group,
                      median_pairwise_correlation = median_cor,
                      max_abs_robust_z_PC1_to_PC5 = max_pca_z,
                      correlation_outlier = median_cor < median(median_cor) - 3 * mad(median_cor, constant = 1.4826),
                      PCA_outlier = max_pca_z > 4, stringsAsFactors = FALSE)
g180_qc$QC_decision <- ifelse(g180_qc$correlation_outlier & g180_qc$PCA_outlier, "REVIEW_REQUIRED",
                              ifelse(g180_qc$correlation_outlier | g180_qc$PCA_outlier, "CAUTION", "PASS"))
write_tsv(g180_qc, file.path(qc_dir, "GSE180687_RMA_sample_qc.tsv"))

annotations <- list(
  GPL10558 = file.path(phase_root, "data", "raw", "GSE76068", "official", "GPL10558.annot.gz"),
  GPL6885 = file.path(phase_root, "data", "raw", "GSE76068", "official", "GPL6885.annot.gz"),
  GPL6244 = file.path(phase_root, "data", "raw", "GSE73571", "official", "GPL6244.annot.gz"),
  GPL1261 = file.path(phase_root, "data", "raw", "GSE180687", "official", "GPL1261.annot.gz"),
  GPL6887 = file.path(phase_root, "data", "raw", "external_platform_annotations", "GPL6887.annot.gz"),
  GPL6884 = file.path(phase_root, "data", "raw", "external_platform_annotations", "GPL6884.annot.gz")
)
annotation_tables <- lapply(annotations, read_geo_annotation)

units <- list(
  GSE76068_tumor = list(dataset = "GSE76068", compartment = "tumor", species = "Homo sapiens", platform = "GPL10558",
    source = "data/processed/GSE76068_tumor_biological_level_proposed_expression.tsv.gz", matrix = read_matrix(file.path(phase_root, "data", "processed", "GSE76068_tumor_biological_level_proposed_expression.tsv.gz"))),
  GSE76068_stroma = list(dataset = "GSE76068", compartment = "stroma", species = "Mus musculus", platform = "GPL6885",
    source = "data/processed/GSE76068_stroma_biological_level_proposed_expression.tsv.gz", matrix = read_matrix(file.path(phase_root, "data", "processed", "GSE76068_stroma_biological_level_proposed_expression.tsv.gz"))),
  GSE73571_tumor = list(dataset = "GSE73571", compartment = "tumor", species = "Homo sapiens", platform = "GPL6244",
    source = "data/processed/GSE73571_tumor_array_level_expression.tsv.gz", matrix = read_matrix(file.path(phase_root, "data", "processed", "GSE73571_tumor_array_level_expression.tsv.gz"))),
  GSE180687_endothelial = list(dataset = "GSE180687", compartment = "endothelial", species = "Mus musculus", platform = "GPL1261",
    source = "13 raw CEL files → affy::rma", matrix = g180_matrix),
  GSE64472_stroma = list(dataset = "GSE64472", compartment = "stroma", species = "Mus musculus", platform = "GPL6887",
    source = "data/processed/GSE64472_stroma_array_level_expression.tsv.gz", matrix = read_matrix(file.path(phase_root, "data", "processed", "GSE64472_stroma_array_level_expression.tsv.gz"))),
  GSE26644_tumor = list(dataset = "GSE26644", compartment = "tumor", species = "Homo sapiens", platform = "GPL6884",
    source = "data/processed/GSE26644_tumor_array_level_expression.tsv.gz", matrix = read_matrix(file.path(phase_root, "data", "processed", "GSE26644_tumor_array_level_expression.tsv.gz"))),
  GSE26644_stroma = list(dataset = "GSE26644", compartment = "stroma", species = "Mus musculus", platform = "GPL6887",
    source = "data/processed/GSE26644_stroma_array_level_expression.tsv.gz", matrix = read_matrix(file.path(phase_root, "data", "processed", "GSE26644_stroma_array_level_expression.tsv.gz")))
)

mapping_summaries <- list()
input_summaries <- list()
for (unit in names(units)) {
  u <- units[[unit]]
  mat <- u$matrix
  if (unit == "GSE73571_tumor") {
    keep_gsm <- manifest$GSM[manifest$dataset == "GSE73571" & manifest$include_exclude == "CORE_INCLUDE"]
    mat <- mat[, keep_gsm, drop = FALSE]
  }
  if (unit %in% c("GSE64472_stroma", "GSE26644_tumor", "GSE26644_stroma")) {
    keep_gsm <- manifest$GSM[manifest$dataset == u$dataset & manifest$compartment == u$compartment & manifest$include_exclude == "CORE_INCLUDE"]
    mat <- mat[, keep_gsm, drop = FALSE]
  }

  if (grepl("^GSE76068", unit)) {
    sample_meta <- data.frame(
      analysis_sample_id = colnames(mat),
      GSM = NA_character_,
      treatment_state = ifelse(grepl("pretreatment", colnames(mat)), "pretreatment",
                               ifelse(grepl("on-treatment_response", colnames(mat)), "response", "escape")),
      resistance_state = ifelse(grepl("pretreatment", colnames(mat)), "baseline",
                                ifelse(grepl("on-treatment_response", colnames(mat)), "sensitive/responding", "acquired resistance/escape")),
      biological_group = ifelse(grepl("pretreatment", colnames(mat)), "pretreatment",
                                ifelse(grepl("on-treatment_response", colnames(mat)), "response", "escape")),
      biological_replicate = sub(".*_(BR[1-4])$", "\\1", colnames(mat)),
      technical_replicate_handling = ifelse(unit == "GSE76068_tumor",
        "approved passing a/b arrays averaged on log2 scale",
        "approved passing arrays averaged; GSM1973644/GSM1973648 excluded and passing mates retained"),
      include_exclude = "CORE_INCLUDE", stringsAsFactors = FALSE
    )
  } else {
    mm <- manifest[match(colnames(mat), manifest$GSM), ]
    sample_meta <- data.frame(analysis_sample_id = colnames(mat), GSM = mm$GSM,
                              treatment_state = mm$treatment_state, resistance_state = mm$resistance_state,
                              biological_group = mm$biological_group,
                              biological_replicate = mm$biological_replicate,
                              technical_replicate_handling = "none declared; one array per biological sample",
                              include_exclude = mm$include_exclude, stringsAsFactors = FALSE)
  }

  write_matrix(mat, file.path(processed_dir, paste0(unit, "_final_probe_expression.tsv")), "probe_id")
  write_tsv(sample_meta, file.path(processed_dir, paste0(unit, "_metadata.tsv")))
  gl <- summarize_gene_level(mat, annotation_tables[[u$platform]], unit)
  write_matrix(gl$main, file.path(processed_dir, paste0(unit, "_gene_level_expression.tsv")))
  write_matrix(gl$median, file.path(processed_dir, paste0(unit, "_gene_level_expression_median_probe.tsv")))
  write_tsv(gl$selected, file.path(processed_dir, paste0(unit, "_selected_probe_per_gene.tsv")))
  mapping_summaries[[unit]] <- cbind(dataset = u$dataset, compartment = u$compartment,
                                     species = u$species, platform = u$platform, gl$summary)

  excluded <- if (unit == "GSE76068_stroma") "GSM1973644;GSM1973648" else if (unit == "GSE73571_tumor")
    paste(manifest$GSM[manifest$dataset == "GSE73571" & manifest$include_exclude == "EXCLUDE"], collapse = ";") else if (unit == "GSE180687_endothelial")
      "none; controls retained for joint RMA/QC" else "none in approved compartment"
  input_summaries[[unit]] <- data.frame(
    unit = unit, dataset = u$dataset, compartment = u$compartment, species = u$species,
    platform = u$platform, input_source = u$source,
    preprocessing = if (unit == "GSE180687_endothelial") "all 13 raw CEL files jointly processed by affy::rma" else "approved Round 1 processed log2 matrix",
    technical_replicate_handling = paste(unique(sample_meta$technical_replicate_handling), collapse = "; "),
    final_sample_count = ncol(mat), biological_replicate_definition = if (grepl("^GSE76068", unit)) "BR1-BR4 matched across states" else "one GEO biological sample",
    excluded_samples = excluded, gene_annotation = paste0(u$platform, " GEO annotation; reliable single-symbol mappings only"),
    final_gene_count = nrow(gl$main), stringsAsFactors = FALSE
  )

  cors <- cor(gl$main)
  med_cor <- vapply(seq_len(nrow(cors)), function(i) median(cors[i, -i]), numeric(1))
  pc <- prcomp(t(gl$main), center = TRUE, scale. = FALSE)
  q <- data.frame(sample_id = colnames(gl$main), biological_group = sample_meta$biological_group,
                  median_pairwise_correlation = med_cor, PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  array_variance = apply(gl$main, 2, var), stringsAsFactors = FALSE)
  write_tsv(q, file.path(qc_dir, paste0(unit, "_formal_input_sample_qc.tsv")))
}

mapping_summary <- do.call(rbind, mapping_summaries)
input_summary <- do.call(rbind, input_summaries)
write_tsv(mapping_summary, file.path(report_dir, "probe_to_gene_mapping_summary.tsv"))
write_tsv(input_summary, file.path(report_dir, "round2a_input_summary.tsv"))

provenance_lines <- c(
  "# Round 2A input provenance", "", paste0("Generated: ", Sys.Date()), "",
  "Round 1 audit inputs were not overwritten. Formal Round 2A files are stored under `data/processed_round2a/`.", "",
  "## Locked preprocessing rules", "",
  "- GSE76068 tumor/stroma: approved deposited log2 matrices collapsed to one biological profile per BR/state; matched BR1–BR4 design is required.",
  "- GSE76068 stroma: GSM1973644 and GSM1973648 excluded; GSM1973643 and GSM1973647 retain BR2/BR4 response.",
  "- GSE180687: all 13 raw GPL1261 CEL files jointly processed with `affy::rma`; controls retained for normalization/QC but no additional core contrast.",
  "- GSE73571, GSE64472 mouse/stroma, and GSE26644 tumor/stroma: approved Round 1 processed matrices reused without redundant raw preprocessing.",
  "- Probe-to-gene primary rule: reliable single-symbol GEO mappings followed by highest across-sample variance probe per gene.",
  "- Sensitivity input: median expression across mapped probes per gene; prepared but not analyzed in Round 2A.", "",
  "## Final inputs", "",
  "| Unit | Species | Platform | Samples | Final genes | Source/preprocessing | Exclusions |",
  "| --- | --- | --- | ---: | ---: | --- | --- |"
)
for (i in seq_len(nrow(input_summary))) {
  x <- input_summary[i, ]
  provenance_lines <- c(provenance_lines, paste0("| ", x$unit, " | ", x$species, " | ", x$platform, " | ",
    x$final_sample_count, " | ", x$final_gene_count, " | ", x$preprocessing, " | ", x$excluded_samples, " |"))
}
provenance_lines <- c(provenance_lines, "", "Detailed sample metadata and selected-probe tables are stored beside each expression matrix.")
writeLines(provenance_lines, file.path(report_dir, "round2a_input_provenance.md"))
writeLines(capture.output(sessionInfo()), file.path(report_dir, "round2a_preprocessing_sessionInfo.txt"))
message("Round 2A preprocessing completed.")
