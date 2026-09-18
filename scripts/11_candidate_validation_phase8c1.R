#!/usr/bin/env Rscript

# Candidate-dataset analyses identified after the discovery program freeze.
# Final evidence roles are assigned by config/candidate_dataset_disposition.tsv.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
bundle_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
analysis_root <- Sys.getenv("VEGF_ANALYSIS_ROOT", unset = bundle_root)
input_root <- Sys.getenv("VEGF_CANDIDATE_INPUT_ROOT", unset = file.path(analysis_root, "data", "raw", "candidate_validation"))
phase8c1_input <- file.path(input_root, "phase8c1")
legacy_library <- Sys.getenv("VEGF_R_LIBRARY", unset = "")
if (nzchar(legacy_library)) .libPaths(c(legacy_library, .libPaths()))
suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(msigdbr)
})

meta_dir <- phase8c1_input
out_root <- Sys.getenv("VEGF_OUTPUT_ROOT", unset = file.path(bundle_root, "outputs", "candidate_validation", "phase8c1"))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_series_matrix <- function(path) {
  z <- readLines(gzfile(path), warn = FALSE)
  a <- match("!series_matrix_table_begin", z) + 1L
  b <- match("!series_matrix_table_end", z) - 1L
  stopifnot(is.finite(a), is.finite(b), b >= a)
  read.delim(text = paste(z[a:b], collapse = "\n"), check.names = FALSE,
             stringsAsFactors = FALSE, quote = "\"", comment.char = "")
}

read_geo_annotation <- function(path) {
  read.delim(gzfile(path), skip = 27, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

read_soft_platform <- function(path) {
  z <- readLines(path, warn = FALSE)
  a <- match("!platform_table_begin", z) + 1L
  b <- match("!platform_table_end", z) - 1L
  stopifnot(is.finite(a), is.finite(b), b >= a)
  read.delim(text = paste(z[a:b], collapse = "\n"), check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "",
             fill = TRUE)
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("[[:space:]]*///.*$", "", x)
  x <- sub("[[:space:]]*//.*$", "", x)
  x <- sub("[;,|].*$", "", x)
  x <- toupper(trimws(x))
  x[x %in% c("", "---", "NA", "NULL")] <- NA_character_
  x
}

collapse_to_symbol <- function(expr, ids, annotation, id_col, symbol_col) {
  idx <- match(as.character(ids), as.character(annotation[[id_col]]))
  symbols <- clean_symbol(annotation[[symbol_col]][idx])
  keep <- !is.na(symbols) & apply(expr, 1, function(v) all(is.finite(v)))
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  iq <- apply(expr, 1, IQR)
  ord <- order(symbols, -iq)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  keep2 <- !duplicated(symbols)
  expr <- expr[keep2, , drop = FALSE]
  rownames(expr) <- symbols[keep2]
  expr
}

frozen <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

get_pathways <- function(species) {
  hallmark <- suppressWarnings(msigdbr(species = species, collection = "H"))
  hallmark$gene_symbol <- toupper(hallmark$gene_symbol)
  split(hallmark$gene_symbol, hallmark$gs_name)
}

run_contrast <- function(dataset, contrast_label, expr, numerator, denominator,
                         species, evidence_role, context_id) {
  samples <- c(denominator, numerator)
  x <- expr[, samples, drop = FALSE]
  group <- factor(c(rep("denominator", length(denominator)),
                    rep("numerator", length(numerator))),
                  levels = c("denominator", "numerator"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- eBayes(contrasts.fit(lmFit(x, design),
                              makeContrasts(numerator - denominator, levels = design)),
                robust = TRUE)
  tt <- topTable(fit, number = Inf, sort.by = "none")
  tt$gene_symbol <- rownames(tt)
  tt <- tt[, c("gene_symbol", setdiff(colnames(tt), "gene_symbol"))]

  ranks <- fit$t[, 1]
  names(ranks) <- rownames(fit$t)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(81465)
  fg <- fgseaMultilevel(pathways = get_pathways(species), stats = ranks,
                        minSize = 15, maxSize = 500, eps = 0)
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";",
                           FUN.VALUE = character(1))
  out <- fg[match(frozen, fg$pathway),
            c("pathway", "NES", "pval", "padj", "size", "leadingEdge")]
  out$dataset <- dataset
  out$context_id <- context_id
  out$contrast <- contrast_label
  out$n_numerator <- length(numerator)
  out$n_denominator <- length(denominator)
  out$direction <- ifelse(out$NES > 0, "positive", "negative")
  out$evidence_role <- evidence_role
  out$identified_date <- "2026-09-11"
  out <- out[, c("dataset", "context_id", "contrast", "pathway", "NES", "pval",
                 "padj", "size", "direction", "n_numerator", "n_denominator",
                 "evidence_role", "identified_date", "leadingEdge")]

  d <- file.path(out_root, dataset)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  safe <- gsub("[^A-Za-z0-9]+", "_", context_id)
  write_tsv(tt, file.path(d, paste0(dataset, "_", safe, "_gene_statistics.tsv")))
  write_tsv(data.frame(sample = samples, group = as.character(group)),
            file.path(d, paste0(dataset, "_", safe, "_sample_manifest.tsv")))
  out
}

results <- list()

# GSE121153: xenografts only; the two single in-vitro arrays are excluded.
m <- read_series_matrix(file.path(meta_dir, "GSE121153_series_matrix.txt.gz"))
ids <- m[[1]]
x <- as.matrix(m[, -1]); storage.mode(x) <- "double"
x <- normalizeBetweenArrays(log2(x + 1), method = "quantile")
a <- read_geo_annotation(file.path(input_root, "platforms", "GPL570.annot.gz"))
x <- collapse_to_symbol(x, ids, a, "ID", "Gene symbol")
results[[length(results) + 1L]] <- run_contrast(
  "GSE121153", "sorafenib-resistant xenografts vs parental xenografts", x,
  paste0("GSM34270", 17:21), c("GSM3427012", "GSM3427013", "GSM3427014"),
  "Homo sapiens", "candidate dataset after program freeze; final external validation",
  "xenograft_resistant_vs_parental")

# GSE59476: day-28 bevacizumab endpoint versus matched vehicle endpoint.
m <- read_series_matrix(file.path(meta_dir, "GSE59476_series_matrix.txt.gz"))
ids <- m[[1]]
x <- as.matrix(m[, -1]); storage.mode(x) <- "double"
x <- normalizeBetweenArrays(x, method = "quantile")
a <- read_soft_platform(file.path(meta_dir, "GPL13497_full.soft"))
x <- collapse_to_symbol(x, ids, a, "ID", "GENE_SYMBOL")
results[[length(results) + 1L]] <- run_contrast(
  "GSE59476", "day-28 bevacizumab-treated xenografts vs day-28 vehicle xenografts", x,
  c("GSM1437602", "GSM1437603"), c("GSM1437600", "GSM1437601"),
  "Homo sapiens", "candidate dataset after program freeze; design-limited sensitivity",
  "bev_day28_vs_vehicle")

# GSE84048: the same resistant-versus-sensitive state contrast is tested
# separately in bulk tumor and tumor-associated endothelial fractions.
m <- read_series_matrix(file.path(meta_dir, "GSE84048_series_matrix.txt.gz"))
ids <- m[[1]]
x <- as.matrix(m[, -1]); storage.mode(x) <- "double"
x <- normalizeBetweenArrays(log2(x + 1), method = "quantile")
a <- read_soft_platform(file.path(meta_dir, "GPL10787_full.soft"))
x <- collapse_to_symbol(x, ids, a, "ID", "GENE_SYMBOL")
results[[length(results) + 1L]] <- run_contrast(
  "GSE84048", "sunitinib-resistant vs sunitinib-sensitive bulk tumors", x,
  paste0("GSM222649", 6:9), paste0("GSM222649", 2:5),
  "Mus musculus", "candidate dataset after program freeze; final external validation",
  "bulk_resistant_vs_sensitive")
results[[length(results) + 1L]] <- run_contrast(
  "GSE84048", "sunitinib-resistant vs sunitinib-sensitive endothelial fractions", x,
  paste0("GSM222648", 4:7), paste0("GSM222648", 0:3),
  "Mus musculus", "candidate dataset after program freeze; final external validation",
  "endothelium_resistant_vs_sensitive")

# GSE66346: KURC1 is the only model with replicated vehicle and resistant states.
m <- read_series_matrix(file.path(meta_dir, "GSE66346_series_matrix.txt.gz"))
ids <- m[[1]]
x <- as.matrix(m[, -1]); storage.mode(x) <- "double"
x <- normalizeBetweenArrays(x, method = "quantile")
a <- read_geo_annotation(file.path(input_root, "platforms", "GPL6244.annot.gz"))
x <- collapse_to_symbol(x, ids, a, "ID", "Gene symbol")
results[[length(results) + 1L]] <- run_contrast(
  "GSE66346", "KURC1 sunitinib-resistant xenografts vs KURC1 vehicle xenografts", x,
  c("GSM1620060", "GSM1620061"), c("GSM1620057", "GSM1620058"),
  "Homo sapiens", "candidate dataset after program freeze; design-limited sensitivity",
  "KURC1_resistant_vs_vehicle")

all_out <- do.call(rbind, results)
write_tsv(all_out, file.path(out_root, "PHASE8C1_NEW_DATASET_FROZEN_PROGRAM_RESULTS.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat("Phase 8C.1 post-lock sensitivity analyses complete\n")
