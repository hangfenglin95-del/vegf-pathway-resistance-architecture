#!/usr/bin/env Rscript

# Candidate-dataset analyses identified after the discovery/pathway freeze.
# Only the six frozen programs are tested. Final evidence roles are assigned
# by config/candidate_dataset_disposition.tsv.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
bundle_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
analysis_root <- Sys.getenv("VEGF_ANALYSIS_ROOT", unset = bundle_root)
input_root <- Sys.getenv("VEGF_CANDIDATE_INPUT_ROOT", unset = file.path(analysis_root, "data", "raw", "candidate_validation"))
legacy_library <- Sys.getenv("VEGF_R_LIBRARY", unset = "")
if (nzchar(legacy_library)) .libPaths(c(legacy_library, .libPaths()))
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

meta_dir <- file.path(input_root, "phase8c2")
out_root <- Sys.getenv("VEGF_OUTPUT_ROOT", unset = file.path(bundle_root, "outputs", "candidate_validation", "phase8c2"))
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
             stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE)
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

collapse_expression <- function(expr, symbols) {
  symbols <- clean_symbol(symbols)
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

collapse_by_annotation <- function(expr, ids, annotation, id_col, symbol_col) {
  idx <- match(as.character(ids), as.character(annotation[[id_col]]))
  collapse_expression(expr, annotation[[symbol_col]][idx])
}

collapse_counts <- function(counts, symbols) {
  symbols <- clean_symbol(symbols)
  keep <- !is.na(symbols) & apply(counts, 1, function(v) all(is.finite(v) & v >= 0))
  counts <- counts[keep, , drop = FALSE]
  symbols <- symbols[keep]
  rowsum(counts, group = symbols, reorder = FALSE)
}

frozen <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

mouse_hallmark <- suppressWarnings(msigdbr(species = "Mus musculus", collection = "H"))
mouse_hallmark$gene_symbol <- toupper(mouse_hallmark$gene_symbol)
pathways <- split(mouse_hallmark$gene_symbol, mouse_hallmark$gs_name)

run_fgsea <- function(dataset, context_id, contrast_label, fit, numerator,
                      denominator, evidence_role, input_type) {
  tt <- topTable(fit, number = Inf, sort.by = "none")
  tt$gene_symbol <- rownames(tt)
  tt <- tt[, c("gene_symbol", setdiff(colnames(tt), "gene_symbol"))]
  ranks <- fit$t[, 1]
  names(ranks) <- rownames(fit$t)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(81465)
  fg <- fgseaMultilevel(pathways = pathways, stats = ranks,
                        minSize = 15, maxSize = 500, eps = 0)
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";",
                           FUN.VALUE = character(1))
  out <- fg[match(frozen, fg$pathway),
            c("pathway", "NES", "pval", "padj", "size", "leadingEdge")]
  stopifnot(nrow(out) == 6L, !any(is.na(out$pathway)))
  out$dataset <- dataset
  out$context_id <- context_id
  out$contrast <- contrast_label
  out$n_numerator <- length(numerator)
  out$n_denominator <- length(denominator)
  out$direction <- ifelse(out$NES > 0, "positive", "negative")
  out$evidence_role <- evidence_role
  out$input_type <- input_type
  out$identified_date <- "2026-09-12"
  out <- out[, c("dataset", "context_id", "contrast", "pathway", "NES", "pval",
                 "padj", "size", "direction", "n_numerator", "n_denominator",
                 "evidence_role", "input_type", "identified_date", "leadingEdge")]
  d <- file.path(out_root, dataset)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  safe <- gsub("[^A-Za-z0-9]+", "_", context_id)
  write_tsv(tt, file.path(d, paste0(dataset, "_", safe, "_gene_statistics.tsv")))
  write_tsv(data.frame(sample = c(denominator, numerator),
                       group = c(rep("denominator", length(denominator)),
                                 rep("numerator", length(numerator)))),
            file.path(d, paste0(dataset, "_", safe, "_sample_manifest.tsv")))
  out
}

run_expression_contrast <- function(dataset, context_id, contrast_label, expr,
                                    numerator, denominator, evidence_role,
                                    input_type) {
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
  run_fgsea(dataset, context_id, contrast_label, fit, numerator, denominator,
            evidence_role, input_type)
}

run_count_contrast <- function(dataset, context_id, contrast_label, counts,
                               numerator, denominator, evidence_role) {
  samples <- c(denominator, numerator)
  x <- round(counts[, samples, drop = FALSE])
  group <- factor(c(rep("denominator", length(denominator)),
                    rep("numerator", length(numerator))),
                  levels = c("denominator", "numerator"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  dge <- DGEList(x)
  keep <- filterByExpr(dge, design)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(contrasts.fit(lmFit(v, design),
                              makeContrasts(numerator - denominator, levels = design)),
                robust = TRUE)
  run_fgsea(dataset, context_id, contrast_label, fit, numerator, denominator,
            evidence_role, "gene-level counts; TMM/voom/limma")
}

results <- list()

# GSE207976: liver sinusoidal endothelial cells isolated from matched
# bevacizumab-sensitive and -resistant HCT116 CRCLM xenografts (3 vs 3).
d207 <- read.delim(gzfile(file.path(meta_dir,
  "GSE207976_LSEC-RE_vs_LSEC-CTL_Allgene_info.txt.gz")), check.names = FALSE)
c207_names <- grep("^Count_", names(d207), value = TRUE)
c207 <- as.matrix(d207[, c207_names]); storage.mode(c207) <- "double"
colnames(c207) <- sub("^Count_[[:space:]]*", "", c207_names)
c207 <- collapse_counts(c207, d207$symbol)
results[[length(results) + 1L]] <- run_count_contrast(
  "GSE207976", "LSEC_resistant_vs_sensitive",
  "LSECs from bevacizumab-resistant vs bevacizumab-sensitive CRCLM xenografts",
  c207, paste0("LSEC-RE", 1:3), paste0("LSEC-C", 1:3),
  "candidate dataset after program freeze; final external validation; stromal endothelial compartment")

# GSE221557: the archive is labelled non-coding RNA profiling, but the three
# deposited transcript-level tables jointly reconstruct six sample-level FPKM
# profiles with broad protein-coding Hallmark coverage.
pair_files <- file.path(meta_dir, paste0("GSE221557_Re", 1:3,
  "_vs_CTL", 1:3, ".lncRNA.transcript_level.Differential_analysis_results.txt.gz"))
pairs <- lapply(seq_along(pair_files), function(i) {
  d <- read.delim(gzfile(pair_files[i]), check.names = FALSE)
  d[, c("transcript_id", "gene_id", "gene_name", paste0("Re", i, "_FPKM"),
        paste0("CTL", i, "_FPKM"))]
})
d221 <- Reduce(function(x, y) merge(x, y, by = c("transcript_id", "gene_id", "gene_name"), all = TRUE), pairs)
x221 <- as.matrix(d221[, unlist(lapply(1:3, function(i) c(paste0("CTL", i, "_FPKM"), paste0("Re", i, "_FPKM"))))])
storage.mode(x221) <- "double"; x221[is.na(x221)] <- 0
x221 <- rowsum(x221, group = clean_symbol(d221$gene_name), reorder = FALSE, na.rm = TRUE)
x221 <- normalizeBetweenArrays(log2(x221 + 0.1), method = "quantile")
results[[length(results) + 1L]] <- run_expression_contrast(
  "GSE221557", "xenograft_resistant_vs_sensitive",
  "bevacizumab-resistant vs bevacizumab-sensitive HCT116 CRCLM xenografts",
  x221, paste0("Re", 1:3, "_FPKM"), paste0("CTL", 1:3, "_FPKM"),
  "candidate dataset after program freeze; design-limited sensitivity with archive-label caveat",
  "transcript-level FPKM aggregated to gene symbol; log2(FPKM+0.1); quantile normalized; limma")

# GSE132568: in-vivo-derived 4T1 axitinib-resistant cells maintained ex vivo;
# the latter is retained as an explicit limitation.
d132 <- read.delim(gzfile(file.path(meta_dir,
  "GSE132568_m4T1-RSEM-estimated-counts.txt.gz")), check.names = FALSE)
map207 <- setNames(clean_symbol(d207$symbol), sub("\\..*$", "", d207$gene_id))
sym132 <- unname(map207[sub("\\..*$", "", d132$geneID)])
c132 <- as.matrix(d132[, -1]); storage.mode(c132) <- "double"
c132 <- collapse_counts(c132, sym132)
results[[length(results) + 1L]] <- run_count_contrast(
  "GSE132568", "AxR_vs_parental",
  "in-vivo-derived axitinib-resistant 4T1 cells vs parental 4T1 cells",
  c132, paste0("m4T1.AxR.", 1:3), paste0("m4T1.P.", 1:3),
  "candidate dataset after program freeze; design-limited sensitivity with ex-vivo-maintenance limitation")

# GSE328515: LibreOffice is used only as a read-only converter for the binary
# legacy .xls deposited by GEO; the original gzip is never modified.
xls_tmp <- tempfile(fileext = ".xls")
csv_dir <- tempfile(pattern = "gse328515_")
dir.create(csv_dir)
in_con <- gzfile(file.path(meta_dir, "GSE328515_Processed-data_counts.xls.gz"), "rb")
out_con <- file(xls_tmp, "wb")
repeat {
  buf <- readBin(in_con, what = "raw", n = 1024 * 1024)
  if (!length(buf)) break
  writeBin(buf, out_con)
}
close(in_con); close(out_con)
soffice <- Sys.which("soffice")
if (!nzchar(soffice)) stop("LibreOffice 'soffice' was not found on PATH; install LibreOffice or add soffice to PATH.")
status <- system2(soffice,
  c("--headless", "--convert-to", "csv", "--outdir", shQuote(csv_dir), shQuote(xls_tmp)),
  stdout = TRUE, stderr = TRUE)
csv_file <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(csv_file) == 1L)
d328 <- read.csv(csv_file, check.names = FALSE)
c328 <- as.matrix(d328[, grep("_count$", names(d328), value = TRUE)]); storage.mode(c328) <- "double"
c328 <- collapse_counts(c328, d328$gene_name)
results[[length(results) + 1L]] <- run_count_contrast(
  "GSE328515", "H22_LR_vs_NR",
  "acquired lenvatinib-resistant vs lenvatinib-sensitive H22 tumors",
  c328, paste0("LR", 1:3, "_count"), paste0("NR", 1:3, "_count"),
  "candidate dataset after program freeze; final external validation; independent in-vivo HCC model")

# GSE78698: long-term escape versus short-term response, separated by the
# sorted tumor and tumor-endothelial compartments.
m786 <- read_series_matrix(file.path(meta_dir, "GSE78698_series_matrix.txt.gz"))
ids786 <- m786[[1]]
x786 <- as.matrix(m786[, -1]); storage.mode(x786) <- "double"
x786 <- normalizeBetweenArrays(x786, method = "quantile")
a786 <- read_geo_annotation(file.path(meta_dir, "GPL6246.annot.gz"))
x786 <- collapse_by_annotation(x786, ids786, a786, "ID", "Gene symbol")
results[[length(results) + 1L]] <- run_expression_contrast(
  "GSE78698", "tumor_long_vs_short_nintedanib",
  "long-term vs short-term nintedanib-exposed sorted tumor cells",
  x786, paste0("GSM207279", 2:6),
  c(paste0("GSM207278", 8:9), paste0("GSM207279", 0:1)),
  "candidate dataset after program freeze; final external validation; sorted tumor compartment",
  "deposited normalized microarray matrix; within-study quantile normalization; limma")
results[[length(results) + 1L]] <- run_expression_contrast(
  "GSE78698", "endothelium_long_vs_short_nintedanib",
  "long-term vs short-term nintedanib-exposed sorted tumor endothelial cells",
  x786, paste0("GSM207280", 3:7), paste0("GSM207280", 0:2),
  "candidate dataset after program freeze; final external validation; sorted endothelial compartment",
  "deposited normalized microarray matrix; within-study quantile normalization; limma")

# GSE80778: two resistance-associated endpoint contrasts against the matched
# vehicle endpoint. These are supportive because treatment and escape effects
# cannot be completely separated.
m807 <- read_series_matrix(file.path(meta_dir, "GSE80778_series_matrix.txt.gz"))
ids807 <- m807[[1]]
x807 <- as.matrix(m807[, -1]); storage.mode(x807) <- "double"
x807 <- normalizeBetweenArrays(x807, method = "quantile")
a807 <- read_soft_platform(file.path(input_root, "phase8c1", "GPL10787_full.soft"))
x807 <- collapse_by_annotation(x807, ids807, a807, "ID", "GENE_SYMBOL")
veh807 <- paste0("GSM21369", 29:31)
results[[length(results) + 1L]] <- run_expression_contrast(
  "GSE80778", "B20_endpoint_vs_vehicle",
  "B20 anti-VEGF resistance-associated endpoint vs vehicle endpoint",
  x807, paste0("GSM21369", 17:20), veh807,
  "candidate dataset after program freeze; design-limited sensitivity; treatment/escape conflation",
  "deposited normalized microarray matrix; within-study quantile normalization; limma")
results[[length(results) + 1L]] <- run_expression_contrast(
  "GSE80778", "nintedanib_endpoint_vs_vehicle",
  "nintedanib resistance-associated endpoint vs vehicle endpoint",
  x807, paste0("GSM21369", 21:24), veh807,
  "candidate dataset after program freeze; design-limited sensitivity; treatment/escape conflation",
  "deposited normalized microarray matrix; within-study quantile normalization; limma")

all_out <- do.call(rbind, results)
stopifnot(length(unique(all_out$dataset)) == 6L,
          length(unique(all_out$context_id)) == 8L,
          nrow(all_out) == 48L)
write_tsv(all_out, file.path(out_root, "PHASE8C2_NEW_DATASET_FROZEN_PROGRAM_RESULTS.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat("Phase 8C.2 post-lock analyses complete:", length(unique(all_out$dataset)),
    "datasets,", length(unique(all_out$context_id)), "contexts,", nrow(all_out), "tests\n")
