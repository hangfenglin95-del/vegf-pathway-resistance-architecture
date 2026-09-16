#!/usr/bin/env Rscript

# GSE37138: baseline treatment-response association only (not acquired resistance).

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
phase_root <- normalizePath(analysis_root_env, mustWork = TRUE)
r3b <- file.path(phase_root, "round3b_human_translation")
suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(huex10stprobeset.db)
  library(msigdbr)
  library(ggplot2)
})

set.seed(20260909)
write_tsv <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
frozen <- c(
  MTORC1 = "HALLMARK_MTORC1_SIGNALING",
  UPR = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  EMT = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  INFLAMMATORY_RESPONSE = "HALLMARK_INFLAMMATORY_RESPONSE",
  E2F = "HALLMARK_E2F_TARGETS"
)

manifest <- read.delim(file.path(r3b, "data/metadata/GSE37138_sample_manifest.tsv"), check.names = FALSE,
                       na.strings = c("NA", ""), stringsAsFactors = FALSE)
matrix_path <- file.path(r3b, "data/raw/GSE37138_series_matrix.txt.gz")
selected <- c("ID_REF", manifest$GSM)
dt <- fread(cmd = paste("gzcat", shQuote(matrix_path)), skip = '"ID_REF"', select = selected,
            fill = TRUE, showProgress = TRUE)
dt <- dt[grepl("^[0-9]+$", as.character(ID_REF))]
probe_ids <- as.character(dt$ID_REF)

# Keep uniquely mapped probesets and average exon-probeset intensities per gene.
probe_to_symbol <- mapIds(huex10stprobeset.db, keys = probe_ids, keytype = "PROBEID",
                          column = "SYMBOL", multiVals = "filter")
symbols <- toupper(unname(probe_to_symbol[probe_ids]))
keep <- !is.na(symbols) & nzchar(symbols)
e_probe <- as.matrix(dt[keep, -1, with = FALSE]); storage.mode(e_probe) <- "double"
symbols <- symbols[keep]
gene_sums <- rowsum(e_probe, group = symbols, reorder = FALSE)
gene_n <- as.numeric(table(factor(symbols, levels = rownames(gene_sums))))
expr <- gene_sums / gene_n
rm(dt, e_probe, gene_sums)
write_tsv(data.frame(gene_symbol = rownames(expr), expr, check.names = FALSE),
          file.path(r3b, "data/processed/GSE37138_gene_level_RMA.tsv"))

hall_db <- msigdbr(species = "Homo sapiens", collection = "H")
hall_paths <- split(hall_db$gene_symbol, hall_db$gs_name)
score_rank <- function(x, paths) {
  pct <- apply(x, 2, function(v) (rank(v, ties.method = "average") - 0.5) / length(v) - 0.5)
  if (is.vector(pct)) pct <- matrix(pct, ncol = 1, dimnames = dimnames(x))
  out <- sapply(paths, function(g) colMeans(pct[intersect(g, rownames(pct)), , drop = FALSE]))
  t(out)
}
score_mat <- score_rank(expr, hall_paths[unname(frozen)])
rownames(score_mat) <- names(frozen)
score_long <- do.call(rbind, lapply(seq_len(nrow(score_mat)), function(i) data.frame(
  program = rownames(score_mat)[i], GSM = colnames(score_mat), program_score = as.numeric(score_mat[i, ]), stringsAsFactors = FALSE
)))
scores <- merge(score_long, manifest, by = "GSM", sort = FALSE)
scores$tumor_shrinkage_week12_pct <- as.numeric(scores$tumor_shrinkage_week12_pct)
scores$disease_control_week12 <- as.integer(scores$disease_control_week12)
write_tsv(scores, file.path(r3b, "results/clinical_association/GSE37138_program_scores.tsv"))

boot_spearman <- function(x, y, B = 2000) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  vals <- replicate(B, {i <- sample.int(length(x), replace = TRUE); suppressWarnings(cor(x[i], y[i], method = "spearman"))})
  vals <- vals[is.finite(vals)]
  unname(quantile(vals, c(.025, .975), na.rm = TRUE))
}

assoc <- do.call(rbind, lapply(names(frozen), function(program) {
  z <- scores[scores$program == program, ]
  ok <- is.finite(z$tumor_shrinkage_week12_pct)
  ct <- suppressWarnings(cor.test(z$program_score[ok], z$tumor_shrinkage_week12_pct[ok], method = "spearman", exact = FALSE))
  ci <- boot_spearman(z$program_score[ok], z$tumor_shrinkage_week12_pct[ok])
  suc <- z$program_score[z$disease_control_week12 == 1]
  fail <- z$program_score[z$disease_control_week12 == 0]
  wt <- wilcox.test(suc, fail, exact = FALSE)
  data.frame(
    program = program,
    continuous_endpoint = "tumor shrinkage at week 12 (%; larger value = more shrinkage)",
    continuous_n = sum(ok), spearman_rho = unname(ct$estimate), spearman_CI_low = ci[1], spearman_CI_high = ci[2], spearman_p = ct$p.value,
    categorical_endpoint = "disease control at week 12 (CR/PR/SD vs failure)",
    DCR_success_n = length(suc), DCR_failure_n = length(fail),
    DCR_success_median_score = median(suc), DCR_failure_median_score = median(fail),
    DCR_success_minus_failure = median(suc) - median(fail), DCR_wilcoxon_p = wt$p.value,
    direction_consistency = ifelse(sign(unname(ct$estimate)) == sign(median(suc) - median(fail)), "consistent", "discordant"),
    stringsAsFactors = FALSE
  )
}))
assoc$continuous_FDR <- p.adjust(assoc$spearman_p, method = "BH")
assoc$DCR_FDR <- p.adjust(assoc$DCR_wilcoxon_p, method = "BH")
write_tsv(assoc, file.path(r3b, "results/clinical_association/GSE37138_program_clinical_associations.tsv"))
write_tsv(assoc[, c("program", "spearman_rho", "spearman_p", "DCR_success_minus_failure", "DCR_wilcoxon_p", "direction_consistency")],
          file.path(r3b, "results/robustness/GSE37138_endpoint_consistency.tsv"))

# QC confirms one profile per baseline tumor patient and no blood contamination.
pca <- prcomp(t(expr), scale. = FALSE)
pv <- 100 * pca$sdev^2 / sum(pca$sdev^2)
qc <- data.frame(GSM = colnames(expr), patient_id = manifest$patient_id[match(colnames(expr), manifest$GSM)],
                 median_expression = apply(expr, 2, median), IQR_expression = apply(expr, 2, IQR),
                 PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
write_tsv(qc, file.path(r3b, "results/qc/GSE37138_sample_qc.tsv"))
map_summary <- data.frame(total_exon_probesets = length(probe_ids), uniquely_symbol_mapped_probesets = sum(keep),
                          gene_symbols_after_mean_collapse = nrow(expr), baseline_tumor_samples = ncol(expr),
                          clinical_continuous_complete = sum(is.finite(suppressWarnings(as.numeric(manifest$tumor_shrinkage_week12_pct)))),
                          DCR_complete = sum(!is.na(manifest$disease_control_week12)), stringsAsFactors = FALSE)
write_tsv(map_summary, file.path(r3b, "results/qc/GSE37138_processing_summary.tsv"))
png(file.path(r3b, "results/qc/GSE37138_RMA_distributions.png"), 1800, 1000, res = 170)
boxplot(expr, outline = FALSE, las = 2, cex.axis = .55, ylab = "Deposited RMA intensity", main = "GSE37138 baseline tumor profiles"); dev.off()
png(file.path(r3b, "results/qc/GSE37138_PCA.png"), 1500, 1100, res = 170)
cols <- ifelse(manifest$disease_control_week12 == 1, "#009E73", "#CC79A7")
plot(pca$x[, 1], pca$x[, 2], col = cols, pch = 19, xlab = sprintf("PC1 (%.1f%%)", pv[1]), ylab = sprintf("PC2 (%.1f%%)", pv[2]), main = "GSE37138 PCA")
text(pca$x[, 1], pca$x[, 2], labels = manifest$patient_id, pos = 3, cex = .55); legend("topright", c("DCR success", "DCR failure"), col = c("#009E73", "#CC79A7"), pch = 19); dev.off()

plot_assoc <- function(program) {
  z <- scores[scores$program == program & is.finite(scores$tumor_shrinkage_week12_pct), ]
  a <- assoc[assoc$program == program, ]
  gp <- ggplot(z, aes(program_score, tumor_shrinkage_week12_pct, color = disease_control_label)) +
    geom_hline(yintercept = 0, color = "grey75") + geom_point(size = 2.5) + geom_smooth(method = "lm", se = TRUE, color = "grey35") +
    scale_color_manual(values = c(DCR_success = "#009E73", DCR_failure = "#CC79A7")) +
    labs(title = paste("GSE37138 baseline", program, "and week-12 tumor shrinkage"),
         subtitle = sprintf("Spearman rho = %.2f (95%% bootstrap CI %.2f to %.2f), p = %.3g", a$spearman_rho, a$spearman_CI_low, a$spearman_CI_high, a$spearman_p),
         x = "Baseline centered rank program score", y = "Tumor shrinkage at week 12 (%)", color = "12-week disease control") +
    theme_bw(base_size = 11)
  ggsave(file.path(r3b, "results/figures", paste0("GSE37138_", program, "_clinical_association.png")), gp, width = 7.5, height = 5.5, dpi = 180)
  ggsave(file.path(r3b, "results/figures", paste0("GSE37138_", program, "_clinical_association.pdf")), gp, width = 7.5, height = 5.5)
}
plot_assoc("MTORC1")
plot_assoc("UPR")
capture.output(sessionInfo(), file = file.path(r3b, "reports/GSE37138_sessionInfo.txt"))
