#!/usr/bin/env Rscript

# GSE79671: paired human longitudinal analysis of frozen Round 2B programs.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
phase_root <- normalizePath(analysis_root_env, mustWork = TRUE)
r3b <- file.path(phase_root, "round3b_human_translation")
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
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
roles <- c(
  MTORC1 = "PRIMARY_SHARED_CANDIDATE", UPR = "SECONDARY_SHARED_CANDIDATE",
  IFNG = "PRIMARY_CONTEXT_DEPENDENT", IFNA = "PRIMARY_CONTEXT_DEPENDENT",
  EMT = "PRIMARY_CONTEXT_DEPENDENT", INFLAMMATORY_RESPONSE = "PRIMARY_CONTEXT_DEPENDENT",
  E2F = "CAUTIONARY_DIVERGENT_CANDIDATE"
)

raw_path <- file.path(r3b, "data/raw/GSE79671_CountMatrix.txt.gz")
manifest <- read.delim(file.path(r3b, "data/metadata/GSE79671_sample_manifest.tsv"), check.names = FALSE)
dt <- fread(cmd = paste("gzcat", shQuote(raw_path)), showProgress = FALSE)
setnames(dt, 1, "ensembl_gene")
count_titles <- setdiff(names(dt), "ensembl_gene")
manifest <- manifest[match(count_titles, manifest$sample_title), ]
stopifnot(!anyNA(manifest$GSM), identical(manifest$sample_title, count_titles))

symbols <- mapIds(org.Hs.eg.db, keys = dt$ensembl_gene, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep_map <- !is.na(symbols) & nzchar(symbols)
counts_probe <- as.matrix(dt[keep_map, -1, with = FALSE])
storage.mode(counts_probe) <- "integer"
counts <- rowsum(counts_probe, group = toupper(unname(symbols[keep_map])), reorder = FALSE)
rm(dt, counts_probe)

# A single frozen expression filter is used for all four contrasts.
y_all <- DGEList(counts = counts)
y_all <- calcNormFactors(y_all, method = "TMM")
global_keep <- rowSums(cpm(y_all) >= 1) >= 5
y_all <- y_all[global_keep, , keep.lib.sizes = FALSE]
y_all <- calcNormFactors(y_all, method = "TMM")
logcpm <- cpm(y_all, log = TRUE, prior.count = 0.5)
write_tsv(data.frame(gene_symbol = rownames(logcpm), logcpm, check.names = FALSE),
          file.path(r3b, "data/processed/GSE79671_TMM_logCPM.tsv"))

fit_voom <- function(sample_idx, design, contrast_vector, contrast_id) {
  stopifnot(qr(design)$rank == ncol(design))
  y <- DGEList(counts = y_all$counts[, sample_idx, drop = FALSE])
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrasts = matrix(contrast_vector, ncol = 1,
                                               dimnames = list(colnames(design), contrast_id)))
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  out <- data.frame(gene_symbol = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr,
                    t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, B = tt$B,
                    contrast = contrast_id, stringsAsFactors = FALSE)
  out <- out[order(out$t, decreasing = TRUE), ]
  write_tsv(out, file.path(r3b, "results/differential", paste0(contrast_id, "_gene_statistics.tsv")))
  list(table = out, ranks = sort(setNames(out$t, out$gene_symbol), decreasing = TRUE), voom = v)
}

paired <- manifest$paired_status == "complete_pair"

# H1 and H2: subject-paired models (~ subject + time), fitted separately.
fits <- list()
for (grp in c("Responder", "NonResponder")) {
  idx <- which(paired & manifest$responder_status == grp)
  md <- droplevels(manifest[idx, ])
  md$subject_id <- factor(md$subject_id)
  md$before_after <- factor(md$before_after, levels = c("Pre", "Post"))
  design <- model.matrix(~ 0 + subject_id + before_after, md)
  id <- if (grp == "Responder") "H1_Responder_Post_vs_Pre" else "H2_NonResponder_Post_vs_Pre"
  write_tsv(data.frame(sample = md$sample_title, design, check.names = FALSE),
            file.path(r3b, "results/differential", paste0(id, "_design_matrix.tsv")))
  cv <- setNames(rep(0, ncol(design)), colnames(design)); cv["before_afterPost"] <- 1
  fits[[id]] <- fit_voom(idx, design, cv, id)
}

# H3: full-rank patient-fixed-effect model with group-specific post slopes.
idx3 <- which(paired)
md3 <- droplevels(manifest[idx3, ])
md3$subject_id <- factor(md3$subject_id)
patient_design <- model.matrix(~ 0 + subject_id, md3)
design3 <- cbind(
  patient_design,
  NonResponder_Post = as.integer(md3$responder_status == "NonResponder" & md3$before_after == "Post"),
  Responder_Post = as.integer(md3$responder_status == "Responder" & md3$before_after == "Post")
)
stopifnot(qr(design3)$rank == ncol(design3))
write_tsv(data.frame(sample = md3$sample_title, design3, check.names = FALSE),
          file.path(r3b, "results/differential/H3_Responder_change_vs_NonResponder_change_design_matrix.tsv"))
cv3 <- setNames(rep(0, ncol(design3)), colnames(design3))
cv3[c("Responder_Post", "NonResponder_Post")] <- c(1, -1)
fits$H3_Responder_change_vs_NonResponder_change <- fit_voom(
  idx3, design3, cv3, "H3_Responder_change_vs_NonResponder_change"
)

# H4: all available pretreatment samples, with batch adjustment when estimable.
idx4 <- which(manifest$before_after == "Pre")
md4 <- droplevels(manifest[idx4, ])
md4$batch <- factor(md4$batch)
md4$responder_status <- factor(md4$responder_status, levels = c("NonResponder", "Responder"))
design4 <- model.matrix(~ batch + responder_status, md4)
stopifnot(qr(design4)$rank == ncol(design4))
write_tsv(data.frame(sample = md4$sample_title, design4, check.names = FALSE),
          file.path(r3b, "results/differential/H4_Responder_Pre_vs_NonResponder_Pre_design_matrix.tsv"))
cv4 <- setNames(rep(0, ncol(design4)), colnames(design4)); cv4["responder_statusResponder"] <- 1
fits$H4_Responder_Pre_vs_NonResponder_Pre <- fit_voom(
  idx4, design4, cv4, "H4_Responder_Pre_vs_NonResponder_Pre"
)

# Full ranked Hallmark and Reactome GSEA for H1-H4.
hall_db <- msigdbr(species = "Homo sapiens", collection = "H")
hall_paths <- split(hall_db$gene_symbol, hall_db$gs_name)
react_db <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
react_paths <- split(react_db$gene_symbol, react_db$gs_name)
gsea_tables <- list()
for (id in names(fits)) {
  ranks <- fits[[id]]$ranks
  set.seed(20260909)
  h <- as.data.frame(fgseaMultilevel(hall_paths, ranks, minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000))
  h <- h[order(h$padj, -abs(h$NES)), c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")]
  h$leadingEdge <- vapply(h$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  write_tsv(h, file.path(r3b, "results/gsea", paste0(id, "_hallmark_full.tsv")))
  set.seed(20260909)
  rr <- as.data.frame(fgseaMultilevel(react_paths, ranks, minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000))
  rr <- rr[order(rr$padj, -abs(rr$NES)), c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")]
  rr$leadingEdge <- vapply(rr$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  write_tsv(rr, file.path(r3b, "results/gsea", paste0(id, "_reactome_full.tsv")))
  gsea_tables[[id]] <- h
}

predefined <- data.frame(program = names(frozen), pathway = unname(frozen), predefined_role = unname(roles), stringsAsFactors = FALSE)
for (id in names(gsea_tables)) {
  z <- gsea_tables[[id]]
  m <- match(predefined$pathway, z$pathway)
  predefined[[paste0(id, "_NES")]] <- z$NES[m]
  predefined[[paste0(id, "_p")]] <- z$pval[m]
  predefined[[paste0(id, "_FDR")]] <- z$padj[m]
}
predefined$direction <- apply(predefined[grep("_NES$", names(predefined))], 1, function(v) paste(ifelse(v > 0, "positive", "negative"), collapse = ";"))
predefined$interpretation <- ifelse(
  predefined$program %in% c("MTORC1", "UPR"),
  "Judge H1 longitudinal direction first, then H3 and patient-level deltas",
  ifelse(predefined$program == "E2F", "Cautionary only; no post-hoc promotion", "Context dependence: compare H1, H2, H3 and patient deltas")
)
write_tsv(predefined, file.path(r3b, "results/predefined_programs/GSE79671_predefined_hallmark.tsv"))

# Fixed sample-level score: mean within-sample percentile rank, centered on zero.
score_rank <- function(expr, paths) {
  pct <- apply(expr, 2, function(v) (rank(v, ties.method = "average") - 0.5) / length(v) - 0.5)
  if (is.vector(pct)) pct <- matrix(pct, ncol = 1, dimnames = dimnames(expr))
  out <- sapply(paths, function(g) colMeans(pct[intersect(g, rownames(pct)), , drop = FALSE]))
  t(out)
}
score_mat <- score_rank(logcpm, hall_paths[unname(frozen)])
rownames(score_mat) <- names(frozen)
score_long <- do.call(rbind, lapply(seq_len(nrow(score_mat)), function(i) data.frame(
  program = rownames(score_mat)[i], sample_title = colnames(score_mat), score = as.numeric(score_mat[i, ]), stringsAsFactors = FALSE
)))
score_long <- merge(score_long, manifest[, c("GSM", "sample_title", "subject_id", "responder_status", "before_after", "paired_status", "batch")], by = "sample_title", sort = FALSE)
write_tsv(score_long, file.path(r3b, "results/patient_level/GSE79671_patient_program_scores.tsv"))

deltas <- list()
for (program in names(frozen)) {
  z <- score_long[score_long$program == program & score_long$paired_status == "complete_pair", ]
  wide <- reshape(z[, c("subject_id", "responder_status", "batch", "before_after", "score")],
                  idvar = c("subject_id", "responder_status", "batch"), timevar = "before_after", direction = "wide")
  names(wide)[names(wide) == "score.Pre"] <- "pre_score"
  names(wide)[names(wide) == "score.Post"] <- "post_score"
  wide$delta <- wide$post_score - wide$pre_score
  wide$program <- program
  deltas[[program]] <- wide[, c("subject_id", "responder_status", "batch", "program", "pre_score", "post_score", "delta")]
}
deltas <- do.call(rbind, deltas); rownames(deltas) <- NULL
write_tsv(deltas, file.path(r3b, "results/patient_level/GSE79671_patient_program_deltas.tsv"))

delta_stats <- do.call(rbind, lapply(names(frozen), function(program) {
  z <- deltas[deltas$program == program, ]
  r <- z$delta[z$responder_status == "Responder"]
  n <- z$delta[z$responder_status == "NonResponder"]
  data.frame(program = program, responder_n = length(r), responder_median_delta = median(r),
             responder_positive = sum(r > 0), responder_signed_rank_p = wilcox.test(r, mu = 0, exact = FALSE)$p.value,
             nonresponder_n = length(n), nonresponder_median_delta = median(n),
             nonresponder_positive = sum(n > 0), nonresponder_signed_rank_p = wilcox.test(n, mu = 0, exact = FALSE)$p.value,
             responder_minus_nonresponder_median_delta = median(r) - median(n),
             interaction_rank_sum_p = wilcox.test(r, n, exact = FALSE)$p.value,
             stringsAsFactors = FALSE)
}))
write_tsv(delta_stats, file.path(r3b, "results/interaction/GSE79671_program_delta_statistics.tsv"))

# QC tables and required visualizations.
lib_size <- colSums(counts)
expressed <- colSums(counts >= 10)
pca <- prcomp(t(logcpm), scale. = FALSE)
pv <- 100 * pca$sdev^2 / sum(pca$sdev^2)
qc <- data.frame(manifest[, c("GSM", "sample_title", "subject_id", "responder_status", "before_after", "paired_status", "batch")],
                 library_size = lib_size[manifest$sample_title], expressed_gene_count_ge10 = expressed[manifest$sample_title],
                 PC1 = pca$x[manifest$sample_title, 1], PC2 = pca$x[manifest$sample_title, 2], stringsAsFactors = FALSE)
write_tsv(qc, file.path(r3b, "results/qc/GSE79671_sample_qc.tsv"))

png(file.path(r3b, "results/qc/GSE79671_library_size_expressed_genes.png"), 1800, 900, res = 170)
par(mfrow = c(1, 2), mar = c(9, 4, 3, 1)); barplot(lib_size / 1e6, las = 2, cex.names = .55, ylab = "Million fragments", main = "Library size")
barplot(expressed, las = 2, cex.names = .55, ylab = "Genes with count >=10", main = "Expressed genes"); dev.off()
png(file.path(r3b, "results/qc/GSE79671_logCPM_distributions.png"), 1800, 1000, res = 170)
boxplot(logcpm, outline = FALSE, las = 2, cex.axis = .55, ylab = "TMM logCPM", main = "GSE79671 expression distributions"); dev.off()
cols <- ifelse(manifest$responder_status == "Responder", "#D55E00", "#0072B2")
pchs <- ifelse(manifest$before_after == "Post", 17, 16)
png(file.path(r3b, "results/qc/GSE79671_PCA_response_time_batch.png"), 1700, 1100, res = 170)
plot(pca$x[, 1], pca$x[, 2], col = cols, pch = pchs, xlab = sprintf("PC1 (%.1f%%)", pv[1]), ylab = sprintf("PC2 (%.1f%%)", pv[2]), main = "GSE79671 PCA")
text(pca$x[, 1], pca$x[, 2], labels = paste0(manifest$subject_id, "-", manifest$before_after, "-", manifest$batch), pos = 3, cex = .55)
legend("topright", c("Responder Pre", "Responder Post", "NonResponder Pre", "NonResponder Post"), col = c("#D55E00", "#D55E00", "#0072B2", "#0072B2"), pch = c(16, 17, 16, 17), cex = .75); dev.off()
png(file.path(r3b, "results/qc/GSE79671_MDS.png"), 1500, 1100, res = 170)
plotMDS(y_all, col = cols, pch = pchs, labels = paste0(manifest$subject_id, "-", manifest$before_after), main = "GSE79671 MDS"); dev.off()
cormat <- cor(logcpm, method = "spearman")
png(file.path(r3b, "results/qc/GSE79671_sample_correlation.png"), 1800, 1700, res = 180)
heatmap(cormat, scale = "none", col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), margins = c(8, 8), main = "GSE79671 sample Spearman correlation"); dev.off()

plot_trajectory <- function(program) {
  z <- score_long[score_long$program == program & score_long$paired_status == "complete_pair", ]
  z$before_after <- factor(z$before_after, levels = c("Pre", "Post"))
  gp <- ggplot(z, aes(before_after, score, group = subject_id, color = responder_status)) +
    geom_line(alpha = .7) + geom_point(size = 2) + facet_wrap(~ responder_status) +
    scale_color_manual(values = c(Responder = "#D55E00", NonResponder = "#0072B2")) +
    labs(title = paste("GSE79671", program, "patient trajectories"), x = NULL,
         y = "Centered mean-percentile rank score") + theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(r3b, "results/figures", paste0("GSE79671_", program, "_patient_trajectory.png")), gp, width = 8, height = 4.7, dpi = 180)
  ggsave(file.path(r3b, "results/figures", paste0("GSE79671_", program, "_patient_trajectory.pdf")), gp, width = 8, height = 4.7)
}
for (p in c("MTORC1", "UPR", "IFNG", "EMT")) plot_trajectory(p)

gp_delta <- ggplot(deltas, aes(responder_status, delta, color = responder_status)) +
  geom_hline(yintercept = 0, color = "grey60") + geom_boxplot(outlier.shape = NA, width = .5) +
  geom_jitter(width = .12, size = 1.8, alpha = .8) + facet_wrap(~ program, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(Responder = "#D55E00", NonResponder = "#0072B2")) +
  labs(title = "GSE79671 patient-level program changes", x = NULL, y = "Post - Pre score") +
  theme_bw(base_size = 10) + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(r3b, "results/figures/GSE79671_interaction_program_deltas.png"), gp_delta, width = 12, height = 6.5, dpi = 180)
ggsave(file.path(r3b, "results/figures/GSE79671_interaction_program_deltas.pdf"), gp_delta, width = 12, height = 6.5)

# Leave-one-patient-out sensitivity for MTORC1 and UPR.
loo_rows <- list()
for (grp in c("Responder", "NonResponder")) {
  subjects <- unique(manifest$subject_id[paired & manifest$responder_status == grp])
  for (drop_subject in c("none", subjects)) {
    idx <- which(paired & manifest$responder_status == grp & (drop_subject == "none" | manifest$subject_id != drop_subject))
    md <- droplevels(manifest[idx, ]); md$subject_id <- factor(md$subject_id); md$before_after <- factor(md$before_after, levels = c("Pre", "Post"))
    design <- model.matrix(~ 0 + subject_id + before_after, md)
    cv <- setNames(rep(0, ncol(design)), colnames(design)); cv["before_afterPost"] <- 1
    fit <- fit_voom(idx, design, cv, paste0("TEMP_", grp, "_", drop_subject))
    set.seed(20260909)
    fg <- as.data.frame(fgseaMultilevel(hall_paths[unname(frozen[c("MTORC1", "UPR")])], fit$ranks,
                                       minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000))
    for (program in c("MTORC1", "UPR")) {
      zz <- fg[fg$pathway == frozen[program], ]
      loo_rows[[length(loo_rows) + 1L]] <- data.frame(group = grp, excluded_patient = drop_subject,
        n_pairs = length(unique(md$subject_id)), program = program, NES = zz$NES, pval = zz$pval, padj = zz$padj,
        direction = ifelse(zz$NES > 0, "positive", "negative"), stringsAsFactors = FALSE)
    }
    unlink(file.path(r3b, "results/differential", paste0("TEMP_", grp, "_", drop_subject, "_gene_statistics.tsv")))
  }
}
write_tsv(do.call(rbind, loo_rows), file.path(r3b, "results/robustness/GSE79671_leave_one_patient_out.tsv"))

# Batch estimability and transcriptomic-remodeling summary.
batch_tab <- as.data.frame(xtabs(~ batch + responder_status + before_after, manifest))
write_tsv(batch_tab,
          file.path(r3b, "results/qc/GSE79671_batch_response_time_crosstab.tsv"))
design_diag <- data.frame(
  contrast = c("H1", "H2", "H3", "H4"),
  n_samples = c(sum(paired & manifest$responder_status == "Responder"), sum(paired & manifest$responder_status == "NonResponder"), sum(paired), sum(manifest$before_after == "Pre")),
  design_columns = c(7, 11, ncol(design3), ncol(design4)),
  design_rank = c(7, 11, qr(design3)$rank, qr(design4)$rank),
  batch_handling = c("batch constant within pair and aliased with subject; subject fixed effect", "batch constant within pair and aliased with subject; subject fixed effect", "patient fixed effect removes time-invariant batch", "batch included explicitly"),
  stringsAsFactors = FALSE)
write_tsv(design_diag, file.path(r3b, "results/robustness/GSE79671_design_and_batch_diagnostics.tsv"))

remodeling <- do.call(rbind, lapply(c("H1_Responder_Post_vs_Pre", "H2_NonResponder_Post_vs_Pre"), function(id) {
  z <- fits[[id]]$table
  data.frame(contrast = id, genes_tested = nrow(z), FDR_lt_0_05 = sum(z$adj.P.Val < .05),
             nominal_p_lt_0_05 = sum(z$P.Value < .05), median_abs_t = median(abs(z$t)),
             percentile95_abs_logFC = unname(quantile(abs(z$logFC), .95)), stringsAsFactors = FALSE)
}))
write_tsv(remodeling, file.path(r3b, "results/differential/GSE79671_transcriptomic_remodeling_summary.tsv"))

capture.output(sessionInfo(), file = file.path(r3b, "reports/GSE79671_sessionInfo.txt"))
