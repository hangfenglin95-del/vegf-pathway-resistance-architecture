#!/usr/bin/env Rscript

# Render the figures affected by the expanded external-validation cohort from
# bundled, frozen source data. This does not rerun differential expression or GSEA.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
bundle_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source_root <- file.path(bundle_root, "data", "figure_source_data")
main_out <- file.path(bundle_root, "outputs", "figures")
supp_out <- file.path(bundle_root, "outputs", "figures")
dir.create(main_out, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_out, recursive = TRUE, showWarnings = FALSE)

read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))

blue <- "#0072B2"
orange <- "#E69F00"
green <- "#009E73"
vermillion <- "#D55E00"
purple <- "#CC79A7"
grey <- "#6B7280"
lightgrey <- "#E5E7EB"
dark <- "#1F2937"
heat_cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)

save_dual <- function(stem, width, height, draw_fun, supplementary = FALSE) {
  outdir <- if (supplementary) supp_out else main_out
  png(file.path(outdir, paste0(stem, ".png")), width = width * 600, height = height * 600,
      res = 600, pointsize = 18, bg = "white", type = "cairo")
  draw_fun(); dev.off()
  cairo_pdf(file.path(outdir, paste0(stem, ".pdf")), width = width, height = height,
            pointsize = 18, family = "sans", bg = "white")
  draw_fun(); dev.off()
  svg(file.path(outdir, paste0(stem, ".svg")), width = width, height = height,
      pointsize = 18, family = "sans", bg = "white", onefile = TRUE)
  draw_fun(); dev.off()
}

panel_label <- function(label, x = 0.01, y = 0.98) {
  usr <- par("usr")
  text(usr[1] + x * diff(usr[1:2]), usr[3] + y * diff(usr[3:4]), label,
       adj = c(0, 1), font = 2, cex = 1.12, xpd = NA, col = dark)
}

draw_box <- function(xleft, ybottom, xright, ytop, text_value, fill = "white",
                     border = grey, cex = 0.85, font = 1) {
  rect(xleft, ybottom, xright, ytop, col = fill, border = border, lwd = 1.2)
  text((xleft + xright) / 2, (ybottom + ytop) / 2, text_value,
       cex = cex, font = font, col = dark, family = "sans", xpd = NA)
}

pretty_path <- function(x) {
  labels <- c(
    HALLMARK_MTORC1_SIGNALING = "MTORC1 signaling",
    HALLMARK_UNFOLDED_PROTEIN_RESPONSE = "Unfolded protein response",
    HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-gamma response",
    HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-alpha response",
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "Epithelial-mesenchymal transition",
    HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
  )
  unname(labels[x])
}

short_context <- c(
  GSE64052_786O_sorafenib_resistant_vs_untreated = "64052\n786-O Sor.",
  GSE64052_786O_sunitinib_resistant_vs_untreated = "64052\n786-O Sun.",
  GSE64052_A498_sorafenib_resistant_vs_untreated = "64052\nA498 Sor.",
  GSE86525_HT29_bevacizumab_resistant_vs_control = "86525\nHT-29 Bev.",
  GSE249415_B20_vs_IgG = "249415\nB20",
  GSE121153_xenograft_resistant_vs_parental = "121153\nHuh7",
  GSE84048_bulk_resistant_vs_sensitive = "84048\nBulk",
  GSE84048_endothelium_resistant_vs_sensitive = "84048\nEndo.",
  GSE207976_LSEC_resistant_vs_sensitive = "207976\nLSEC",
  GSE328515_H22_LR_vs_NR = "328515\nH22",
  GSE78698_tumor_long_vs_short_nintedanib = "78698\nTumor",
  GSE78698_endothelium_long_vs_short_nintedanib = "78698\nEndo."
)

# Figure 1 -------------------------------------------------------------------
draw_figure1 <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1)); panel_label("A")
  title("Central question", adj = 0.08, font.main = 2, cex.main = 1.05)
  draw_box(0.08, 0.57, 0.92, 0.82,
           "Does VEGF-pathway resistance\nconverge on one\ntranscriptional endpoint?",
           fill = "#F3F4F6", border = dark, cex = 0.88, font = 2)
  arrows(0.50, 0.56, 0.50, 0.43, length = 0.08, lwd = 1.4, col = grey)
  draw_box(0.08, 0.14, 0.44, 0.39, "Single shared\nendpoint", fill = "#FDECEC", border = vermillion, cex = 0.95)
  draw_box(0.56, 0.14, 0.92, 0.39, "Multiple structured\narchitectures", fill = "#E8F4FA", border = blue, cex = 0.72, font = 2)
  text(0.74, 0.08, "Supported model", col = blue, cex = 0.8, font = 2)

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1)); panel_label("B")
  title("Evidence hierarchy", adj = 0.08, font.main = 2, cex.main = 1.05)
  ys <- c(0.83, 0.66, 0.49, 0.29, 0.09); hs <- c(0.07, 0.075, 0.085, 0.10, 0.075)
  labs <- c("DISCOVERY\n5 datasets / 8 contexts",
            "CROSS-CONTEXT\nPATHWAY ARCHITECTURE",
            "FROZEN RECURRENT /\nDIVERGENT PROGRAMS",
            "FINAL EXTERNAL VALIDATION\n8 datasets / 12 contexts\n+ supportive serial model",
            "HUMAN RESPONSE HISTORY\nGSE79671 · 16 pairs")
  fills <- c("#DCEAF4", "#E8F4FA", "#DFF2EA", "#E8F4FA", "#FCE8D5")
  for (i in seq_along(ys)) {
    draw_box(0.18, ys[i] - hs[i], 0.82, ys[i] + hs[i], labs[i], fill = fills[i], cex = 0.64, font = 2)
    if (i < length(ys)) arrows(0.5, ys[i] - hs[i] - 0.008, 0.5, ys[i + 1] + hs[i + 1] + 0.008, length = 0.04, col = grey)
  }

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1)); panel_label("C")
  title("Dataset roles and independence", adj = 0.08, font.main = 2, cex.main = 1.05)
  draw_box(0.06, 0.68, 0.94, 0.88,
           "DISCOVERY · 5 independent datasets\nGSE76068 · GSE73571 · GSE180687 · GSE64472 · GSE26644",
           fill = "#DCEAF4", border = blue, cex = 0.62, font = 2)
  draw_box(0.06, 0.36, 0.94, 0.64,
           "FINAL EXTERNAL VALIDATION · 8 independent datasets\nGSE64052 · GSE86525 · GSE249415 · GSE121153\nGSE84048 · GSE207976 · GSE328515 · GSE78698",
           fill = "#E8F4FA", border = blue, cex = 0.60, font = 2)
  draw_box(0.06, 0.17, 0.47, 0.31,
           "SUPPORTIVE ONLY\nGSE81465 · GSE45161", fill = "#F3F4F6", border = grey, cex = 0.60, font = 2)
  draw_box(0.53, 0.17, 0.94, 0.31,
           "HUMAN EVIDENCE\nGSE79671 · GSE37138", fill = "#FFF7ED", border = orange, cex = 0.60, font = 2)
  text(0.50, 0.08, "Nested contexts remain grouped within their GEO series.", cex = 0.64, col = grey)

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1)); panel_label("D")
  title("Analysis strategy", adj = 0.08, font.main = 2, cex.main = 1.05)
  xs <- seq(0.1, 0.9, length.out = 5)
  steps <- c("Within-study\nDE", "Ranked\nGSEA", "Pathway\narchitecture", "Frozen-panel\nexternal\nvalidation", "Human\ntranslation")
  for (i in seq_along(xs)) {
    draw_box(xs[i] - 0.075, 0.48, xs[i] + 0.075, 0.70, steps[i], fill = "white", border = blue, cex = 0.49, font = 2)
    if (i < length(xs)) arrows(xs[i] + 0.08, 0.59, xs[i + 1] - 0.08, 0.59, length = 0.045, col = grey)
  }
  text(0.5, 0.30, "Six programs frozen before all external-dataset tests", cex = 0.70, col = dark)
  text(0.5, 0.18, "Dataset-level recurrence primary\ncontext-level results secondary", cex = 0.60, col = vermillion, font = 2)
  mtext("Evidence hierarchy; nested contexts do not increase independent-dataset counts", outer = TRUE,
        side = 1, line = -1, cex = 0.68, col = grey)
}

# Figure 4 -------------------------------------------------------------------
draw_figure4 <- function() {
  layout(matrix(c(1, 1, 2, 3, 4, 4), nrow = 3, byrow = TRUE), heights = c(1.30, 1, 0.84))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)
  val <- read_tab(file.path(source_root, "Figure_4", "Figure4A_ExpandedValidationNES.tsv"))
  val <- val[order(as.numeric(val$pathway_order), as.numeric(val$context_order)), ]
  paths <- unique(val$pathway); contrasts <- unique(val$contrast)
  mat <- matrix(NA_real_, length(paths), length(contrasts), dimnames = list(paths, contrasts))
  for (i in seq_len(nrow(val))) mat[val$pathway[i], val$contrast[i]] <- as.numeric(val$NES[i])

  par(mar = c(7.8, 11.2, 3.2, 1))
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), col = heat_cols, zlim = c(-3.8, 3.8),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(1, at = seq_len(ncol(mat)), labels = unname(short_context[contrasts]), las = 2, tick = FALSE, cex.axis = 0.61, line = -0.3)
  axis(2, at = seq_len(nrow(mat)), labels = pretty_path(paths), las = 2, tick = FALSE, cex.axis = 0.72, line = -0.5)
  box(col = grey, lwd = 0.8)
  abline(h = seq(1.5, nrow(mat) - 0.5, by = 1), col = "white", lwd = 0.4)
  abline(v = seq(1.5, ncol(mat) - 0.5, by = 1), col = "white", lwd = 0.35)
  abline(v = c(3.5, 4.5, 5.5, 7.5, 8.5, 9.5, 10.5), col = dark, lwd = 0.8)
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
    text(j, i, sprintf("%.2f", mat[i, j]), cex = 0.48,
         col = if (abs(mat[i, j]) > 2.2) "white" else dark)
  }
  title("A  Frozen-program NES across 12 formal validation contexts", adj = 0, font.main = 2, cex.main = 1)
  mtext("Blue: negative; red: positive · heavier separators denote independent GEO series", side = 3, adj = 1, cex = 0.60, col = grey, line = 0.2)

  ds <- read_tab(file.path(source_root, "Figure_4", "Figure4BC_RecurrentDatasetSummary.tsv"))
  dataset_order <- unique(ds$dataset)
  plot_dataset <- function(pathway, title_text, positive_col) {
    d <- ds[ds$pathway == pathway, ]; d <- d[match(dataset_order, d$dataset), ]
    y <- rev(seq_len(nrow(d)))
    x <- as.numeric(d$median_NES)
    par(mar = c(4.6, 6.5, 3, 1))
    lim <- max(2.6, max(abs(x)) + 0.55)
    plot(c(-lim, lim), c(0.5, length(y) + 0.5), type = "n", yaxt = "n", xlab = "Median context NES", ylab = "", main = title_text)
    axis(2, at = y, labels = d$dataset, las = 2, tick = FALSE, cex.axis = 0.72)
    abline(v = 0, col = grey, lty = 2)
    segments(0, y, x, y, col = lightgrey, lwd = 3)
    cols <- ifelse(d$dataset_direction == "positive", positive_col, ifelse(d$dataset_direction == "negative", vermillion, purple))
    pch <- ifelse(d$dataset_direction == "mixed", 23, 21)
    points(x, y, pch = pch, bg = cols, col = "white", cex = 1.35)
    text(x, y, sprintf("%.2f", x), pos = ifelse(x >= 0, 4, 2), cex = 0.61)
    legend("bottomright", legend = c("positive", "negative", "mixed"), pt.bg = c(positive_col, vermillion, purple),
           pch = c(21, 21, 23), bty = "n", cex = 0.59)
  }
  plot_dataset("HALLMARK_MTORC1_SIGNALING", "B  MTORC1: 6 positive / 2 negative datasets", blue)
  plot_dataset("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "C  UPR: 5 positive / 2 negative / 1 mixed", green)

  temp <- read_tab(file.path(source_root, "Figure_4", "Figure4D_GSE81465Temporal.tsv"))
  stages <- c("IgG", "G1", "G4", "G9")
  par(mar = c(5.8, 5, 3, 1))
  plot(c(1, 4), range(as.numeric(temp$program_score_mean), na.rm = TRUE) * 1.35,
       type = "n", xaxt = "n", xlab = "", ylab = "Array-level mean rank score", main = "D  Supportive GSE81465 serial states")
  axis(1, at = 1:4, labels = c("IgG", "G1 early", "G4 intermediate", "G9 late"), cex.axis = 0.72)
  for (j in seq_along(c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"))) {
    path <- c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")[j]
    d <- temp[temp$pathway == path, ]; d <- d[match(stages, d$stage), ]
    x <- 1:4 + c(-0.03, 0.03)[j]; y <- as.numeric(d$program_score_mean)
    lines(x, y, col = c(blue, green)[j], lwd = 2)
    points(x, y, pch = 21, bg = c(blue, green)[j], col = "white", cex = 1.15)
  }
  abline(h = 0, col = lightgrey)
  legend("topleft", legend = c("MTORC1", "UPR"), col = c(blue, green), lwd = 2, pch = 19, bty = "n", cex = 0.68)
  mtext("Supportive only: replicate independence unresolved", side = 1, line = 4.1, cex = 0.60, col = grey)
}

# Supplementary Figure S3 ----------------------------------------------------
draw_s3 <- function() {
  d <- read_tab(file.path(source_root, "Figure_4", "Figure4A_ExpandedValidationNES.tsv"))
  d <- d[order(as.numeric(d$pathway_order), as.numeric(d$context_order)), ]
  paths <- unique(d$pathway); contrasts <- unique(d$contrast)
  mat <- matrix(NA_real_, length(paths), length(contrasts), dimnames = list(paths, contrasts))
  fdr <- matrix(NA_real_, length(paths), length(contrasts), dimnames = list(paths, contrasts))
  for (i in seq_len(nrow(d))) {
    mat[d$pathway[i], d$contrast[i]] <- -log10(max(as.numeric(d$padj[i]), 1e-300))
    fdr[d$pathway[i], d$contrast[i]] <- as.numeric(d$padj[i])
  }
  par(mar = c(8.5, 11.5, 3.5, 1), family = "sans")
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), col = colorRampPalette(c("white", "#FDD0A2", vermillion))(101),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(1, at = seq_len(ncol(mat)), labels = unname(short_context[contrasts]), las = 2, tick = FALSE, cex.axis = 0.63)
  axis(2, at = seq_len(nrow(mat)), labels = pretty_path(paths), las = 2, tick = FALSE, cex.axis = 0.72)
  abline(v = c(3.5, 4.5, 5.5, 7.5, 8.5, 9.5, 10.5), col = dark, lwd = 0.8)
  fmt_fdr <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) text(j, i, fmt_fdr(fdr[i, j]), cex = 0.53)
  box(col = grey)
  title("Formal external-validation FDR values", adj = 0, font.main = 2)
  mtext("Cells show frozen Benjamini-Hochberg FDR; heavier separators denote independent GEO series", side = 3, adj = 1, cex = 0.64, col = grey)
}

save_dual("Figure1", 10.8, 8.2, draw_figure1)
save_dual("Figure4", 11.5, 10.8, draw_figure4)
save_dual("FigureS3_Validation_FDR_matrix", 11.5, 7.2, draw_s3, supplementary = TRUE)
cat("Rendered Phase 10B Figure 1, Figure 4, and Supplementary Figure S3.\n")
