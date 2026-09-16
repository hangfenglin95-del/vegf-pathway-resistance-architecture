#!/usr/bin/env Rscript

# Deterministic rendering of publication figures from frozen source-data TSVs.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the project root can be resolved.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
bundle_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source_root <- file.path(bundle_root, "data", "figure_source_data")
figure_root <- file.path(bundle_root, "outputs", "figures")
supp_root <- file.path(bundle_root, "outputs", "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_root, recursive = TRUE, showWarnings = FALSE)

render_set <- strsplit(Sys.getenv("VEGF_FIGURES", "all"), ",", fixed = TRUE)[[1]]
should_render <- function(id) "all" %in% render_set || id %in% render_set

read_tab <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
}

blue <- "#0072B2"
orange <- "#E69F00"
green <- "#009E73"
vermillion <- "#D55E00"
purple <- "#CC79A7"
sky <- "#56B4E9"
grey <- "#6B7280"
lightgrey <- "#E5E7EB"
dark <- "#1F2937"
heat_cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)

save_dual <- function(stem, width, height, draw_fun, supplementary = FALSE) {
  outdir <- if (supplementary) supp_root else figure_root
  stem_out <- stem
  png(file.path(outdir, paste0(stem_out, ".png")), width = width * 600, height = height * 600,
      res = 600, pointsize = 18, bg = "white", type = "cairo")
  draw_fun()
  dev.off()
  cairo_pdf(file.path(outdir, paste0(stem_out, ".pdf")), width = width, height = height,
            pointsize = 18, family = "sans", bg = "white")
  draw_fun()
  dev.off()
  svg(file.path(outdir, paste0(stem_out, ".svg")), width = width, height = height,
      pointsize = 18, family = "sans", bg = "white", onefile = TRUE)
  draw_fun()
  dev.off()
}

panel_label <- function(label, x = 0.01, y = 0.98) {
  usr <- par("usr")
  text(usr[1] + x * diff(usr[1:2]), usr[3] + y * diff(usr[3:4]), label,
       adj = c(0, 1), font = 2, cex = 1.2, xpd = NA, col = dark)
}

draw_box <- function(xleft, ybottom, xright, ytop, text_value, fill = "white",
                     border = grey, cex = 0.9, font = 1, line_height = 1) {
  rect(xleft, ybottom, xright, ytop, col = fill, border = border, lwd = 1.2)
  text((xleft + xright) / 2, (ybottom + ytop) / 2, text_value,
       cex = cex, font = font, col = dark, family = "sans", xpd = NA)
}

heatmap_panel <- function(mat, row_labels, col_labels, main = "", zlim = NULL,
                          cex_row = 0.65, cex_col = 0.65, show_values = FALSE,
                          scale_label = "NES") {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (is.null(zlim)) {
    lim <- max(abs(mat), na.rm = TRUE)
    zlim <- c(-lim, lim)
  }
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), col = heat_cols,
        zlim = zlim, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(1, at = seq_len(ncol(mat)), labels = col_labels, las = 2, tick = FALSE,
       cex.axis = cex_col, line = -0.5)
  axis(2, at = seq_len(nrow(mat)), labels = row_labels, las = 2, tick = FALSE,
       cex.axis = cex_row, line = -0.5)
  box(col = grey, lwd = 0.8)
  abline(h = seq(1.5, nrow(mat) - 0.5, by = 1), col = "white", lwd = 0.35)
  abline(v = seq(1.5, ncol(mat) - 0.5, by = 1), col = "white", lwd = 0.35)
  if (show_values) {
    for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
      text(j, i, sprintf("%.2f", mat[i, j]), cex = 0.64,
           col = if (abs(mat[i, j]) > 1.8) "white" else dark)
    }
  }
  title(main = main, adj = 0, font.main = 2, cex.main = 1)
  mtext("Blue: negative; red: positive", side = 3, adj = 1,
        cex = 0.60, col = grey, line = 0.1)
}

pretty_path <- function(x) {
  y <- gsub("HALLMARK_", "", x, fixed = TRUE)
  y <- gsub("_", " ", y, fixed = TRUE)
  tools::toTitleCase(tolower(y))
}

short_contrast <- c(
  "GSE76068_tumor_Escape_vs_Response" = "76068 Tumor",
  "GSE76068_stroma_Escape_vs_Response" = "76068 Stroma",
  "GSE73571_tumor_Resistant_vs_Sensitive" = "73571 Tumor",
  "GSE180687_endothelial_Resistant_vs_Sensitive" = "180687 Endo.",
  "GSE64472_stroma_Cediranib_resistant_vs_sensitive" = "64472 Cedi.",
  "GSE64472_stroma_Vandetanib_resistant_vs_sensitive" = "64472 Vand.",
  "GSE26644_tumor_Bevacizumab_resistant_vs_vehicle" = "26644 Tumor",
  "GSE26644_stroma_Bevacizumab_resistant_vs_vehicle" = "26644 Stroma",
  "GSE81465_G9_vs_G1" = "81465 G9 vs G1",
  "GSE64052_786O_sorafenib_resistant_vs_untreated" = "64052 786-O Sor.",
  "GSE64052_786O_sunitinib_resistant_vs_untreated" = "64052 786-O Sun.",
  "GSE64052_A498_sorafenib_resistant_vs_untreated" = "64052 A498 Sor.",
  "GSE86525_HT29_bevacizumab_resistant_vs_control" = "86525 HT-29 Bev."
)

# Figure 1 -------------------------------------------------------------------
draw_figure1 <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  panel_label("A")
  title("Central question", adj = 0.08, font.main = 2, cex.main = 1.05)
  draw_box(0.08, 0.57, 0.92, 0.82,
           "Does VEGF-pathway resistance\nconverge on one\ntranscriptional endpoint?",
           fill = "#F3F4F6", border = dark, cex = 0.88, font = 2)
  arrows(0.50, 0.56, 0.50, 0.43, length = 0.08, lwd = 1.4, col = grey)
  draw_box(0.08, 0.14, 0.44, 0.39, "Single shared\nendpoint", fill = "#FDECEC", border = vermillion, cex = 0.95)
  draw_box(0.56, 0.14, 0.92, 0.39, "Multiple structured\narchitectures", fill = "#E8F4FA", border = blue, cex = 0.72, font = 2)
  text(0.74, 0.08, "Supported model", col = blue, cex = 0.8, font = 2)

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  panel_label("B")
  title("Evidence hierarchy", adj = 0.08, font.main = 2, cex.main = 1.05)
  ys <- c(0.83, 0.66, 0.49, 0.29, 0.09)
  hs <- c(0.07, 0.075, 0.085, 0.10, 0.075)
  hierarchy <- read_tab(file.path(source_root, "Figure_1", "Figure1B_EvidenceHierarchy.tsv"))
  discovery_scope <- hierarchy$scope[hierarchy$evidence_layer == "Discovery"][1]
  validation_scope <- hierarchy$scope[hierarchy$evidence_layer == "Primary experimental validation"][1]
  labs <- c(paste0("DISCOVERY\n", sub(" resistance", "", discovery_scope)),
            "CROSS-CONTEXT\nPATHWAY ARCHITECTURE",
            "FROZEN RECURRENT /\nDIVERGENT PROGRAMS",
            paste0("PRIMARY VALIDATION\n", validation_scope, "\n+ supportive serial model"),
            "HUMAN RESPONSE HISTORY\nGSE79671 · 16 pairs")
  fills <- c("#DCEAF4", "#E8F4FA", "#DFF2EA", "#E8F4FA", "#FCE8D5")
  for (i in seq_along(ys)) {
    draw_box(0.18, ys[i] - hs[i], 0.82, ys[i] + hs[i], labs[i], fill = fills[i], cex = 0.66, font = 2)
    if (i < length(ys)) arrows(0.5, ys[i] - hs[i] - 0.008, 0.5, ys[i + 1] + hs[i + 1] + 0.008, length = 0.04, col = grey)
  }

  tax <- read_tab(file.path(source_root, "Figure_1", "Figure1C_DatasetTaxonomy.tsv"))
  tax <- tax[tax$analysis_role != "legacy", ]
  par(mar = c(2, 1, 2.5, 0.3))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  panel_label("C")
  title("Dataset taxonomy and independence", adj = 0.08, font.main = 2, cex.main = 1.05)
  headers <- c("Dataset", "Evidence role", "Cancer/model\ncontext", "Compartment")
  xpos <- c(0.03, 0.22, 0.51, 0.79)
  text(xpos, 0.92, headers, adj = 0, font = 2, cex = 0.68, col = dark)
  abline(h = 0.82, col = dark, lwd = 1)
  yy <- seq(0.77, 0.12, length.out = nrow(tax))
  for (i in seq_len(nrow(tax))) {
    raw_role <- tax$analysis_role[i]
    role <- if (grepl("Discovery", raw_role, ignore.case = TRUE)) "Discovery" else if (raw_role == "validation") "Validation" else if (raw_role == "validation-supportive") "Val. support" else if (raw_role == "human-longitudinal") "Human long." else "Baseline support"
    context_map <- c("GSE76068" = "ccRCC PDX", "GSE73571" = "HCC xenograft", "GSE180687" = "Ovarian xenograft", "GSE64472" = "NSCLC xenograft", "GSE26644" = "NSCLC xenograft", "GSE81465" = "GBM xenograft", "GSE64052" = "RCC xenograft", "GSE86525" = "CRC xenograft", "GSE45161" = "GBM xenograft", "GSE79671" = "Recurrent GBM", "GSE37138" = "NSCLC biopsy")
    cancer <- unname(context_map[tax$dataset[i]])
    compartment_map <- c("GSE76068" = "tumor / stroma", "GSE73571" = "tumor", "GSE180687" = "endothelial", "GSE64472" = "mouse stroma", "GSE26644" = "tumor / stroma", "GSE81465" = "tumor", "GSE64052" = "tumor", "GSE86525" = "tumor", "GSE45161" = "tumor", "GSE79671" = "bulk tumor", "GSE37138" = "baseline tumor")
    comp <- unname(compartment_map[tax$dataset[i]])
    text(xpos, yy[i], c(tax$dataset[i], role, cancer, comp), adj = 0, cex = 0.64, col = dark)
  }
  text(0.04, 0.04, "Contrasts sharing a GEO dataset/system are displayed as nested evidence.", adj = 0, cex = 0.60, col = grey)

  par(mar = c(1, 1, 2.5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  panel_label("D")
  title("Analysis strategy", adj = 0.08, font.main = 2, cex.main = 1.05)
  xs <- seq(0.1, 0.9, length.out = 5)
  steps <- c("Within-study\nDE", "Ranked\nGSEA", "Pathway\narchitecture", "Frozen\nprogram\nvalidation", "Human\ntranslation")
  for (i in seq_along(xs)) {
    draw_box(xs[i] - 0.075, 0.48, xs[i] + 0.075, 0.70, steps[i], fill = "white", border = blue, cex = 0.50, font = 2)
    if (i < length(xs)) arrows(xs[i] + 0.08, 0.59, xs[i + 1] - 0.08, 0.59, length = 0.045, col = grey)
  }
  text(0.5, 0.30, "Hallmark primary · Reactome supporting · no pathway reselection",
       cex = 0.72, col = dark)
  text(0.5, 0.18, "Discovery ≠ validation ≠ human longitudinal\n≠ baseline association",
       cex = 0.60, col = vermillion, font = 2)
  mtext("Evidence hierarchy; contrast counts are not independent-study counts", outer = TRUE,
        side = 1, line = -1, cex = 0.68, col = grey)
}
if (should_render("Figure1")) save_dual("Figure1", 10, 7.8, draw_figure1)

# Figure 2 -------------------------------------------------------------------
draw_figure2 <- function() {
  layout(matrix(c(1, 2, 1, 3), nrow = 2, byrow = TRUE), widths = c(2.08, 1.17), heights = c(1, 0.72))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)

  dat <- read_tab(file.path(source_root, "Figure_2", "Figure2A_HallmarkRawNES.tsv"))
  paths <- unique(dat$pathway)
  contrasts <- unique(dat$contrast_id)
  mat <- matrix(as.numeric(dat$NES), nrow = length(paths), ncol = length(contrasts), byrow = FALSE,
                dimnames = list(paths, contrasts))
  # Rebuild explicitly to avoid relying on melt row order.
  for (i in seq_along(paths)) for (j in seq_along(contrasts)) {
    mat[i, j] <- as.numeric(dat$NES[dat$pathway == paths[i] & dat$contrast_id == contrasts[j]][1])
  }
  par(mar = c(6.5, 9, 3.0, 0.5))
  heatmap_panel(mat, sapply(paths, pretty_path), unname(short_contrast[contrasts]),
                main = "A  Hallmark NES: discovery contexts", zlim = c(-3.6, 3.6), cex_row = 0.60, cex_col = 0.66)

  corr <- read_tab(file.path(source_root, "Figure_2", "Figure2B_HallmarkCorrelations.tsv"))
  cm <- matrix(NA_real_, nrow = length(contrasts), ncol = length(contrasts), dimnames = list(contrasts, contrasts))
  for (i in seq_len(nrow(corr))) cm[corr$contrast_id[i], corr$contrast_id_2[i]] <- as.numeric(corr$spearman_rho[i])
  par(mar = c(4.5, 4.5, 3.0, 0.5))
  heatmap_panel(cm, unname(short_contrast[contrasts]), unname(short_contrast[contrasts]),
                main = "B  Correlation matrix", zlim = c(-1, 1), cex_row = 0.52, cex_col = 0.64, show_values = TRUE,
                scale_label = "correlation")

  sm <- read_tab(file.path(source_root, "Figure_2", "Figure2D_CorrelationSummary.tsv"))
  getv <- function(id) as.numeric(sm$exact_value[sm$number_id == id])
  par(mar = c(2, 2, 3, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  panel_label("C")
  title("Pairwise correlation summary", adj = 0.08, font.main = 2, cex.main = 1)
  draw_box(0.08, 0.55, 0.45, 0.83, sprintf("Mean rho\n%.3f", getv("ROUND2B_HALLMARK_PAIRWISE_MEAN_RHO")), fill = "#F3F4F6", cex = 0.82, font = 2)
  draw_box(0.55, 0.55, 0.92, 0.83, sprintf("Median rho\n%.3f", getv("ROUND2B_HALLMARK_PAIRWISE_MEDIAN_RHO")), fill = "#F3F4F6", cex = 0.82, font = 2)
  segments(0.15, 0.31, 0.85, 0.31, lwd = 3, col = lightgrey)
  minv <- getv("ROUND2B_HALLMARK_PAIRWISE_MIN_RHO"); maxv <- getv("ROUND2B_HALLMARK_PAIRWISE_MAX_RHO")
  points(c(0.15, 0.85), c(0.31, 0.31), pch = 19, col = c(blue, vermillion), cex = 1.2)
  text(c(0.15, 0.85), c(0.20, 0.20), sprintf("%.3f", c(minv, maxv)), cex = 0.75)
  text(0.5, 0.08, "Observed range", cex = 0.68, col = grey)
}
if (should_render("Figure2")) save_dual("Figure2", 10, 12, draw_figure2)

# Figure 3 -------------------------------------------------------------------
draw_figure3 <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1.45, 1), heights = c(0.9, 1.25))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)
  rec <- read_tab(file.path(source_root, "Figure_3", "Figure3A_RecurrentProgramsRawNES.tsv"))
  paths <- c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
  contrasts <- unique(rec$contrast_id)
  mat <- matrix(NA_real_, 2, length(contrasts), dimnames = list(paths, contrasts))
  for (i in seq_len(nrow(rec))) mat[rec$pathway[i], rec$contrast_id[i]] <- as.numeric(rec$NES[i])
  par(mar = c(8, 10, 3, 1))
  heatmap_panel(mat, c("MTORC1 signaling", "Unfolded protein response"), unname(short_contrast[contrasts]),
                main = "A  Recurrent programs", zlim = c(-3, 3), cex_row = 0.75, cex_col = 0.66, show_values = TRUE)

  sm <- read_tab(file.path(source_root, "Figure_3", "Figure3B_RecurrentProgramSummary.tsv"))
  sm <- sm[match(paths, sm$pathway), ]
  par(mar = c(5.2, 5.5, 3, 1))
  mids <- barplot(as.numeric(sm$median_NES), names.arg = c("MTORC1", "UPR"), col = c(blue, green), border = NA,
                  ylim = c(0, 1.8), ylab = "Median NES", main = "B  Recurrence summary")
  abline(h = 0, col = grey)
  text(mids, as.numeric(sm$median_NES) + 0.12,
       labels = paste0(sm$n_positive, "/", sm$n_available, " positive\n", sm$TierA_positive, "/", as.integer(sm$TierA_positive) + as.integer(sm$TierA_negative), " Tier A"),
       cex = 0.72, font = 2)
  mtext("Contrasts nested by dataset/system", side = 1, line = 3.4, cex = 0.64, col = grey)

  div <- read_tab(file.path(source_root, "Figure_3", "Figure3C_DivergentProgramsRawNES.tsv"))
  dpaths <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INFLAMMATORY_RESPONSE")
  dm <- matrix(NA_real_, length(dpaths), length(contrasts), dimnames = list(dpaths, contrasts))
  for (i in seq_len(nrow(div))) dm[div$pathway[i], div$contrast_id[i]] <- as.numeric(div$NES[i])
  par(mar = c(8, 10, 3, 1))
  heatmap_panel(dm, c("Interferon-gamma response", "Interferon-alpha response", "Epithelial-mesenchymal transition", "Inflammatory response"),
                unname(short_contrast[contrasts]), main = "C  Context-dependent branches", zlim = c(-3.6, 3.6), cex_row = 0.68, cex_col = 0.66, show_values = TRUE)

  ds <- read_tab(file.path(source_root, "Figure_3", "Figure3D_DivergenceSummary.tsv"))
  ds <- ds[match(dpaths, ds$pathway), ]
  par(mar = c(5.6, 10, 3, 1))
  y <- seq_len(nrow(ds))
  plot(c(-4.1, 4.1), c(0.5, 4.5), type = "n", yaxt = "n", xlab = "Signed NES", ylab = "",
       main = "D  NES ranges")
  axis(2, at = y, labels = c("IFN-gamma", "IFN-alpha", "EMT", "Inflammation"), las = 2, tick = FALSE)
  abline(v = 0, col = grey, lty = 2)
  segments(as.numeric(ds$min_negative_NES), y, as.numeric(ds$max_positive_NES), y, lwd = 3, col = lightgrey)
  points(as.numeric(ds$min_negative_NES), y, pch = 19, col = blue, cex = 1.2)
  points(as.numeric(ds$max_positive_NES), y, pch = 19, col = vermillion, cex = 1.2)
  text(3.95, y + 0.13, labels = paste0(ds$n_positive, "+ / ", ds$n_negative, "−"), adj = 1, cex = 0.66, font = 2)
  mtext("Raw signed NES; negative contexts retained", outer = TRUE, side = 1, line = -1, cex = 0.68, col = grey)
}
if (should_render("Figure3")) save_dual("Figure3", 10, 10.5, draw_figure3)

# Figure 4 -------------------------------------------------------------------
draw_figure4 <- function() {
  layout(matrix(c(1, 1, 2, 1, 1, 3, 4, 4, 4), nrow = 3, byrow = TRUE), widths = c(1.2, 1.2, 1), heights = c(1, 1, 1.05))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)
  val <- read_tab(file.path(source_root, "Figure_4", "Figure4A_PredefinedValidationNES.tsv"))
  if (!"final_evidence_status" %in% names(val)) stop("Figure 4 source data lack final_evidence_status")
  val <- val[val$final_evidence_status == "PRIMARY_VALIDATION", ]
  stopifnot(length(unique(val$dataset)) == 2L, length(unique(val$contrast)) == 4L)
  paths <- c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INFLAMMATORY_RESPONSE")
  contrasts <- unique(val$contrast)
  mat <- matrix(NA_real_, length(paths), length(contrasts), dimnames = list(paths, contrasts))
  for (i in seq_len(nrow(val))) mat[val$pathway[i], val$contrast[i]] <- as.numeric(val$NES[i])
  par(mar = c(9, 11, 3, 1))
  heatmap_panel(mat, c("MTORC1 signaling", "Unfolded protein response", "Interferon-gamma response", "Interferon-alpha response", "Epithelial-mesenchymal transition", "Inflammatory response"),
                unname(short_contrast[contrasts]), main = "A  Primary non-discovery validation", zlim = c(-3.8, 3.8), cex_row = 0.62, cex_col = 0.68, show_values = TRUE)

  plot_validation <- function(pathway, title_text, panel, col_point) {
    d <- val[val$pathway == pathway, ]
    d <- d[match(contrasts, d$contrast), ]
    par(mar = c(5.2, 9, 3, 1))
    y <- seq_len(nrow(d))
    lim <- max(4, max(abs(as.numeric(d$NES))) + 0.5)
    plot(c(-lim - 0.4, lim + 0.8), c(0.5, length(y) + 0.5), type = "n", yaxt = "n", xlab = "NES", ylab = "", main = title_text)
    axis(2, at = y, labels = unname(short_contrast[d$contrast]), las = 2, tick = FALSE, cex.axis = 0.72)
    abline(v = 0, col = grey, lty = 2)
    segments(0, y, as.numeric(d$NES), y, col = lightgrey, lwd = 3)
    point_cols <- ifelse(d$dataset == "GSE86525", vermillion, col_point)
    points(as.numeric(d$NES), y, pch = 21, bg = point_cols, col = "white", cex = 1.4)
    text(as.numeric(d$NES), y, labels = sprintf("%.2f", as.numeric(d$NES)), pos = ifelse(as.numeric(d$NES) >= 0, 4, 2), cex = 0.68)
    mtext("GSE86525 negative", side = 3, line = 0.1, adj = 0.5, cex = 0.48, col = vermillion)
    panel_label(panel)
  }
  plot_validation("HALLMARK_MTORC1_SIGNALING", "B  MTORC1", "", blue)
  plot_validation("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "C  UPR", "", green)

  temp <- read_tab(file.path(source_root, "Figure_4", "Figure4D_GSE81465Temporal.tsv"))
  stages <- c("IgG", "G1", "G4", "G9")
  par(mar = c(5.5, 5, 3, 1))
  plot(c(1, 4), range(as.numeric(temp$program_score_mean), na.rm = TRUE) * 1.35,
       type = "n", xaxt = "n", xlab = "", ylab = "Array-level mean rank score", main = "D  GSE81465 serial xenograft states")
  axis(1, at = 1:4, labels = c("IgG", "G1 early", "G4 intermediate", "G9 late"))
  for (j in seq_along(c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"))) {
    path <- c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")[j]
    d <- temp[temp$pathway == path, ]; d <- d[match(stages, d$stage), ]
    x <- 1:4 + c(-0.03, 0.03)[j]
    y <- as.numeric(d$program_score_mean)
    lines(x, y, col = c(blue, green)[j], lwd = 2)
    points(x, y, pch = 21, bg = c(blue, green)[j], col = "white", cex = 1.25)
  }
  abline(h = 0, col = lightgrey)
  legend("topleft", legend = c("MTORC1 signaling", "Unfolded protein response"), col = c(blue, green), lwd = 2, pch = 19, bty = "n", cex = 0.75)
  text(4, par("usr")[4], "G9 vs G1 direction: positive / positive", adj = c(1, 1.2), cex = 0.68, col = dark)
  mtext("Supportive only: replicate independence unresolved; primary validation = 4 contrasts / 2 datasets (3 GSE64052 contexts nested)",
        side = 1, line = 4.2, cex = 0.64, col = grey)
}
if (should_render("Figure4")) save_dual("Figure4", 10, 11, draw_figure4)

# Figure 5 -------------------------------------------------------------------
trajectory_panel <- function(program, panel, title_text, ylim = NULL) {
  d <- read_tab(file.path(source_root, "Figure_5", paste0("Figure5", panel, "_", program, "_PatientTrajectories.tsv")))
  d$score <- as.numeric(d$score)
  if (is.null(ylim)) ylim <- range(d$score, na.rm = TRUE) + c(-0.02, 0.02)
  plot(c(0.7, 5.3), ylim, type = "n", xaxt = "n", xlab = "", ylab = "Rank-based program score", main = paste0(panel, "  ", title_text))
  axis(1, at = c(1, 2, 4, 5), labels = c("Pre", "Progression", "Pre", "Progression"), cex.axis = 0.58, gap.axis = -1)
  mtext("Responders", side = 1, at = 1.5, line = 2.2, cex = 0.64, col = blue)
  mtext("No response", side = 1, at = 4.5, line = 2.2, cex = 0.50, col = orange)
  abline(v = 3, col = lightgrey)
  for (group in c("Responder", "NonResponder")) {
    subs <- unique(d$subject_id[d$responder_status == group])
    for (subject in subs) {
      q <- d[d$subject_id == subject, ]; q <- q[match(c("Pre", "Post"), q$before_after), ]
      x <- if (group == "Responder") c(1, 2) else c(4, 5)
      col <- if (group == "Responder") adjustcolor(blue, 0.60) else adjustcolor(orange, 0.60)
      lines(x, q$score, col = col, lwd = 1.2)
      points(x, q$score, pch = 21, bg = col, col = "white", cex = 0.85)
    }
  }
}

draw_figure5 <- function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(1, 1, 1), heights = c(1, 1))
  par(family = "sans", fg = dark, col.axis = dark, col.lab = dark)
  par(mar = c(2, 1, 3, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("A  Human longitudinal design", adj = 0, font.main = 2, cex.main = 0.86)
  draw_box(0.03, 0.72, 0.97, 0.88, "GSE79671 · 16 paired patients", fill = "#F3F4F6", font = 2, cex = 0.78)
  arrows(0.50, 0.71, 0.50, 0.61, length = 0.06, col = grey)
  draw_box(0.08, 0.32, 0.46, 0.62, "Initial\nresponders\nn = 6\nPre to\nprogression", fill = "#DCEAF4", border = blue, cex = 0.66, font = 2)
  draw_box(0.54, 0.32, 0.92, 0.62, "No objective\nresponse\nn = 10\nPre to\nprogression", fill = "#FCE8D5", border = orange, cex = 0.66, font = 2)
  text(0.5, 0.22, "Paired within patient", cex = 0.72, font = 2)
  text(0.5, 0.12, "Response group × time interaction", cex = 0.68, col = grey)
  text(0.5, 0.04, "49 vs 0 FDR-significant DEGs", cex = 0.68, col = dark)

  par(mar = c(5, 4.3, 3, 1)); trajectory_panel("MTORC1", "B", "MTORC1 signaling")
  par(mar = c(5, 4.3, 3, 1), cex.main = 0.86); trajectory_panel("UPR", "C", "Unfolded protein response")
  par(mar = c(5, 4.3, 3, 1)); trajectory_panel("EMT", "D", "EMT")
  par(mar = c(5, 4.3, 3, 1)); trajectory_panel("IFNG", "E", "IFN-gamma response")

  d <- read_tab(file.path(source_root, "Figure_5", "Figure5F_PatientDeltas.tsv"))
  stats <- read_tab(file.path(source_root, "Figure_5", "Figure5F_CountBasedInteraction.tsv"))
  programs <- c("MTORC1", "UPR", "IFNG", "IFNA", "EMT", "INFLAMMATORY_RESPONSE")
  labels <- c("MTORC1", "UPR", "IFN-gamma", "IFN-alpha", "EMT", "Inflammation")
  par(mar = c(6.5, 4.3, 3, 1))
  rng <- range(as.numeric(d$delta[d$program %in% programs]), na.rm = TRUE) * 1.15
  plot(c(0.5, 6.5), rng, type = "n", xaxt = "n", xlab = "", ylab = "Patient-level pathway-score change", main = "F  Patient score changes")
  axis(1, at = 1:6, labels = labels, las = 2, cex.axis = 0.72)
  abline(h = 0, col = grey, lty = 2)
  set.seed(1)
  for (i in seq_along(programs)) {
    for (group in c("Responder", "NonResponder")) {
      q <- as.numeric(d$delta[d$program == programs[i] & d$responder_status == group])
      offset <- if (group == "Responder") -0.15 else 0.15
      col <- if (group == "Responder") blue else orange
      x <- rep(i + offset, length(q)) + seq(-0.045, 0.045, length.out = length(q))
      points(x, q, pch = 21, bg = adjustcolor(col, 0.65), col = "white", cex = 0.76)
      segments(i + offset - 0.09, median(q), i + offset + 0.09, median(q), col = col, lwd = 2.5)
    }
    h3 <- as.numeric(stats$H3_Responder_change_vs_NonResponder_change_NES[stats$program == programs[i]])
    text(i, rng[2], sprintf("%.2f", h3), adj = c(0.5, 1.15), cex = 0.70, col = dark)
  }
  legend("bottom", legend = c("Responder", "No response"), pt.bg = c(blue, orange), pch = 21, bty = "n", cex = 0.70, horiz = TRUE)
  mtext("Interaction NES", side = 3, line = 0.1, adj = 0, cex = 0.64, col = grey)
  mtext("Patient-level scores; annotations show group-level interaction GSEA NES", outer = TRUE, side = 1, line = -1, cex = 0.60, col = grey)
}
if (should_render("Figure5")) save_dual("Figure5", 10, 11, draw_figure5)

# Figure 6 -------------------------------------------------------------------
draw_figure6 <- function() {
  par(mar = c(1, 1, 1, 1), family = "sans", fg = dark)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.96, "Integrated model of VEGF-pathway resistance architecture", font = 2, cex = 1.24, col = dark)
  text(0.5, 0.92, "Conceptual synthesis · associations, not causal mechanisms", cex = 0.74, col = grey)
  draw_box(0.30, 0.78, 0.70, 0.88, "VEGF-pathway inhibition\nand adaptive pressure", fill = "#F3F4F6", border = dark, cex = 0.92, font = 2)
  arrows(0.5, 0.77, 0.5, 0.72, length = 0.05, lty = 2, col = grey, lwd = 1.4)
  text(0.54, 0.745, "associated states", cex = 0.66, col = grey, adj = 0)
  draw_box(0.055, 0.43, 0.465, 0.70,
           "RECURRENT EXPERIMENTAL\nPROGRAMS\n\nMTORC1 signaling\nUnfolded protein response\n\nrecurrent · not universal",
           fill = "#DCEAF4", border = blue, cex = 0.68, font = 2)
  draw_box(0.535, 0.43, 0.945, 0.70,
           "CONTEXT-DEPENDENT\nBRANCHES\n\nInterferon responses\nEpithelial-mesenchymal transition\nInflammatory response\n\nbidirectional / magnitude-dependent",
           fill = "#DFF2EA", border = green, cex = 0.66, font = 2)
  segments(0.5, 0.72, 0.26, 0.70, lty = 2, col = grey, lwd = 1.2)
  segments(0.5, 0.72, 0.74, 0.70, lty = 2, col = grey, lwd = 1.2)
  rect(0.08, 0.15, 0.92, 0.38, col = "#FFF7ED", border = orange, lwd = 1.7)
  text(0.50, 0.355, "HUMAN RESPONSE HISTORY", col = orange, font = 2, cex = 0.86)
  draw_box(0.14, 0.20, 0.46, 0.30, "Initial objective response\nto progression", fill = "white", border = blue, cex = 0.72)
  draw_box(0.54, 0.20, 0.86, 0.30, "No objective response\nto progression", fill = "white", border = orange, cex = 0.70)
  text(0.50, 0.26, "≠", cex = 1.7, font = 2, col = vermillion)
  segments(0.26, 0.43, 0.31, 0.38, lty = 2, col = grey)
  segments(0.74, 0.43, 0.69, 0.38, lty = 2, col = grey)
  text(0.50, 0.10, "The same clinical endpoint of treatment failure can arise\nthrough different transcriptional states.",
       cex = 0.72, font = 2, col = dark)
  text(0.50, 0.04, "MTORC1/UPR directions in responder progression oppose their dominant experimental direction.",
       cex = 0.72, col = vermillion)
  text(0.015, 0.91, "A", adj = c(0, 1), font = 2, cex = 1.0)
  mtext("Conceptual synthesis; dashed connectors denote association, not causation", side = 1, line = -1, cex = 0.68, col = grey)
}
if (should_render("Figure6")) save_dual("Figure6", 10, 8.0, draw_figure6)

# Supporting GSE37138 figure --------------------------------------------------
draw_gse37138 <- function() {
  d <- read_tab(file.path(source_root, "Supplementary", "GSE37138_MTORC1_UPR_Scores.tsv"))
  assoc <- read_tab(file.path(source_root, "Supplementary", "GSE37138_MTORC1_UPR_Associations.tsv"))
  layout(matrix(1:2, nrow = 1))
  for (program in c("MTORC1", "UPR")) {
    q <- d[d$program == program & !is.na(d$tumor_shrinkage_week12_pct), ]
    x <- as.numeric(q$program_score); y <- as.numeric(q$tumor_shrinkage_week12_pct)
    row <- assoc[assoc$program == program, ]
    par(mar = c(4.5, 4.5, 3, 1))
    plot(x, y, pch = 21, bg = if (program == "MTORC1") blue else green, col = "white", cex = 1.0,
         xlab = paste(program, "baseline program score"), ylab = "Week-12 tumor shrinkage (%)",
         main = paste(program, "supporting baseline association"))
    abline(lm(y ~ x), col = grey, lty = 2)
    legend("topright", legend = sprintf("Spearman rho = %.2f\nP = %.3f\nFDR = %.3f\nn = %s",
                                        as.numeric(row$spearman_rho), as.numeric(row$spearman_p), as.numeric(row$continuous_FDR), row$continuous_n),
           bty = "n", cex = 0.72)
  }
  mtext("Single-arm bevacizumab + erlotinib; supporting association, not predictive validation",
        outer = TRUE, side = 1, line = -1, cex = 0.68, col = grey)
}
if (should_render("FigureS4")) save_dual("FigureS4_GSE37138_baseline_associations", 11, 5.5, draw_gse37138, supplementary = TRUE)

# Supplementary Figure S1: descriptive PCA moved out of the main figure.
draw_s1 <- function() {
  pca <- read_tab(file.path(source_root, "Figure_2", "Figure2C_HallmarkPCA.tsv"))
  cols <- c("tumor" = blue, "stroma" = orange, "endothelial" = green)
  par(mar = c(5, 5, 3, 1), family = "sans")
  plot(as.numeric(pca$PC1), as.numeric(pca$PC2), pch = 21, bg = cols[pca$compartment], col = "white", lwd = 0.8,
       cex = 1.8, xlab = "PC1 score", ylab = "PC2 score", main = "Discovery Hallmark pathway-space PCA")
  abline(h = 0, v = 0, col = lightgrey, lty = 2)
  text(as.numeric(pca$PC1), as.numeric(pca$PC2), labels = pca$contrast_display, pos = 3, cex = 0.68, col = dark)
  legend("bottomright", legend = names(cols), pt.bg = cols, pch = 21, bty = "n", cex = 0.8)
  mtext("Descriptive ordination of the same raw 50-pathway NES matrix shown in Figure 2", side = 1, line = 4, cex = 0.68, col = grey)
}
if (should_render("FigureS1")) save_dual("FigureS1_Discovery_Hallmark_PCA", 9, 7, draw_s1, supplementary = TRUE)

# Supplementary Figure S2: frozen recurrence robustness fields.
draw_s2 <- function() {
  d <- read_tab(file.path(source_root, "Figure_3", "Figure3B_RecurrentProgramSummary.tsv"))
  d <- d[match(c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"), d$pathway), ]
  mat <- rbind(as.numeric(d$direction_consistency), as.numeric(d$TierA_direction_consistency), as.numeric(d$lodo_min_direction_consistency))
  colnames(mat) <- c("MTORC1", "UPR"); rownames(mat) <- c("All contexts", "Tier A", "Minimum LODO")
  par(mar = c(5, 5, 3, 1), family = "sans")
  bp <- barplot(mat, beside = TRUE, col = c(blue, sky, orange), border = NA, ylim = c(0, 1.08),
                ylab = "Positive-direction consistency", main = "Robustness of recurrent-program direction")
  abline(h = c(0.5, 1), col = lightgrey, lty = c(2, 1))
  text(bp, mat + 0.035, sprintf("%.2f", mat), cex = 0.75)
  legend("bottomleft", legend = rownames(mat), fill = c(blue, sky, orange), bty = "n", cex = 0.8)
  mtext("LODO, leave one dataset out; summaries are descriptive", side = 1, line = 4, cex = 0.68, col = grey)
}
if (should_render("FigureS2")) save_dual("FigureS2_Recurrent_program_robustness", 8, 6, draw_s2, supplementary = TRUE)

# Supplementary Figure S3: full frozen validation FDR display.
draw_s3 <- function() {
  d <- read_tab(file.path(source_root, "Figure_4", "Figure4A_PredefinedValidationNES.tsv"))
  if (!"final_evidence_status" %in% names(d)) stop("Figure S3 source data lack final_evidence_status")
  d <- d[d$final_evidence_status == "PRIMARY_VALIDATION", ]
  stopifnot(length(unique(d$dataset)) == 2L, length(unique(d$contrast)) == 4L)
  paths <- unique(d$pathway); contrasts <- unique(d$contrast)
  mat <- matrix(NA_real_, length(paths), length(contrasts), dimnames = list(paths, contrasts))
  for (i in seq_len(nrow(d))) mat[d$pathway[i], d$contrast[i]] <- -log10(max(as.numeric(d$padj[i]), 1e-300))
  par(mar = c(10, 12, 3, 1), family = "sans")
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), col = colorRampPalette(c("white", "#FDD0A2", vermillion))(101), axes = FALSE,
        xlab = "", ylab = "", useRaster = TRUE)
  axis(1, at = seq_len(ncol(mat)), labels = unname(short_contrast[contrasts]), las = 2, tick = FALSE, cex.axis = 0.75)
  axis(2, at = seq_len(nrow(mat)), labels = sapply(paths, pretty_path), las = 2, tick = FALSE, cex.axis = 0.75)
  fmt_fdr <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) text(j, i, fmt_fdr(as.numeric(d$padj[d$pathway == paths[i] & d$contrast == contrasts[j]][1])), cex = 0.68)
  box(col = grey); title("Primary non-discovery validation FDR values", adj = 0, font.main = 2)
  mtext("Cells show frozen Benjamini–Hochberg FDR; color intensity is −log10(FDR)", side = 3, adj = 1, cex = 0.68, col = grey)
}
if (should_render("FigureS3")) save_dual("FigureS3_Validation_FDR_matrix", 10, 7, draw_s3, supplementary = TRUE)

# Supplementary Figure S5: all frozen human group-level contrasts.
draw_s5 <- function() {
  d <- read_tab(file.path(source_root, "Figure_5", "Figure5F_CountBasedInteraction.tsv"))
  programs <- c("MTORC1", "UPR", "IFNG", "IFNA", "EMT", "INFLAMMATORY_RESPONSE")
  d <- d[match(programs, d$program), ]
  mat <- cbind(as.numeric(d$H1_Responder_Post_vs_Pre_NES), as.numeric(d$H2_NonResponder_Post_vs_Pre_NES),
               as.numeric(d$H3_Responder_change_vs_NonResponder_change_NES), as.numeric(d$H4_Responder_Pre_vs_NonResponder_Pre_NES))
  rownames(mat) <- c("MTORC1", "UPR", "IFN-gamma", "IFN-alpha", "EMT", "Inflammation")
  colnames(mat) <- c("H1: responder\nprogression vs pre", "H2: no objective response\nprogression vs pre", "H3: interaction", "H4: baseline\ngroup contrast")
  par(mar = c(9.5, 9, 4.5, 1), family = "sans")
  heatmap_panel(mat, rownames(mat), colnames(mat), main = "GSE79671 group-level count-based GSEA", zlim = c(-3.5, 3.5), cex_row = 0.8, cex_col = 0.62, show_values = TRUE)
  mtext("These group-level NES values are distinct from patient-level rank-score changes", side = 1, line = 8.2, cex = 0.66, col = grey)
}
if (should_render("FigureS5")) save_dual("FigureS5_GSE79671_group_level_GSEA", 9, 8, draw_s5, supplementary = TRUE)

cat("Rendered six main publication figures and five supplementary figures.\n")
