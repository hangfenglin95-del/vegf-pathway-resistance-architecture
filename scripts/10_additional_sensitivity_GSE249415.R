#!/usr/bin/env Rscript

# Additional-dataset sensitivity analysis performed after the discovery and
# validation programs had been frozen.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the project root can be resolved.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
input_env <- Sys.getenv("VEGF_GSE249415_INPUT", "")
if (!nzchar(input_env)) stop("Set VEGF_GSE249415_INPUT to GSE249415_raw_HS.txt.gz.")
input <- normalizePath(input_env, mustWork = TRUE)
bundle_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(msigdbr)
})

output_root_env <- Sys.getenv("VEGF_OUTPUT_ROOT", "")
output_root <- if (nzchar(output_root_env)) normalizePath(output_root_env, mustWork = FALSE) else file.path(bundle_root, "outputs")
out_dir <- file.path(output_root, "additional_sensitivity", "GSE249415")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

raw <- read.delim(gzfile(input), check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "")
sample_cols <- c("control.1", "control.2", "control.3", "B20.1", "B20.2", "B20.3")
expr <- as.matrix(raw[, sample_cols])
storage.mode(expr) <- "double"
expr <- normalizeBetweenArrays(log2(expr + 1), method = "quantile")

symbols <- toupper(trimws(raw$SYMBOL))
symbols[symbols %in% c("", "---", "NA", "NULL")] <- NA_character_
symbols[grepl("///|//|;|,", symbols)] <- NA_character_
keep <- !is.na(symbols) & apply(expr, 1, function(v) all(is.finite(v)))
expr <- expr[keep, , drop = FALSE]
symbols <- symbols[keep]
iq <- apply(expr, 1, IQR)
ord <- order(symbols, -iq)
expr <- expr[ord, , drop = FALSE]
symbols <- symbols[ord]
expr <- expr[!duplicated(symbols), , drop = FALSE]
rownames(expr) <- symbols[!duplicated(symbols)]

group <- factor(c(rep("control", 3), rep("B20", 3)), levels = c("control", "B20"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- eBayes(contrasts.fit(lmFit(expr, design), makeContrasts(B20 - control, levels = design)),
              robust = TRUE)
tt <- topTable(fit, number = Inf, sort.by = "none")
tt$gene_symbol <- rownames(tt)
tt <- tt[, c("gene_symbol", setdiff(colnames(tt), "gene_symbol"))]
write_tsv(tt, file.path(out_dir, "GSE249415_B20_vs_IgG_gene_statistics.tsv"))

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
ranks <- fit$t[, 1]
names(ranks) <- rownames(fit$t)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
set.seed(81465)
fg <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500,
                      eps = 0)
fg <- as.data.frame(fg)
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))

frozen <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
out <- fg[match(frozen, fg$pathway), c("pathway", "NES", "pval", "padj", "size", "leadingEdge")]
out$dataset <- "GSE249415"
out$contrast <- "B20-treated adaptive-resistance endpoint vs IgG control"
out$n_numerator <- 3L
out$n_denominator <- 3L
out$direction <- ifelse(out$NES > 0, "positive", "negative")
out$evidence_role <- "additional-dataset sensitivity; not discovery and not primary validation"
out$identified_date <- "2026-09-11"
out <- out[, c("dataset", "contrast", "pathway", "NES", "pval", "padj", "size",
               "direction", "n_numerator", "n_denominator", "evidence_role",
               "identified_date", "leadingEdge")]
write_tsv(out, file.path(out_dir, "GSE249415_frozen_program_sensitivity.tsv"))
write_tsv(data.frame(sample = sample_cols, group = as.character(group)),
          file.path(out_dir, "GSE249415_sample_manifest.tsv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("GSE249415 additional sensitivity analysis complete\n")
