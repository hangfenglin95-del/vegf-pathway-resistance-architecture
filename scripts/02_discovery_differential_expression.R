#!/usr/bin/env Rscript

# Within-dataset discovery limma models and complete ranked lists.
# This script deliberately contains no cross-dataset comparison or integration.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
root <- normalizePath(analysis_root_env, mustWork = TRUE)
suppressPackageStartupMessages(library(limma))

in_dir <- file.path(root, "data", "processed_round2a")
out_dir <- file.path(root, "results", "round2a", "differential")
qc_dir <- file.path(root, "results", "round2a", "qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_expr <- function(unit, suffix = "gene_level_expression.tsv") {
  x <- read.delim(file.path(in_dir, paste0(unit, "_", suffix)), check.names = FALSE)
  rn <- x[[1]]; m <- as.matrix(x[-1]); storage.mode(m) <- "double"; rownames(m) <- rn; m
}
read_meta <- function(unit) read.delim(file.path(in_dir, paste0(unit, "_metadata.tsv")), check.names = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

specs <- list(
  GSE76068_tumor = list(dataset="GSE76068", compartment="tumor", species="Homo sapiens", drug="sunitinib", tier="Tier A", type="blocked",
    contrasts=c(Response_vs_Pretreatment="stateresponse-statepretreatment", Escape_vs_Response="stateescape-stateresponse", Escape_vs_Pretreatment="stateescape-statepretreatment")),
  GSE76068_stroma = list(dataset="GSE76068", compartment="stroma", species="Mus musculus", drug="sunitinib", tier="Tier A", type="blocked",
    contrasts=c(Response_vs_Pretreatment="stateresponse-statepretreatment", Escape_vs_Response="stateescape-stateresponse", Escape_vs_Pretreatment="stateescape-statepretreatment")),
  GSE73571_tumor = list(dataset="GSE73571", compartment="tumor", species="Homo sapiens", drug="sorafenib", tier="Tier A", type="groups",
    contrasts=c(Resistant_vs_Sensitive="groupresistant-groupsensitive")),
  GSE180687_endothelial = list(dataset="GSE180687", compartment="endothelial", species="Mus musculus", drug="B20 anti-VEGF-A antibody", tier="Tier A", type="groups180",
    contrasts=c(Resistant_vs_Sensitive="groupresistant-groupsensitive")),
  GSE64472_stroma = list(dataset="GSE64472", compartment="stroma", species="Mus musculus", drug="cediranib/vandetanib", tier="Tier A", type="groups644",
    contrasts=c(Cediranib_resistant_vs_sensitive="groupcediranib_resistant-groupcediranib_sensitive", Vandetanib_resistant_vs_sensitive="groupvandetanib_resistant-groupvandetanib_sensitive")),
  GSE26644_tumor = list(dataset="GSE26644", compartment="tumor", species="Homo sapiens", drug="bevacizumab", tier="Tier B", type="groups266",
    contrasts=c(Bevacizumab_resistant_vs_vehicle="groupresistant-groupvehicle")),
  GSE26644_stroma = list(dataset="GSE26644", compartment="stroma", species="Mus musculus", drug="bevacizumab", tier="Tier B", type="groups266",
    contrasts=c(Bevacizumab_resistant_vs_vehicle="groupresistant-groupvehicle"))
)

make_design <- function(meta, type) {
  if (type == "blocked") {
    state <- factor(meta$treatment_state, levels=c("pretreatment","response","escape"))
    block <- factor(meta$biological_replicate, levels=c("BR1","BR2","BR3","BR4"))
    return(model.matrix(~0 + state + block))
  }
  if (type == "groups") {
    group <- factor(ifelse(grepl("acquired resistant", meta$treatment_state), "resistant", "sensitive"), levels=c("sensitive","resistant"))
  } else if (type == "groups180") {
    group <- factor(ifelse(grepl("control", meta$biological_group), "control", ifelse(grepl("sensitive", meta$biological_group), "sensitive", "resistant")), levels=c("control","sensitive","resistant"))
  } else if (type == "groups644") {
    group <- factor(ifelse(grepl("cediranib_sensitive", meta$biological_group), "cediranib_sensitive",
                    ifelse(grepl("cediranib_acquired", meta$biological_group), "cediranib_resistant",
                    ifelse(grepl("vandetanib_sensitive", meta$biological_group), "vandetanib_sensitive", "vandetanib_resistant"))))
  } else if (type == "groups266") {
    group <- factor(ifelse(grepl("vehicle", meta$biological_group), "vehicle", "resistant"), levels=c("vehicle","resistant"))
  }
  model.matrix(~0 + group)
}

design_summaries <- list(); de_summaries <- list(); sensitivity <- list(); all_stats <- list()
for (unit in names(specs)) {
  s <- specs[[unit]]; expr <- read_expr(unit); med <- read_expr(unit, "gene_level_expression_median_probe.tsv"); meta <- read_meta(unit)
  stopifnot(identical(colnames(expr), meta$analysis_sample_id), identical(rownames(expr), rownames(med)))
  design <- make_design(meta, s$type)
  if (qr(design)$rank != ncol(design)) stop("Non-full-rank design: ", unit)
  contrast <- makeContrasts(contrasts=unname(s$contrasts), levels=design)
  colnames(contrast) <- names(s$contrasts)
  write_tsv(data.frame(sample_id=meta$analysis_sample_id, design, check.names=FALSE), file.path(out_dir, paste0(unit, "_design_matrix.tsv")))
  write_tsv(data.frame(coefficient=rownames(contrast), contrast, check.names=FALSE), file.path(out_dir, paste0(unit, "_contrast_matrix.tsv")))
  fit <- eBayes(contrasts.fit(lmFit(expr, design), contrast), trend=TRUE, robust=TRUE)
  fit_med <- eBayes(contrasts.fit(lmFit(med, design), contrast), trend=TRUE, robust=TRUE)
  selected <- read.delim(file.path(in_dir, paste0(unit, "_selected_probe_per_gene.tsv")), check.names=FALSE)
  gid <- setNames(selected$gene_id, selected$gene_symbol)
  for (cn in colnames(contrast)) {
    stem <- paste0(unit, "_", cn)
    tt <- topTable(fit, coef=cn, number=Inf, sort.by="none")
    tt$gene_symbol <- rownames(tt); tt$gene_id <- unname(gid[tt$gene_symbol])
    tt <- tt[, c("gene_symbol","gene_id","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
    tt <- tt[order(tt$P.Value, -abs(tt$logFC)), ]
    tt$dataset <- s$dataset; tt$compartment <- s$compartment; tt$species <- s$species; tt$drug <- s$drug
    tt$contrast <- cn; tt$evidence_tier <- if (grepl("Response_vs_Pretreatment|Escape_vs_Pretreatment", cn)) "Temporal/supportive" else s$tier
    parts <- strsplit(cn, "_vs_", fixed=TRUE)[[1]]
    tt$numerator <- parts[1]; tt$denominator <- parts[2]
    write_tsv(tt, file.path(out_dir, paste0(stem, "_gene_stats.tsv")))
    write_tsv(tt, file.path(out_dir, paste0(stem, "_gene_statistics.tsv")))
    rnk <- tt[is.finite(tt$t), c("gene_symbol","t")]; rnk <- rnk[order(rnk$t, decreasing=TRUE), ]
    write.table(rnk, file.path(out_dir, paste0(stem, ".rnk")), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(rnk, file.path(out_dir, paste0(stem, "_ranked_genes.rnk")), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    all_stats[[stem]] <- tt
    de_summaries[[stem]] <- data.frame(contrast_id=stem, dataset=s$dataset, compartment=s$compartment,
      evidence_tier=unique(tt$evidence_tier), n_samples=nrow(design), n_genes=nrow(tt),
      nominal_p_lt_0_05=sum(tt$P.Value < .05), fdr_lt_0_05=sum(tt$adj.P.Val < .05),
      fdr_lt_0_05_up=sum(tt$adj.P.Val < .05 & tt$logFC > 0), fdr_lt_0_05_down=sum(tt$adj.P.Val < .05 & tt$logFC < 0),
      abs_logFC_ge_1=sum(abs(tt$logFC) >= 1), stringsAsFactors=FALSE)
    tt2 <- topTable(fit_med, coef=cn, number=Inf, sort.by="none")
    common <- intersect(rownames(tt2), tt$gene_symbol)
    primary_t <- setNames(tt$t, tt$gene_symbol)[common]; median_t <- setNames(tt2[common, "t"], common)
    sensitivity[[stem]] <- data.frame(contrast_id=stem, common_genes=length(common),
      moderated_t_spearman=cor(primary_t, median_t, method="spearman", use="complete.obs"),
      top100_sign_concordance=mean(sign(primary_t[names(sort(abs(primary_t), decreasing=TRUE))[1:min(100,length(common))]]) == sign(median_t[names(sort(abs(primary_t), decreasing=TRUE))[1:min(100,length(common))]])),
      stringsAsFactors=FALSE)
  }
  design_summaries[[unit]] <- data.frame(unit=unit, n_samples=nrow(design), n_coefficients=ncol(design), design_rank=qr(design)$rank,
    residual_df=nrow(design)-qr(design)$rank, full_rank=TRUE, stringsAsFactors=FALSE)
}

# Predeclared sensitivity check for the GSE180687 sample flagged by joint RMA QC.
unit <- "GSE180687_endothelial"; expr <- read_expr(unit); meta <- read_meta(unit); keep <- meta$analysis_sample_id != "GSM5468046"
design <- make_design(meta[keep, ], "groups180")
fit <- eBayes(contrasts.fit(lmFit(expr[, keep], design), makeContrasts(groupresistant-groupsensitive, levels=design)), trend=TRUE, robust=TRUE)
loo <- topTable(fit, number=Inf, sort.by="none"); main <- all_stats[["GSE180687_endothelial_Resistant_vs_Sensitive"]]
main_t <- setNames(main$t, main$gene_symbol)[rownames(loo)]
loo_summary <- data.frame(excluded_sample="GSM5468046", reason="joint RMA PCA and correlation flag", retained_sensitive_n=3, retained_resistant_n=4,
  moderated_t_spearman=cor(main_t, loo$t, method="spearman"), top100_main_sign_concordance=mean(sign(main_t[names(sort(abs(main_t), decreasing=TRUE))[1:100]]) == sign(loo[names(sort(abs(main_t), decreasing=TRUE))[1:100],"t"])))
write_tsv(data.frame(gene_symbol=rownames(loo), loo), file.path(qc_dir, "GSE180687_leave_GSM5468046_out_gene_stats.tsv"))
write_tsv(loo_summary, file.path(qc_dir, "GSE180687_leave_GSM5468046_out_summary.tsv"))

write_tsv(do.call(rbind, design_summaries), file.path(root, "reports", "round2a_design_diagnostics.tsv"))
write_tsv(do.call(rbind, de_summaries), file.path(root, "reports", "round2a_differential_summary.tsv"))
write_tsv(do.call(rbind, de_summaries), file.path(out_dir, "round2a_DE_summary.tsv"))
write_tsv(do.call(rbind, sensitivity), file.path(root, "reports", "round2a_probe_collapse_sensitivity.tsv"))
writeLines(capture.output(sessionInfo()), file.path(root, "reports", "round2a_differential_sessionInfo.txt"))
message("Round 2A differential analysis completed: ", length(de_summaries), " contrasts")
