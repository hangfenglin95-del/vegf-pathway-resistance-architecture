#!/usr/bin/env Rscript

# Rebuild role-dependent validation summaries from frozen DE/GSEA outputs only.
# This script does not preprocess expression, fit DE models, or run GSEA.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
phase_root <- normalizePath(analysis_root_env, mustWork = TRUE)
r3 <- file.path(phase_root, "round3a_validation")
out_root_env <- Sys.getenv("ROUND3A_SUMMARY_OUTPUT_ROOT", "")
out_root <- if (nzchar(out_root_env)) normalizePath(out_root_env, mustWork = TRUE) else r3

write_tsv <- function(x, p) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

role_path <- normalizePath(file.path(dirname(script_path), "..", "config", "validation_contrast_roles.tsv"), mustWork = TRUE)
if (!file.exists(role_path)) stop("Missing role config: ", role_path)
role <- read.delim(role_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("dataset", "contrast", "final_evidence_status", "primary_validation_eligible", "evidence_tier")
if (!all(required %in% names(role))) stop("Role config missing required fields")
if (anyDuplicated(role$contrast)) stop("Role config contains duplicate contrasts")
primary_role <- role[role$final_evidence_status == "PRIMARY_VALIDATION", , drop = FALSE]
if (length(unique(primary_role$dataset)) != 2L) stop("primary_dataset_count != 2")
if (nrow(primary_role) != 4L) stop("primary_contrast_count != 4")
if (any(primary_role$dataset == "GSE81465")) stop("GSE81465 entered PRIMARY_VALIDATION")
g814_role <- role[role$dataset == "GSE81465", , drop = FALSE]
if (nrow(g814_role) != 1L || g814_role$final_evidence_status != "SUPPORTIVE_ONLY" ||
    g814_role$primary_validation_eligible != "NO" || g814_role$evidence_tier != "Supportive-only") {
  stop("GSE81465 supportive-only assertion failed")
}

frozen <- c(
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_E2F_TARGETS")
frozen_roles <- setNames(c("PRIMARY_SHARED_CANDIDATE", "SECONDARY_SHARED_CANDIDATE",
  rep("PRIMARY_DIVERGENT_CANDIDATE", 4), "CAUTIONARY_DIVERGENT_CANDIDATE"), frozen)

meta <- data.frame(
  contrast = c(
    "GSE81465_G9_vs_G1",
    "GSE64052_786O_sorafenib_resistant_vs_untreated",
    "GSE64052_786O_sunitinib_resistant_vs_untreated",
    "GSE64052_A498_sorafenib_resistant_vs_untreated",
    "GSE86525_HT29_bevacizumab_resistant_vs_control",
    "GSE45161_NSC11R_vs_parental_control", "GSE45161_NSC11R_vs_parental_bevacizumab",
    "GSE45161_U87R_vs_parental_control", "GSE45161_U87R_vs_parental_bevacizumab"),
  dataset = c("GSE81465", rep("GSE64052", 3), "GSE86525", rep("GSE45161", 4)),
  drug = c("bevacizumab", "sorafenib", "sunitinib", "sorafenib", "bevacizumab",
           "none at harvest", "bevacizumab", "none at harvest", "bevacizumab"),
  cancer_type = c("glioblastoma", rep("renal cell carcinoma", 3), "colorectal cancer", rep("glioblastoma", 4)),
  same_experimental_system_id = c(
    "GSE81465_U87_bevacizumab_generational",
    "GSE64052_786O_sorafenib", "GSE64052_786O_sunitinib", "GSE64052_A498_sorafenib",
    "GSE86525_HT29_bevacizumab", rep("GSE45161_NSC11_derivative", 2), rep("GSE45161_U87_derivative", 2)),
  stringsAsFactors = FALSE)

pre_list <- list()
for (id in meta$contrast) {
  gsea_path <- file.path(r3, "results", "gsea", paste0(id, "_hallmark_full.tsv"))
  if (!file.exists(gsea_path)) stop("Missing frozen Hallmark result: ", gsea_path)
  g <- read.delim(gsea_path, check.names = FALSE, stringsAsFactors = FALSE)
  z <- g[match(frozen, g$pathway), c("pathway", "NES", "pval", "padj"), drop = FALSE]
  if (anyNA(z$pathway)) stop("Frozen pathway missing for ", id)
  m <- meta[meta$contrast == id, , drop = FALSE]
  rr <- role[role$contrast == id, , drop = FALSE]
  final_status <- if (nrow(rr)) rr$final_evidence_status else "SUPPORTIVE_ONLY"
  evidence_tier <- if (nrow(rr)) rr$evidence_tier else "Supportive-only"
  z$dataset <- m$dataset
  z$contrast <- id
  z$drug <- m$drug
  z$cancer_type <- m$cancer_type
  z$evidence_tier <- evidence_tier
  z$same_experimental_system_id <- m$same_experimental_system_id
  z$direction <- ifelse(z$NES > 0, "positive", "negative")
  z$predefined_role <- unname(frozen_roles[z$pathway])
  z$validation_interpretation <- ifelse(final_status == "PRIMARY_VALIDATION",
    "eligible primary validation contrast",
    ifelse(m$dataset == "GSE81465",
      "supportive descriptive evidence; replicate provenance unresolved; not biological-replication validation",
      "supportive directional evidence; excluded from primary denominator"))
  z$final_evidence_status <- final_status
  pre_list[[id]] <- z[, c("pathway", "dataset", "contrast", "drug", "cancer_type", "evidence_tier",
    "same_experimental_system_id", "NES", "pval", "padj", "direction", "predefined_role",
    "validation_interpretation", "final_evidence_status")]
}
pre <- do.call(rbind, pre_list)
rownames(pre) <- NULL
eligible <- pre[pre$final_evidence_status == "PRIMARY_VALIDATION", , drop = FALSE]
if (length(unique(eligible$dataset)) != 2L || length(unique(eligible$contrast)) != 4L) stop("Final primary table is not 2 datasets / 4 contrasts")
if (any(eligible$dataset == "GSE81465")) stop("GSE81465 entered final primary table")
write_tsv(pre, file.path(out_root, "results/predefined_validation/predefined_hallmark_validation.tsv"))

path_summary <- function(path) {
  z <- eligible[eligible$pathway == path, , drop = FALSE]
  pos_datasets <- length(unique(z$dataset[z$NES > 0]))
  pos_systems <- length(unique(z$same_experimental_system_id[z$NES > 0]))
  neg_systems <- length(unique(z$same_experimental_system_id[z$NES < 0]))
  frac <- mean(z$NES > 0); med <- median(z$NES)
  strong <- frac >= .70 && pos_systems >= 2 && med > .75 && any(z$NES > 0 & z$pval < .05)
  moderate <- frac > .50 && med > 0 && pos_datasets >= 2
  rating <- if (strong) "STRONGLY_VALIDATED" else if (moderate) "MODERATELY_VALIDATED" else if (frac >= .50 || med > 0) "WEAKLY_SUPPORTED" else "NOT_VALIDATED"
  data.frame(pathway = path, n_available = nrow(z), n_positive = sum(z$NES > 0), n_negative = sum(z$NES < 0),
    positive_fraction = frac, median_NES = med,
    number_nominal_significant_positive = sum(z$NES > 0 & z$pval < .05),
    number_FDR_significant_positive = sum(z$NES > 0 & z$padj < .05),
    number_independent_systems_positive = pos_systems, number_independent_systems_negative = neg_systems,
    number_datasets_positive = pos_datasets, number_datasets_total = length(unique(z$dataset)),
    validation_rating = rating, stringsAsFactors = FALSE)
}
mtor <- path_summary(frozen[1]); upr <- path_summary(frozen[2])
write_tsv(mtor, file.path(out_root, "results/predefined_validation/MTORC1_validation_summary.tsv"))
write_tsv(upr, file.path(out_root, "results/predefined_validation/UPR_validation_summary.tsv"))

div_list <- list()
for (p in frozen[3:7]) {
  z <- eligible[eligible$pathway == p, , drop = FALSE]
  sp <- length(unique(z$same_experimental_system_id[z$NES > 0])); sn <- length(unique(z$same_experimental_system_id[z$NES < 0]))
  rating <- if (sp >= 1 && sn >= 1 && max(sp, sn) >= 2) "DIVERGENCE_VALIDATED" else if (sp >= 1 || sn >= 1) "CONTEXT_DEPENDENCE_SUPPORTED" else "NOT_VALIDATED"
  div_list[[p]] <- data.frame(pathway = p, n_available = nrow(z), n_positive = sum(z$NES > 0), n_negative = sum(z$NES < 0),
    NES_min = min(z$NES), NES_max = max(z$NES), NES_range = diff(range(z$NES)), NES_SD = sd(z$NES),
    independent_systems_positive = sp, independent_systems_negative = sn,
    datasets_positive = length(unique(z$dataset[z$NES > 0])), datasets_negative = length(unique(z$dataset[z$NES < 0])),
    validation_rating = rating, role = unname(frozen_roles[p]), stringsAsFactors = FALSE)
}
divergence <- do.call(rbind, div_list)
write_tsv(divergence, file.path(out_root, "results/predefined_validation/divergent_program_validation_summary.tsv"))

theme_map <- data.frame(theme = c("MTORC1_related", "UPR_related", "IFN_related", "EMT_inflammation_related"),
  regex = c("MTOR|NUTRIENT|STARVATION|TRANSLATION|RIBOSOM|PROTEIN_SYNTHESIS", "UNFOLDED|ER_STRESS|PROTEIN_FOLDING",
            "INTERFERON", "EXTRACELLULAR_MATRIX|ECM_|COLLAGEN|INFLAMMAT|CYTOKINE|CHEMOKINE"), stringsAsFactors = FALSE)
write_tsv(theme_map, file.path(out_root, "results/predefined_validation/reactome_theme_mapping.tsv"))
react_list <- list()
for (id in primary_role$contrast) {
  g <- read.delim(file.path(r3, "results/gsea", paste0(id, "_reactome_full.tsv")), check.names = FALSE, stringsAsFactors = FALSE)
  g$evidence_tier <- role$evidence_tier[match(id, role$contrast)]
  for (i in seq_len(nrow(theme_map))) {
    q <- g[grepl(theme_map$regex[i], g$pathway), , drop = FALSE]
    if (nrow(q)) { q$predefined_theme <- theme_map$theme[i]; react_list[[paste(id, i)]] <- q }
  }
}
if (length(react_list)) write_tsv(do.call(rbind, react_list), file.path(out_root, "results/predefined_validation/predefined_reactome_support.tsv"))

primary <- unique(eligible[, c("pathway", "dataset", "contrast", "evidence_tier", "NES", "pval", "padj", "same_experimental_system_id")])
robust_summary <- function(z, label, excluded = "none") data.frame(analysis = label, excluded = excluded,
  pathway = unique(z$pathway), n_contrasts = nrow(z), n_datasets = length(unique(z$dataset)),
  n_systems = length(unique(z$same_experimental_system_id)), n_positive = sum(z$NES > 0),
  positive_fraction = mean(z$NES > 0), median_NES = median(z$NES), stringsAsFactors = FALSE)
lodo <- list(); loc <- list()
for (p in frozen[1:2]) {
  for (d in unique(primary$dataset)) lodo[[paste(p, d)]] <- robust_summary(primary[primary$pathway == p & primary$dataset != d, ], "leave_one_validation_dataset_out", d)
  for (id in unique(primary$contrast)) loc[[paste(p, id)]] <- robust_summary(primary[primary$pathway == p & primary$contrast != id, ], "leave_one_validation_contrast_out", id)
}
lodo <- do.call(rbind, lodo); loc <- do.call(rbind, loc)
smallest <- "GSE86525_HT29_bevacizumab_resistant_vs_control"
small <- do.call(rbind, lapply(frozen[1:2], function(p) robust_summary(primary[primary$pathway == p & primary$contrast != smallest, ], "exclude_smallest_validation_contrast", smallest)))
write_tsv(lodo, file.path(out_root, "results/robustness/leave_one_validation_dataset_out.tsv"))
write_tsv(loc, file.path(out_root, "results/robustness/leave_one_validation_contrast_out.tsv"))
write_tsv(small, file.path(out_root, "results/robustness/exclude_smallest_validation_contrast.tsv"))

disc <- read.delim(file.path(phase_root, "results/round2b/matrices/hallmark_core_NES_matrix.tsv"), check.names = FALSE)
disc_long <- do.call(rbind, lapply(seq_len(nrow(disc)), function(i) data.frame(pathway = disc$pathway[i], NES = as.numeric(disc[i, -1]), stringsAsFactors = FALSE)))
dv <- list()
for (p in frozen) {
  d <- disc_long[disc_long$pathway == p, ]; v <- primary[primary$pathway == p, ]
  dv[[p]] <- data.frame(pathway = p, discovery_n = nrow(d), discovery_positive = sum(d$NES > 0), discovery_positive_fraction = mean(d$NES > 0), discovery_median_NES = median(d$NES),
    validation_n = nrow(v), validation_positive = sum(v$NES > 0), validation_positive_fraction = mean(v$NES > 0), validation_median_NES = median(v$NES),
    validation_datasets = length(unique(v$dataset)), validation_systems = length(unique(v$same_experimental_system_id)), stringsAsFactors = FALSE)
}
write_tsv(do.call(rbind, dv), file.path(out_root, "results/predefined_validation/discovery_vs_validation_predefined_programs.tsv"))

sample_path <- file.path(r3, "data/metadata/validation_sample_manifest.tsv")
sample_manifest <- read.delim(sample_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_manifest$evidence_tier[sample_manifest$dataset == "GSE81465"] <- "Supportive-only"
sample_manifest$include_exclude[sample_manifest$dataset == "GSE81465"] <- "supportive_descriptive"
for (i in seq_len(nrow(role))) {
  hit <- sample_manifest$proposed_contrast == role$contrast[i]
  sample_manifest$evidence_tier[hit] <- role$evidence_tier[i]
  sample_manifest$include_exclude[hit] <- ifelse(role$final_evidence_status[i] == "PRIMARY_VALIDATION", "include_primary", "supportive_descriptive")
}
write_tsv(sample_manifest, file.path(out_root, "data/metadata/validation_sample_manifest.tsv"))

readiness_path <- file.path(r3, "data/metadata/validation_contrast_readiness.tsv")
readiness <- read.delim(readiness_path, check.names = FALSE, stringsAsFactors = FALSE)
for (i in seq_len(nrow(role))) {
  hit <- readiness$contrast_id == role$contrast[i]
  readiness$evidence_tier[hit] <- role$evidence_tier[i]
  readiness$eligible_primary_denominator[hit] <- role$final_evidence_status[i] == "PRIMARY_VALIDATION"
  readiness$readiness[hit] <- ifelse(role$final_evidence_status[i] == "PRIMARY_VALIDATION", "READY_WITH_CAUTION", "SUPPORTIVE_ONLY")
  if (role$dataset[i] == "GSE81465") readiness$caution[hit] <- "Replicate provenance unresolved; supportive serial-state evidence only"
}
write_tsv(readiness, file.path(out_root, "data/metadata/validation_contrast_readiness.tsv"))

decision <- data.frame(primary_datasets = 2, primary_contrasts = 4, primary_systems = 4,
  MTORC1_rating = mtor$validation_rating, UPR_rating = upr$validation_rating,
  divergence_validated_primary_paths = sum(divergence$pathway %in% frozen[3:6] & divergence$validation_rating == "DIVERGENCE_VALIDATED"),
  validation_model = "Final four-contrast validation", final_status = "FINAL_ROLE_ASSIGNMENT",
  reason = "Primary validation restricted to GSE64052 and GSE86525; GSE81465 supportive-only; nesting leverage retained",
  stop = "ANALYSIS_COMPLETE", stringsAsFactors = FALSE)
write_tsv(decision, file.path(out_root, "reports/round3a_final_decision.tsv"))

assertions <- data.frame(
  assertion = c("primary_dataset_count == 2", "primary_contrast_count == 4", "GSE81465 not in PRIMARY_VALIDATION", "GSE81465 role == SUPPORTIVE_ONLY"),
  observed = c(length(unique(eligible$dataset)), length(unique(eligible$contrast)), !any(eligible$dataset == "GSE81465"), g814_role$final_evidence_status),
  status = "PASS", stringsAsFactors = FALSE)
write_tsv(assertions, file.path(out_root, "reports/round3a_role_assertions.tsv"))
message("Round3A role-dependent summary rebuild complete: 2 datasets / 4 contrasts; GSE81465 SUPPORTIVE_ONLY")
