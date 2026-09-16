#!/usr/bin/env Rscript

# External-model validation of the frozen discovery programs.
# Discovery matrices are read-only and never enter the validation denominator.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
phase_root <- normalizePath(analysis_root_env, mustWork = TRUE)
r3 <- file.path(phase_root, "round3a_validation")
suppressPackageStartupMessages({library(limma); library(fgsea); library(msigdbr)})

dirs <- c("data/processed", "data/metadata", "results/qc", "results/differential",
          "results/gsea", "results/predefined_validation", "results/temporal",
          "results/robustness", "results/figures", "reports")
for (d in dirs) dir.create(file.path(r3, d), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE,
                                        quote = FALSE, na = "NA")
read_soft_matrix <- function(path) {
  z <- readLines(gzfile(path), warn = FALSE)
  a <- grep("^!series_matrix_table_begin", z) + 1L
  b <- grep("^!series_matrix_table_end", z) - 1L
  read.delim(textConnection(z[a:b]), check.names = FALSE, stringsAsFactors = FALSE)
}
read_soft_platform <- function(path) {
  z <- readLines(path, warn = FALSE)
  a <- grep("^!platform_table_begin", z) + 1L
  b <- grep("^!platform_table_end", z) - 1L
  read.delim(textConnection(z[a:b]), check.names = FALSE, stringsAsFactors = FALSE,
             quote = "", comment.char = "")
}
read_geo_annot <- function(path) {
  z <- readLines(gzfile(path), warn = FALSE)
  a <- grep("^ID\t", z)[1]
  read.delim(textConnection(z[a:length(z)]), check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "---", "NA", "null")] <- NA_character_
  x[grepl("///|//|;|,", x)] <- NA_character_
  toupper(x)
}
collapse_gene <- function(expr, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols) & apply(expr, 1, function(v) all(is.finite(v)))
  expr <- expr[keep, , drop = FALSE]; symbols <- symbols[keep]
  iq <- apply(expr, 1, IQR, na.rm = TRUE)
  ord <- order(symbols, -iq)
  expr <- expr[ord, , drop = FALSE]; symbols <- symbols[ord]
  first <- !duplicated(symbols)
  out <- expr[first, , drop = FALSE]
  rownames(out) <- symbols[first]
  out
}
sample_cor_min <- function(x) min(cor(x, use = "pairwise.complete.obs")[upper.tri(cor(x))], na.rm = TRUE)

# The project-level role table is the single source for every primary/supportive
# decision in Round 3A. Scientific values are not encoded in this table.
role_config_path <- normalizePath(file.path(dirname(script_path), "..", "config", "validation_contrast_roles.tsv"), mustWork = TRUE)
if (!file.exists(role_config_path)) stop("Missing final validation role config: ", role_config_path)
validation_roles <- read.delim(role_config_path, check.names = FALSE, stringsAsFactors = FALSE)
required_role_fields <- c("dataset", "contrast", "final_evidence_status", "primary_validation_eligible", "evidence_tier")
if (!all(required_role_fields %in% names(validation_roles))) stop("Validation role config has missing required fields")
if (anyDuplicated(validation_roles$contrast)) stop("Validation role config contains duplicate contrast definitions")
role_row <- function(contrast) {
  z <- validation_roles[validation_roles$contrast == contrast, , drop = FALSE]
  if (nrow(z) != 1L) stop("Expected exactly one role definition for ", contrast)
  z
}
role_tier <- function(contrast) role_row(contrast)$evidence_tier[[1]]
role_status <- function(contrast) role_row(contrast)$final_evidence_status[[1]]
role_primary <- function(contrast) identical(role_status(contrast), "PRIMARY_VALIDATION")

# -----------------------------------------------------------------------------
# Frozen manifest: manual interpretation after GEO/GSM/publication cross-check.
# -----------------------------------------------------------------------------
manifest_rows <- list()
add_manifest <- function(dataset, gsm, title, organism, platform, cancer_type, model,
                         drug, drug_class, treatment_stage, resistance_state, compartment,
                         biological_replicate, technical_replicate, proposed_contrast,
                         comparator_type, evidence_tier, include_exclude, exclusion_reason,
                         notes, same_experimental_system_id) {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    dataset, GSM = gsm, title, organism, platform, cancer_type, model, drug, drug_class,
    treatment_stage, resistance_state, compartment, biological_replicate,
    technical_replicate, proposed_contrast, comparator_type, evidence_tier,
    include_exclude, exclusion_reason, notes, same_experimental_system_id,
    stringsAsFactors = FALSE)
}

g814 <- c("GSM2208762","GSM2208763","GSM2208764","GSM2208765","GSM2208766","GSM2208767",
          "GSM2208768","GSM2208769","GSM2208770","GSM2208771","GSM2208772","GSM2208773")
t814 <- c("Glioblastoma xenograft - bevacizumab, replicate 1","Glioblastoma xenograft - bevacizumab, replicate 2","Glioblastoma xenograft - bevacizumab, replicate 3",
          "Glioblastoma xenograft + bevacizumab at one generation, replicate 1","Glioblastoma xenograft + bevacizumab at one generation, replicate 2","Glioblastoma xenograft + bevacizumab at one generation, replicate 3",
          "Glioblastoma xenograft + bevacizumab at four generations, replicate 1","Glioblastoma xenograft + bevacizumab at four generations, replicate 2","Glioblastoma xenograft + bevacizumab at four generations, replicate 3",
          "Glioblastoma xenograft + bevacizumab at nine generations, replicate 1","Glioblastoma xenograft + bevacizumab at nine generations, replicate 2","Glioblastoma xenograft + bevacizumab at nine generations, replicate 3")
stage814 <- rep(c("IgG control", "generation 1 / early bevacizumab", "generation 4 / intermediate", "generation 9 / late resistant"), each = 3)
state814 <- rep(c("control", "early responsive/sensitive", "intermediate adaptation", "late resistant"), each = 3)
for (i in seq_along(g814)) add_manifest("GSE81465", g814[i], t814[i], "Homo sapiens", "GPL10558",
  "glioblastoma", "U87-derived serial generational xenograft", "bevacizumab", "VEGF ligand blockade",
  stage814[i], state814[i], "bulk xenograft; human-array tumor-enriched signal",
  "GEO calls independent biological; publication Methods calls n=3 technical replicates/group",
  paste0("replicate ", ((i-1) %% 3)+1),
  if (i >= 10) "GSE81465_G9_vs_G1" else if (i >= 4 && i <= 6) "GSE81465_G9_vs_G1" else "GSE81465_temporal_only",
  if (i >= 4) "longitudinal resistance-development stage" else "IgG control", role_tier("GSE81465_G9_vs_G1"),
  "supportive_descriptive", "",
  "Resistance stage supported by serial passage design; replicate-type conflict retained as a limitation.",
  "GSE81465_U87_bevacizumab_generational")

g640 <- paste0("GSM", 1563498:1563525)
t640 <- c("33B3_Untreated_786-O","33B9_Untreated_786-O","33C6_Untreated_786-O","33D3_Untreated_786-O","33D9_Untreated_786-O",
          "33B6_Resistant_786-O","33C3_Resistant_786-O","33C9_Resistant_786-O","33D6_Resistant_786-O","33E3_Resistant_786-O","35B7_Resistant_A498",
          "35A1_Untreated_786-O","35A2_Untreated_786-O","35A3_Untreated_786-O","35A4_Untreated_786-O","35A5_Untreated_786-O",
          "35D4_Resistant_786-O","35D5_Resistant_786-O","35D7_Resistant_786-O","35D8_Resistant_786-O",
          "35A6_Untreated_A498","35A8_Untreated_A498","35A9_Untreated_A498","35B1_Untreated_A498",
          "35B8_Resistant_A498","35B9_Resistant_A498","35C1_Resistant_A498","35C2_Resistant_A498")
drug640 <- c(rep("sorafenib", 11), rep("sunitinib", 9), rep("sorafenib", 8))
model640 <- ifelse(grepl("A498", t640), "A498 xenograft", "786-O xenograft")
state640 <- ifelse(grepl("Resistant", t640), "resistant during continued therapy", "untreated baseline")
sys640 <- ifelse(model640 == "A498 xenograft", "GSE64052_A498_sorafenib",
                 ifelse(drug640 == "sunitinib", "GSE64052_786O_sunitinib", "GSE64052_786O_sorafenib"))
cid640 <- paste0(sys640, "_resistant_vs_untreated")
for (i in seq_along(g640)) add_manifest("GSE64052", g640[i], t640[i], "Homo sapiens", "GPL570",
  "clear-cell renal cell carcinoma", model640[i], drug640[i], "multikinase VEGFR-targeting TKI",
  ifelse(state640[i] == "untreated baseline", "pretreatment 12-mm tumor", "resistant endpoint under continued drug"),
  state640[i], "bulk xenograft; human-array tumor-enriched signal", paste0("mouse tumor ", sub("_.*", "", t640[i])), "not reported",
  cid640[i], "untreated baseline", role_tier(cid640[i]), ifelse(role_primary(cid640[i]), "include_primary", "supportive_descriptive"), "",
  "Resistance was growth despite therapy; untreated vs resistant retains exposure/time/tumor-size confounding.", sys640[i])

g865 <- paste0("GSM230499", 2:7)
t865 <- c("HT-29_Control_rep1","HT-29_Control_rep2","HT-29_Control_rep3",
          "HT-29_Bevacizumab-resistant rep1","HT-29_Bevacizumab-resistant rep2","HT-29_Bevacizumab-resistant rep3")
for (i in seq_along(g865)) add_manifest("GSE86525", g865[i], t865[i], "Homo sapiens", "GPL16699",
  "colorectal adenocarcinoma", "HT-29 subcutaneous xenograft", "bevacizumab", "VEGF ligand blockade",
  if (i <= 3) "untreated control" else paste0("continued treatment to endpoint; ", c("3 weeks","4 weeks","5 weeks")[i-3]),
  if (i <= 3) "control" else "bevacizumab-treated resistant/escape endpoint",
  "bulk xenograft; human-array tumor-enriched signal", paste0("replicate ", ((i-1) %% 3)+1), "not reported",
  "GSE86525_HT29_bevacizumab_resistant_vs_control", "untreated control", role_tier("GSE86525_HT29_bevacizumab_resistant_vs_control"), ifelse(role_primary("GSE86525_HT29_bevacizumab_resistant_vs_control"), "include_primary", "supportive_descriptive"), "",
  if (i <= 3) "Non-treated." else paste0("Treatment duration differs by endpoint: ", c("3","4","5")[i-3], " weeks; 5 mg/kg twice weekly."),
  "GSE86525_HT29_bevacizumab")

g451 <- paste0("GSM10984", sprintf("%02d", 4:13))
t451 <- c("NSC11 control1","NSC11 control2","NSC11 Avastin1","NSC11 Avastin2","NSC11R Control","NSC11R Avastin",
          "U87 Control","U87 Avastin","U87R Control","U87R Avastin")
model451 <- ifelse(grepl("NSC11", t451), "NSC11 orthotopic xenograft", "U87 orthotopic xenograft")
drug451 <- ifelse(grepl("Avastin", t451), "bevacizumab", "none at harvest")
state451 <- ifelse(grepl("R ", t451), "resistant derivative", "parental")
cid451 <- ifelse(grepl("NSC11", t451), ifelse(grepl("Avastin", t451), "GSE45161_NSC11R_vs_parental_bevacizumab", "GSE45161_NSC11R_vs_parental_control"),
                 ifelse(grepl("Avastin", t451), "GSE45161_U87R_vs_parental_bevacizumab", "GSE45161_U87R_vs_parental_control"))
for (i in seq_along(g451)) add_manifest("GSE45161", g451[i], t451[i], "Homo sapiens", "GPL9324",
  "glioblastoma", model451[i], drug451[i], "VEGF ligand blockade", ifelse(drug451[i] == "bevacizumab", "on-treatment endpoint", "untreated endpoint"),
  state451[i], "bulk orthotopic xenograft; custom human Entrez CDF", ifelse(i <= 4, paste0("replicate ", ((i-1) %% 2)+1), "single sample"),
  "not reported", cid451[i], "matched parental state", "Supportive-only", "supportive_only",
  "Resistant derivative arm has n=1 per treatment state; excluded from primary denominator.",
  "Publication confirms acquired-resistant derivatives, but GEO contains only one resistant-derivative sample per treatment state.",
  paste0("GSE45161_", ifelse(grepl("NSC11", t451[i]), "NSC11", "U87"), "_derivative"))

manifest <- do.call(rbind, manifest_rows)
write_tsv(manifest, file.path(r3, "data/metadata/validation_sample_manifest.tsv"))

final_role_contrasts <- c("GSE81465_G9_vs_G1",
    "GSE64052_786O_sorafenib_resistant_vs_untreated",
    "GSE64052_786O_sunitinib_resistant_vs_untreated",
    "GSE64052_A498_sorafenib_resistant_vs_untreated",
    "GSE86525_HT29_bevacizumab_resistant_vs_control")
final_role_primary <- vapply(final_role_contrasts, role_primary, logical(1))
readiness <- data.frame(
  contrast_id = c(final_role_contrasts,
    "GSE45161_NSC11R_vs_parental_control","GSE45161_NSC11R_vs_parental_bevacizumab",
    "GSE45161_U87R_vs_parental_control","GSE45161_U87R_vs_parental_bevacizumab"),
  dataset = c("GSE81465",rep("GSE64052",3),"GSE86525",rep("GSE45161",4)),
  same_experimental_system_id = c("GSE81465_U87_bevacizumab_generational",
    "GSE64052_786O_sorafenib","GSE64052_786O_sunitinib","GSE64052_A498_sorafenib",
    "GSE86525_HT29_bevacizumab",rep("GSE45161_NSC11_derivative",2),rep("GSE45161_U87_derivative",2)),
  evidence_tier = c(vapply(final_role_contrasts, role_tier, character(1)), rep("Supportive-only",4)),
  n_resistant = c(3,5,4,5,3,1,1,1,1), n_comparator = c(3,5,5,4,3,2,2,1,1),
  eligible_primary_denominator = c(final_role_primary,rep(FALSE,4)),
  readiness = c(ifelse(final_role_primary,"READY_WITH_CAUTION","SUPPORTIVE_ONLY"),rep("SUPPORTIVE_DIRECTIONAL_ONLY",4)),
  caution = c("GEO calls biological replicates but publication Methods calls technical replicates",
    rep("Tier B: resistance is confounded with exposure, elapsed time, and tumor size",3),
    "Tier B: treated endpoint vs untreated; treated replicates reached endpoint at 3, 4, and 5 weeks",
    rep("Resistant derivative arm n=1; no inferential claim",4)), stringsAsFactors = FALSE)
write_tsv(readiness, file.path(r3, "data/metadata/validation_contrast_readiness.tsv"))

# -----------------------------------------------------------------------------
# Independent preprocessing. No matrices are merged or batch-corrected.
# -----------------------------------------------------------------------------
raw <- file.path(r3, "data/raw")
exprs <- list()

# GSE81465: only dataset without an official processed matrix; neqc from IDAT.
idat_dir <- file.path(raw, "GSE81465")
idats <- list.files(idat_dir, "idat$", full.names = TRUE)
x814 <- read.idat(idats, file.path(idat_dir, "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"), verbose = FALSE)
x814 <- neqc(x814)
colnames(x814$E) <- sub("_.*", "", basename(idats))
e814 <- collapse_gene(x814$E, clean_symbol(x814$genes$Symbol))
exprs$GSE81465 <- e814

# GSE64052: official processed MAS5 signal, log2 plus within-study quantile normalization.
m640 <- read_soft_matrix(file.path(raw, "GSE64052/GSE64052_series_matrix.txt.gz"))
rownames(m640) <- m640[[1]]; e640p <- as.matrix(m640[-1]); storage.mode(e640p) <- "double"
e640p <- normalizeBetweenArrays(log2(e640p + 1), method = "quantile")
a570 <- read_geo_annot(file.path(raw, "GPL570.annot.gz"))
sym570 <- clean_symbol(a570$`Gene symbol`[match(rownames(e640p), a570$ID)])
e640 <- collapse_gene(e640p, sym570); exprs$GSE64052 <- e640

# GSE86525: official processed background-corrected signal, log2 plus quantile normalization.
m865 <- read_soft_matrix(file.path(raw, "GSE86525/GSE86525_series_matrix.txt.gz"))
rownames(m865) <- m865[[1]]; e865p <- as.matrix(m865[-1]); storage.mode(e865p) <- "double"
e865p <- normalizeBetweenArrays(log2(e865p + 1), method = "quantile")
a16699 <- read_soft_platform(file.path(raw, "GPL16699_full.soft"))
sym16699 <- clean_symbol(a16699$GENE_SYMBOL[match(rownames(e865p), a16699$ID)])
e865 <- collapse_gene(e865p, sym16699); exprs$GSE86525 <- e865

# GSE45161: official log2, median-normalized custom Entrez CDF matrix.
m451 <- read_soft_matrix(file.path(raw, "GSE45161/GSE45161_series_matrix.txt.gz"))
rownames(m451) <- m451[[1]]; e451p <- as.matrix(m451[-1]); storage.mode(e451p) <- "double"
a9324 <- read_soft_platform(file.path(raw, "GPL9324_full.soft"))
id451 <- sub("_at$", "", rownames(e451p))
sym451 <- clean_symbol(a9324$ORF[match(id451, a9324$EntrezID)])
e451 <- collapse_gene(e451p, sym451); exprs$GSE45161 <- e451

for (d in names(exprs)) {
  write_tsv(data.frame(gene_symbol = rownames(exprs[[d]]), exprs[[d]], check.names = FALSE),
            file.path(r3, "data/processed", paste0(d, "_gene_expression.tsv")))
}

# -----------------------------------------------------------------------------
# Per-study QC and figures.
# -----------------------------------------------------------------------------
group_map <- list(
  GSE81465 = setNames(rep(c("IgG","G1","G4","G9"), each=3), g814),
  GSE64052 = setNames(paste(model640, drug640, ifelse(grepl("Resistant",t640),"resistant","untreated"),sep="|"), g640),
  GSE86525 = setNames(ifelse(seq_along(g865)<=3,"control","resistant"),g865),
  GSE45161 = setNames(t451,g451))
qc_rows <- list()
for (d in names(exprs)) {
  e <- exprs[[d]]; gr <- unname(group_map[[d]][colnames(e)])
  pc <- prcomp(t(e), scale. = FALSE)
  pv <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  cors <- cor(e, use = "pairwise.complete.obs")
  sample_qc <- data.frame(dataset=d, sample=colnames(e), group=gr,
    missing_fraction=colMeans(!is.finite(e)), median_expression=apply(e,2,median),
    IQR_expression=apply(e,2,IQR), PC1=pc$x[,1], PC2=pc$x[,2],
    mean_correlation=rowMeans(cors-diag(ncol(e)),na.rm=TRUE)*ncol(e)/(ncol(e)-1), stringsAsFactors=FALSE)
  write_tsv(sample_qc,file.path(r3,"results/qc",paste0(d,"_sample_qc.tsv")))
  q <- data.frame(dataset=d,n_samples=ncol(e),n_genes=nrow(e),missing_fraction=mean(!is.finite(e)),
    min_pairwise_correlation=sample_cor_min(e),PC1_variance_percent=pv[1],PC2_variance_percent=pv[2],
    outlier_flag=ifelse(any(sample_qc$mean_correlation < median(sample_qc$mean_correlation)-3*mad(sample_qc$mean_correlation)),"REVIEW","NONE"),
    preprocessing=if (d=="GSE81465") "raw IDAT neqc" else if(d=="GSE45161") "official processed log2 custom-CDF" else "official processed signal; log2+quantile within study",
    stringsAsFactors=FALSE)
  qc_rows[[d]] <- q; write_tsv(q,file.path(r3,"results/qc",paste0(d,"_qc_summary.tsv")))
  png(file.path(r3,"results/qc",paste0(d,"_distributions.png")),1600,1100,res=170)
  boxplot(e,las=2,outline=FALSE,main=paste(d,"expression distributions"),ylab="normalized log2 expression",cex.axis=.65); dev.off()
  png(file.path(r3,"results/qc",paste0(d,"_PCA.png")),1400,1100,res=170)
  cols <- as.integer(factor(gr)); plot(pc$x[,1],pc$x[,2],pch=19,col=cols,xlab=sprintf("PC1 (%.1f%%)",pv[1]),ylab=sprintf("PC2 (%.1f%%)",pv[2]),main=paste(d,"PCA"));
  text(pc$x[,1],pc$x[,2],labels=colnames(e),pos=3,cex=.58); legend("topright",legend=levels(factor(gr)),col=seq_along(levels(factor(gr))),pch=19,cex=.55); dev.off()
  png(file.path(r3,"results/qc",paste0(d,"_correlation_heatmap.png")),1500,1400,res=170)
  heatmap(cors,Rowv=NA,Colv=NA,scale="none",col=colorRampPalette(c("#2166AC","white","#B2182B"))(101),margins=c(10,10),main=paste(d,"sample correlation")); dev.off()
  png(file.path(r3,"results/qc",paste0(d,"_clustering.png")),1600,1000,res=170)
  plot(hclust(as.dist(1-cors)),main=paste(d,"hierarchical clustering"),xlab="",sub=""); dev.off()
}
write_tsv(do.call(rbind,qc_rows),file.path(r3,"results/qc/all_validation_dataset_qc_summary.tsv"))

# -----------------------------------------------------------------------------
# Differential analysis: primary contrasts and supportive directional contrasts.
# -----------------------------------------------------------------------------
contrast_specs <- list(
  GSE81465_G9_vs_G1=list(dataset="GSE81465",num=g814[10:12],den=g814[4:6],drug="bevacizumab",cancer="glioblastoma",tier=role_tier("GSE81465_G9_vs_G1"),system="GSE81465_U87_bevacizumab_generational",primary=role_primary("GSE81465_G9_vs_G1")),
  GSE64052_786O_sorafenib_resistant_vs_untreated=list(dataset="GSE64052",num=g640[6:10],den=g640[1:5],drug="sorafenib",cancer="renal cell carcinoma",tier=role_tier("GSE64052_786O_sorafenib_resistant_vs_untreated"),system="GSE64052_786O_sorafenib",primary=role_primary("GSE64052_786O_sorafenib_resistant_vs_untreated")),
  GSE64052_786O_sunitinib_resistant_vs_untreated=list(dataset="GSE64052",num=g640[17:20],den=g640[12:16],drug="sunitinib",cancer="renal cell carcinoma",tier=role_tier("GSE64052_786O_sunitinib_resistant_vs_untreated"),system="GSE64052_786O_sunitinib",primary=role_primary("GSE64052_786O_sunitinib_resistant_vs_untreated")),
  GSE64052_A498_sorafenib_resistant_vs_untreated=list(dataset="GSE64052",num=g640[c(11,25:28)],den=g640[21:24],drug="sorafenib",cancer="renal cell carcinoma",tier=role_tier("GSE64052_A498_sorafenib_resistant_vs_untreated"),system="GSE64052_A498_sorafenib",primary=role_primary("GSE64052_A498_sorafenib_resistant_vs_untreated")),
  GSE86525_HT29_bevacizumab_resistant_vs_control=list(dataset="GSE86525",num=g865[4:6],den=g865[1:3],drug="bevacizumab",cancer="colorectal cancer",tier=role_tier("GSE86525_HT29_bevacizumab_resistant_vs_control"),system="GSE86525_HT29_bevacizumab",primary=role_primary("GSE86525_HT29_bevacizumab_resistant_vs_control")),
  GSE45161_NSC11R_vs_parental_control=list(dataset="GSE45161",num=g451[5],den=g451[1:2],drug="none at harvest",cancer="glioblastoma",tier="Supportive-only",system="GSE45161_NSC11_derivative",primary=FALSE),
  GSE45161_NSC11R_vs_parental_bevacizumab=list(dataset="GSE45161",num=g451[6],den=g451[3:4],drug="bevacizumab",cancer="glioblastoma",tier="Supportive-only",system="GSE45161_NSC11_derivative",primary=FALSE),
  GSE45161_U87R_vs_parental_control=list(dataset="GSE45161",num=g451[9],den=g451[7],drug="none at harvest",cancer="glioblastoma",tier="Supportive-only",system="GSE45161_U87_derivative",primary=FALSE),
  GSE45161_U87R_vs_parental_bevacizumab=list(dataset="GSE45161",num=g451[10],den=g451[8],drug="bevacizumab",cancer="glioblastoma",tier="Supportive-only",system="GSE45161_U87_derivative",primary=FALSE))

for (id in intersect(names(contrast_specs), validation_roles$contrast)) {
  contrast_specs[[id]]$tier <- role_tier(id)
  contrast_specs[[id]]$primary <- role_primary(id)
}
primary_ids <- names(contrast_specs)[vapply(contrast_specs, function(s) isTRUE(s$primary), logical(1))]
if (length(primary_ids) != 4L) stop("Final primary contrast assertion failed: expected 4")
if (length(unique(vapply(contrast_specs[primary_ids], function(s) s$dataset, character(1)))) != 2L) stop("Final primary dataset assertion failed: expected 2")
if ("GSE81465_G9_vs_G1" %in% primary_ids) stop("GSE81465 entered PRIMARY_VALIDATION")
if (!identical(role_status("GSE81465_G9_vs_G1"), "SUPPORTIVE_ONLY") || !identical(role_tier("GSE81465_G9_vs_G1"), "Supportive-only")) stop("GSE81465 supportive-only assertion failed")

rank_cache <- list(); de_summary <- list()
for (id in names(contrast_specs)) {
  s <- contrast_specs[[id]]; e <- exprs[[s$dataset]]; num <- intersect(s$num,colnames(e)); den <- intersect(s$den,colnames(e));
  sube <- e[,c(den,num),drop=FALSE]; grp <- factor(c(rep("den",length(den)),rep("num",length(num))),levels=c("den","num"))
  if (length(num)>=2 && length(den)>=2) {
    design <- model.matrix(~grp); fit <- eBayes(lmFit(sube,design),trend=TRUE,robust=TRUE)
    tt <- topTable(fit,coef="grpnum",number=Inf,sort.by="none")
    stat <- data.frame(gene_symbol=rownames(tt),logFC=tt$logFC,AveExpr=tt$AveExpr,t=tt$t,P.Value=tt$P.Value,adj.P.Val=tt$adj.P.Val,B=tt$B,stringsAsFactors=FALSE)
    ranks <- setNames(stat$t,stat$gene_symbol)
  } else {
    lf <- rowMeans(sube[,grp=="num",drop=FALSE])-rowMeans(sube[,grp=="den",drop=FALSE])
    stat <- data.frame(gene_symbol=rownames(sube),logFC=lf,AveExpr=rowMeans(sube),t=lf,P.Value=NA_real_,adj.P.Val=NA_real_,B=NA_real_,stringsAsFactors=FALSE)
    ranks <- setNames(lf,rownames(sube))
  }
  stat$dataset <- s$dataset; stat$contrast <- id; stat$drug <- s$drug; stat$cancer_type <- s$cancer; stat$evidence_tier <- s$tier; stat$same_experimental_system_id <- s$system
  stat <- stat[order(stat$t,decreasing=TRUE),]
  write_tsv(stat,file.path(r3,"results/differential",paste0(id,"_gene_statistics.tsv")))
  ranks <- ranks[is.finite(ranks)]; ranks <- sort(ranks[!duplicated(names(ranks))],decreasing=TRUE); rank_cache[[id]] <- ranks
  de_summary[[id]] <- data.frame(contrast_id=id,dataset=s$dataset,evidence_tier=s$tier,primary_denominator=s$primary,
    n_resistant=length(num),n_comparator=length(den),n_ranked_genes=length(ranks),positive_statistic="resistant/later-resistant higher",stringsAsFactors=FALSE)
}
write_tsv(do.call(rbind,de_summary),file.path(r3,"results/differential/differential_analysis_summary.tsv"))

# Additional GSE81465 temporal contrasts are descriptive and excluded from denominator.
temporal_specs <- list(GSE81465_G1_vs_IgG=list(num=g814[4:6],den=g814[1:3]),
                       GSE81465_G4_vs_IgG=list(num=g814[7:9],den=g814[1:3]),
                       GSE81465_G9_vs_IgG=list(num=g814[10:12],den=g814[1:3]))
for (id in names(temporal_specs)) {
  s <- temporal_specs[[id]]; sube <- e814[,c(s$den,s$num)]; grp <- factor(rep(c("den","num"),each=3),levels=c("den","num"))
  fit <- eBayes(lmFit(sube,model.matrix(~grp)),trend=TRUE,robust=TRUE); tt <- topTable(fit,coef="grpnum",number=Inf,sort.by="none")
  rank_cache[[id]] <- sort(setNames(tt$t,rownames(tt)),decreasing=TRUE)
}

# -----------------------------------------------------------------------------
# Ranked GSEA using the same MSigDB and engine as Round 2B.
# -----------------------------------------------------------------------------
hall_db <- msigdbr(db_species="HS",species="human",collection="H")
react_db <- msigdbr(db_species="HS",species="human",collection="C2",subcollection="CP:REACTOME")
pathways <- list(hallmark=split(hall_db$gene_symbol,hall_db$gs_name),reactome=split(react_db$gene_symbol,react_db$gs_name))
gsea_cache <- list(); gsea_summary <- list()
for (id in names(rank_cache)) for (coll in names(pathways)) {
  set.seed(20260909); fg <- as.data.frame(fgseaMultilevel(pathways[[coll]],rank_cache[[id]],minSize=15,maxSize=500,eps=0))
  fg$leadingEdge <- vapply(fg$leadingEdge,paste,collapse=";",FUN.VALUE=character(1)); fg <- fg[order(fg$padj,-abs(fg$NES)),]
  if (id %in% names(contrast_specs)) {
    s <- contrast_specs[[id]]; fg$dataset <- s$dataset; fg$drug <- s$drug; fg$cancer_type <- s$cancer; fg$evidence_tier <- s$tier; fg$same_experimental_system_id <- s$system
  } else { fg$dataset <- "GSE81465"; fg$drug <- "bevacizumab"; fg$cancer_type <- "glioblastoma"; fg$evidence_tier <- "Temporal-only"; fg$same_experimental_system_id <- "GSE81465_U87_bevacizumab_generational" }
  fg$contrast <- id; fg <- fg[,c("pathway","dataset","contrast","drug","cancer_type","evidence_tier","same_experimental_system_id","NES","pval","padj","log2err","ES","size","leadingEdge")]
  write_tsv(fg,file.path(r3,"results/gsea",paste0(id,"_",coll,"_full.tsv")))
  gsea_cache[[paste(id,coll,sep="|")]] <- fg
  gsea_summary[[paste(id,coll)]] <- data.frame(contrast=id,collection=coll,tested=nrow(fg),nominal_p_lt_0_05=sum(fg$pval<.05),FDR_lt_0_05=sum(fg$padj<.05),stringsAsFactors=FALSE)
}
write_tsv(do.call(rbind,gsea_summary),file.path(r3,"results/gsea/gsea_run_summary.tsv"))
write_tsv(data.frame(collection=c("Hallmark","Reactome"),msigdb_version=c(unique(hall_db$db_version),unique(react_db$db_version)),
  gene_sets=c(length(pathways$hallmark),length(pathways$reactome)),engine="fgseaMultilevel",min_size=15,max_size=500,rank="limma moderated t; supportive n=1 contrasts use logFC",stringsAsFactors=FALSE),
  file.path(r3,"reports/round3a_gsea_gene_set_manifest.tsv"))

frozen <- c("HALLMARK_MTORC1_SIGNALING","HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_E2F_TARGETS")
roles <- setNames(c("PRIMARY_SHARED_CANDIDATE","SECONDARY_SHARED_CANDIDATE",rep("PRIMARY_DIVERGENT_CANDIDATE",4),"CAUTIONARY_DIVERGENT_CANDIDATE"),frozen)
pre <- list()
for (id in names(contrast_specs)) {
  z <- gsea_cache[[paste(id,"hallmark",sep="|")]]; z <- z[match(frozen,z$pathway),]
  z$direction <- ifelse(z$NES>0,"positive","negative"); z$predefined_role <- unname(roles[z$pathway])
  final_status <- if (id %in% validation_roles$contrast) role_status(id) else "SUPPORTIVE_ONLY"
  z$validation_interpretation <- ifelse(final_status == "PRIMARY_VALIDATION","eligible primary validation contrast","supportive descriptive evidence; excluded from primary denominator")
  z$final_evidence_status <- final_status
  pre[[id]] <- z[,c("pathway","dataset","contrast","drug","cancer_type","evidence_tier","same_experimental_system_id","NES","pval","padj","direction","predefined_role","validation_interpretation","final_evidence_status")]
}
pre <- do.call(rbind,pre); rownames(pre) <- NULL
write_tsv(pre,file.path(r3,"results/predefined_validation/predefined_hallmark_validation.tsv"))
eligible <- pre[pre$final_evidence_status == "PRIMARY_VALIDATION",]
if (length(unique(eligible$dataset)) != 2L || length(unique(eligible$contrast)) != 4L) stop("Final primary validation universe is not 2 datasets / 4 contrasts")
if (any(eligible$dataset == "GSE81465")) stop("GSE81465 entered the eligible primary table")

path_summary <- function(path, shared=TRUE) {
  z <- eligible[eligible$pathway==path,]
  pos_datasets <- length(unique(z$dataset[z$NES>0])); pos_systems <- length(unique(z$same_experimental_system_id[z$NES>0]))
  neg_systems <- length(unique(z$same_experimental_system_id[z$NES<0])); frac <- mean(z$NES>0); med <- median(z$NES)
  strong <- frac>=.70 && pos_systems>=2 && med>.75 && any(z$NES>0 & z$pval<.05)
  moderate <- frac>.50 && med>0 && pos_datasets>=2
  rating <- if(strong) "STRONGLY_VALIDATED" else if(moderate) "MODERATELY_VALIDATED" else if(frac>=.50 || med>0) "WEAKLY_SUPPORTED" else "NOT_VALIDATED"
  data.frame(pathway=path,n_available=nrow(z),n_positive=sum(z$NES>0),n_negative=sum(z$NES<0),positive_fraction=frac,median_NES=med,
    number_nominal_significant_positive=sum(z$NES>0 & z$pval<.05),number_FDR_significant_positive=sum(z$NES>0 & z$padj<.05),
    number_independent_systems_positive=pos_systems,number_independent_systems_negative=neg_systems,
    number_datasets_positive=pos_datasets,number_datasets_total=length(unique(z$dataset)),validation_rating=rating,stringsAsFactors=FALSE)
}
mtor <- path_summary(frozen[1]); upr <- path_summary(frozen[2])
write_tsv(mtor,file.path(r3,"results/predefined_validation/MTORC1_validation_summary.tsv"))
write_tsv(upr,file.path(r3,"results/predefined_validation/UPR_validation_summary.tsv"))

div_rows <- list()
for (p in frozen[3:7]) {
  z <- eligible[eligible$pathway==p,]; np <- sum(z$NES>0); nn <- sum(z$NES<0)
  sp <- length(unique(z$same_experimental_system_id[z$NES>0])); sn <- length(unique(z$same_experimental_system_id[z$NES<0]))
  if (sp>=1 && sn>=1 && max(sp,sn)>=2) rating <- "DIVERGENCE_VALIDATED"
  else if (sp>=1 || sn>=1) rating <- "CONTEXT_DEPENDENCE_SUPPORTED"
  else rating <- "NOT_VALIDATED"
  div_rows[[p]] <- data.frame(pathway=p,n_available=nrow(z),n_positive=np,n_negative=nn,NES_min=min(z$NES),NES_max=max(z$NES),NES_range=diff(range(z$NES)),NES_SD=sd(z$NES),
    independent_systems_positive=sp,independent_systems_negative=sn,datasets_positive=length(unique(z$dataset[z$NES>0])),datasets_negative=length(unique(z$dataset[z$NES<0])),
    validation_rating=rating,role=unname(roles[p]),stringsAsFactors=FALSE)
}
divergence <- do.call(rbind,div_rows); write_tsv(divergence,file.path(r3,"results/predefined_validation/divergent_program_validation_summary.tsv"))

# Related Reactome checks, transparently mapped before reading direction.
theme_map <- data.frame(theme=c("MTORC1_related","UPR_related","IFN_related","EMT_inflammation_related"),
  regex=c("MTOR|NUTRIENT|STARVATION|TRANSLATION|RIBOSOM|PROTEIN_SYNTHESIS","UNFOLDED|ER_STRESS|PROTEIN_FOLDING","INTERFERON","EXTRACELLULAR_MATRIX|ECM_|COLLAGEN|INFLAMMAT|CYTOKINE|CHEMOKINE"),stringsAsFactors=FALSE)
write_tsv(theme_map,file.path(r3,"results/predefined_validation/reactome_theme_mapping.tsv"))
react_support <- list()
for (id in names(contrast_specs)[vapply(contrast_specs,function(s)s$primary,logical(1))]) {
  z <- gsea_cache[[paste(id,"reactome",sep="|")]]
  for (i in seq_len(nrow(theme_map))) {
    q <- z[grepl(theme_map$regex[i],z$pathway),]; if(nrow(q)){q$predefined_theme <- theme_map$theme[i];react_support[[paste(id,i)]]<-q}
  }
}
if(length(react_support)) write_tsv(do.call(rbind,react_support),file.path(r3,"results/predefined_validation/predefined_reactome_support.tsv"))

# -----------------------------------------------------------------------------
# GSE81465 temporal trajectory: sample-level standardized program score plus GSEA.
# -----------------------------------------------------------------------------
stage <- factor(group_map$GSE81465[colnames(e814)],levels=c("IgG","G1","G4","G9"))
traj <- list()
for (p in frozen) {
  genes <- intersect(pathways$hallmark[[p]],rownames(e814)); zz <- t(scale(t(e814[genes,,drop=FALSE])))
  sc <- colMeans(zz,na.rm=TRUE)
  for (st in levels(stage)) {
    v <- sc[stage==st]; tid <- if(st=="G1") "GSE81465_G1_vs_IgG" else if(st=="G4") "GSE81465_G4_vs_IgG" else if(st=="G9") "GSE81465_G9_vs_IgG" else NA_character_
    gz <- if(is.na(tid)) NULL else gsea_cache[[paste(tid,"hallmark",sep="|")]]
    rr <- if(is.null(gz)) NULL else gz[gz$pathway==p,]
    traj[[paste(p,st)]] <- data.frame(pathway=p,stage=st,generation=c(IgG=0,G1=1,G4=4,G9=9)[st],n=length(v),program_score_mean=mean(v),program_score_SE=sd(v)/sqrt(length(v)),
      vs_IgG_NES=if(is.null(rr)) 0 else rr$NES,vs_IgG_pval=if(is.null(rr)) NA else rr$pval,vs_IgG_padj=if(is.null(rr)) NA else rr$padj,stringsAsFactors=FALSE)
  }
}
traj <- do.call(rbind,traj); write_tsv(traj,file.path(r3,"results/temporal/GSE81465_predefined_program_trajectory.tsv"))

# -----------------------------------------------------------------------------
# Robustness: final-universe leave-one-dataset/contrast-out and smallest contrast.
# -----------------------------------------------------------------------------
primary <- unique(eligible[,c("pathway","dataset","contrast","evidence_tier","NES","pval","padj","same_experimental_system_id")])
robust_summary <- function(z,label,excluded="none") data.frame(analysis=label,excluded=excluded,pathway=unique(z$pathway),n_contrasts=nrow(z),n_datasets=length(unique(z$dataset)),
  n_systems=length(unique(z$same_experimental_system_id)),n_positive=sum(z$NES>0),positive_fraction=mean(z$NES>0),median_NES=median(z$NES),stringsAsFactors=FALSE)
lodo <- list()
for (p in frozen[1:2]) for (d in unique(primary$dataset)) {
  z <- primary[primary$pathway==p & primary$dataset!=d,]
  lodo[[paste(p,d)]] <- robust_summary(z,"leave_one_validation_dataset_out",d)
}
lodo <- do.call(rbind,lodo); write_tsv(lodo,file.path(r3,"results/robustness/leave_one_validation_dataset_out.tsv"))
locontrast <- list()
for (p in frozen[1:2]) for (id in unique(primary$contrast)) {
  z <- primary[primary$pathway==p & primary$contrast!=id,]
  locontrast[[paste(p,id)]] <- robust_summary(z,"leave_one_validation_contrast_out",id)
}
locontrast <- do.call(rbind,locontrast); write_tsv(locontrast,file.path(r3,"results/robustness/leave_one_validation_contrast_out.tsv"))
smallest_id <- "GSE86525_HT29_bevacizumab_resistant_vs_control"
small <- do.call(rbind,lapply(frozen[1:2],function(p) robust_summary(primary[primary$pathway==p & primary$contrast!=smallest_id,],"exclude_smallest_validation_contrast",smallest_id)))
write_tsv(small,file.path(r3,"results/robustness/exclude_smallest_validation_contrast.tsv"))

# Discovery versus validation kept separate.
disc <- read.delim(file.path(phase_root,"results/round2b/matrices/hallmark_core_NES_matrix.tsv"),check.names=FALSE)
disc_long <- do.call(rbind,lapply(seq_len(nrow(disc)),function(i)data.frame(pathway=disc$pathway[i],NES=as.numeric(disc[i,-1]),contrast=names(disc)[-1],stringsAsFactors=FALSE)))
dv <- list()
for(p in frozen){d<-disc_long[disc_long$pathway==p,];v<-primary[primary$pathway==p,];dv[[p]]<-data.frame(pathway=p,
  discovery_n=nrow(d),discovery_positive=sum(d$NES>0),discovery_positive_fraction=mean(d$NES>0),discovery_median_NES=median(d$NES),
  validation_n=nrow(v),validation_positive=sum(v$NES>0),validation_positive_fraction=mean(v$NES>0),validation_median_NES=median(v$NES),
  validation_datasets=length(unique(v$dataset)),validation_systems=length(unique(v$same_experimental_system_id)),stringsAsFactors=FALSE)}
dv <- do.call(rbind,dv); write_tsv(dv,file.path(r3,"results/predefined_validation/discovery_vs_validation_predefined_programs.tsv"))

# -----------------------------------------------------------------------------
# Primary figures.
# -----------------------------------------------------------------------------
valmat <- reshape(primary[,c("pathway","contrast","NES")],idvar="pathway",timevar="contrast",direction="wide")
rownames(valmat)<-valmat$pathway; valmat<-as.matrix(valmat[,-1]); colnames(valmat)<-sub("^NES\\.","",colnames(valmat)); valmat<-valmat[frozen,,drop=FALSE]
dmat <- as.matrix(disc[match(frozen,disc$pathway),-1]); rownames(dmat)<-frozen
heat_draw <- function(m,title,boundary=NULL){
  lim<-max(2,max(abs(m),na.rm=TRUE)); pal<-colorRampPalette(c("#2166AC","white","#B2182B"))(101)
  image(seq_len(ncol(m)),seq_len(nrow(m)),t(m[nrow(m):1,,drop=FALSE]),col=pal,zlim=c(-lim,lim),axes=FALSE,xlab="",ylab="",main=title)
  axis(1,at=seq_len(ncol(m)),labels=colnames(m),las=2,cex.axis=.58);axis(2,at=seq_len(nrow(m)),labels=rev(sub("HALLMARK_","",rownames(m))),las=2,cex.axis=.7)
  abline(h=seq(.5,nrow(m)+.5,1),v=seq(.5,ncol(m)+.5,1),col=adjustcolor("grey70",.5)); if(!is.null(boundary))abline(v=boundary+.5,lwd=4,col="black")
}
for(ext in c("png","pdf")){
  p<-file.path(r3,"results/figures",paste0("round3a_predefined_program_validation_heatmap.",ext)); if(ext=="png")png(p,1900,1200,res=170) else pdf(p,14,8)
  par(mar=c(13,17,4,2));heat_draw(valmat,"Round 3A predefined-program validation (NES)");dev.off()
  comb<-cbind(dmat,valmat);p<-file.path(r3,"results/figures",paste0("discovery_vs_validation_heatmap.",ext));if(ext=="png")png(p,2600,1300,res=170) else pdf(p,18,9)
  par(mar=c(14,17,5,2));heat_draw(comb,"Frozen programs: DISCOVERY | INDEPENDENT VALIDATION",ncol(dmat));mtext("DISCOVERY",side=3,at=(ncol(dmat)+1)/2,line=.5,font=2);mtext("VALIDATION",side=3,at=ncol(dmat)+(ncol(valmat)+1)/2,line=.5,font=2);dev.off()
}
plot_path <- function(p,out,title){
  d<-disc_long[disc_long$pathway==p,];d$set<-"Discovery";v<-primary[primary$pathway==p,c("contrast","NES","evidence_tier")];v$set<-paste("Validation",v$evidence_tier);v$evidence_tier<-NULL;z<-rbind(d[,c("contrast","NES","set")],v)
  z<-z[order(z$NES),];cols<-c(Discovery="#777777","Validation Tier B"="#1B9E77")
  draw<-function(){par(mar=c(5,18,4,2));y<-seq_len(nrow(z));plot(z$NES,y,pch=19,col=cols[z$set],yaxt="n",ylab="",xlab="NES (descriptive; not a pooled effect)",main=title,xlim=range(c(-2.7,2.7,z$NES)));axis(2,at=y,labels=z$contrast,las=2,cex.axis=.55);abline(v=0,lty=2);legend("bottomright",legend=names(cols),col=cols,pch=19,cex=.72)}
  png(paste0(out,".png"),1900,1500,res=180);draw();dev.off();pdf(paste0(out,".pdf"),13,10);draw();dev.off()
}
plot_path(frozen[1],file.path(r3,"results/figures/MTORC1_validation_plot"),"MTORC1: discovery and independent validation NES")
plot_path(frozen[2],file.path(r3,"results/figures/UPR_validation_plot"),"UPR: discovery and independent validation NES")
draw_traj<-function(){par(mfrow=c(3,3),mar=c(4,4,4,1),oma=c(0,0,1,0));for(p in frozen){z<-traj[traj$pathway==p,];yl<-range(c(z$program_score_mean-z$program_score_SE,z$program_score_mean+z$program_score_SE));plot(z$generation,z$program_score_mean,type="b",pch=19,xaxt="n",xlab="generation",ylab="mean standardized program score",main=sub("HALLMARK_","",p),cex.main=.82,cex.lab=.82,ylim=yl);arrows(z$generation,z$program_score_mean-z$program_score_SE,z$generation,z$program_score_mean+z$program_score_SE,angle=90,code=3,length=.04);axis(1,at=z$generation,labels=z$stage)}}
png(file.path(r3,"results/figures/GSE81465_predefined_program_trajectory.png"),2100,1800,res=180);draw_traj();dev.off()
pdf(file.path(r3,"results/figures/GSE81465_predefined_program_trajectory.pdf"),14,12);draw_traj();dev.off()

# Status and compact machine-readable final decision.
div_primary <- divergence[divergence$pathway %in% frozen[3:6],]
div_n <- sum(div_primary$validation_rating=="DIVERGENCE_VALIDATED")
if(mtor$validation_rating=="STRONGLY_VALIDATED" && div_n>=3) {model<-"Validation Model 1";status<-"ROUND3A_PASS_STRONG_VALIDATION"} else
if(mtor$validation_rating=="STRONGLY_VALIDATED" || mtor$validation_rating=="MODERATELY_VALIDATED") {model<-"Validation Model 2";status<-"ROUND3A_PASS_BACKBONE_VALIDATED"} else
if(div_n>=2) {model<-"Validation Model 3";status<-"ROUND3A_PASS_HETEROGENEITY_VALIDATED"} else
if(mtor$validation_rating=="WEAKLY_SUPPORTED" || any(div_primary$validation_rating=="DIVERGENCE_VALIDATED")) {model<-"Validation Model 3";status<-"ROUND3A_PASS_WITH_CAUTION"} else
{model<-"Validation Model 4";status<-"ROUND3A_HOLD"}
decision <- data.frame(primary_datasets=length(unique(primary$dataset)),primary_contrasts=length(unique(primary$contrast)),primary_systems=length(unique(primary$same_experimental_system_id)),
  MTORC1_rating=mtor$validation_rating,UPR_rating=upr$validation_rating,divergence_validated_primary_paths=div_n,validation_model=model,final_status=status,stop="ANALYSIS_COMPLETE",stringsAsFactors=FALSE)
write_tsv(decision,file.path(r3,"reports/round3a_final_decision.tsv"))
writeLines(capture.output(sessionInfo()),file.path(r3,"reports/round3a_sessionInfo.txt"))
message("Round 3A analysis complete: ",status," / ",model)
