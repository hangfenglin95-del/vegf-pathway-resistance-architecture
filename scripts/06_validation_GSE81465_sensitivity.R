#!/usr/bin/env Rscript

# Validation sensitivity analyses and interpretation reports.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
phase_root <- normalizePath(analysis_root_env, mustWork = TRUE)
r3 <- file.path(phase_root, "round3a_validation")
suppressPackageStartupMessages({library(limma); library(fgsea); library(msigdbr)})
write_tsv <- function(x,p) write.table(x,p,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
read_expr <- function(p){z<-read.delim(p,check.names=FALSE);rownames(z)<-z$gene_symbol;as.matrix(z[,-1])}
frozen <- c("HALLMARK_MTORC1_SIGNALING","HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_E2F_TARGETS")

# GSE81465 late-sample leave-one-out addresses the QC deviation of GSM2208773.
# In report-only mode, reuse the locked sensitivity output; do not rerun DE/GSEA.
reports_only <- identical(Sys.getenv("ROUND3A_REPORTS_ONLY", "0"), "1")
if (reports_only) {
  sens <- read.delim(file.path(r3,"results/robustness/GSE81465_leave_one_late_sample_out.tsv"), check.names=FALSE)
} else {
  e <- read_expr(file.path(r3,"data/processed/GSE81465_gene_expression.tsv"))
  g1 <- paste0("GSM220876",5:7); g9 <- paste0("GSM220877",1:3)
  hall <- msigdbr(db_species="HS",species="human",collection="H")
  paths <- split(hall$gene_symbol,hall$gs_name)
  sens <- list()
  for(drop in c("none",g9)){
    num <- if(drop=="none") g9 else setdiff(g9,drop); sube<-e[,c(g1,num),drop=FALSE]
    grp<-factor(c(rep("early",length(g1)),rep("late",length(num))),levels=c("early","late"))
    fit<-eBayes(lmFit(sube,model.matrix(~grp)),trend=TRUE,robust=TRUE);tt<-topTable(fit,coef="grplate",number=Inf,sort.by="none")
    ranks<-sort(setNames(tt$t,rownames(tt)),decreasing=TRUE);set.seed(20260909)
    fg<-as.data.frame(fgseaMultilevel(paths,ranks,minSize=15,maxSize=500,eps=0))
    q<-fg[match(frozen,fg$pathway),c("pathway","NES","pval","padj")];q$excluded_late_sample<-drop;q$n_late<-length(num);q$n_early<-length(g1)
    sens[[drop]]<-q
  }
  sens<-do.call(rbind,sens);sens<-sens[,c("pathway","excluded_late_sample","n_late","n_early","NES","pval","padj")]
  write_tsv(sens,file.path(r3,"results/robustness/GSE81465_leave_one_late_sample_out.tsv"))
}


capture.output(sessionInfo(), file = file.path(r3, "reports", "GSE81465_sensitivity_sessionInfo.txt"))
message("GSE81465 leave-one-late-sample-out sensitivity analysis complete")
