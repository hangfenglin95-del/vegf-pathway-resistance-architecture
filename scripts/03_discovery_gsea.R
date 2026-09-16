#!/usr/bin/env Rscript

# Per-contrast discovery GSEA and within-GSE76068 temporal classification.
# No cross-context matrices, clustering, meta-analysis, or shared/divergent calls are produced.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]), mustWork = TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
root <- normalizePath(analysis_root_env, mustWork = TRUE)
suppressPackageStartupMessages({library(fgsea); library(msigdbr); library(babelgene)})

de_dir <- file.path(root, "results", "round2a", "differential")
gsea_dir <- file.path(root, "results", "round2a", "gsea")
fig_dir <- file.path(root, "results", "round2a", "figures")
map_dir <- file.path(root, "data", "processed_round2a")
dir.create(gsea_dir, recursive=TRUE, showWarnings=FALSE); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep="\t",row.names=FALSE,quote=FALSE,na="NA")

hall_db <- msigdbr(db_species="HS", species="human", collection="H")
react_db <- msigdbr(db_species="HS", species="human", collection="C2", subcollection="CP:REACTOME")
pathways <- list(
  hallmark=split(hall_db$gene_symbol, hall_db$gs_name),
  reactome=split(react_db$gene_symbol, react_db$gs_name)
)
manifest <- data.frame(collection=c("hallmark","reactome"), msigdb_version=c(unique(hall_db$db_version),unique(react_db$db_version)),
  gene_set_count=c(length(pathways$hallmark),length(pathways$reactome)), min_size=15,max_size=500,
  engine="fgseaMultilevel", rank_statistic="limma moderated t", stringsAsFactors=FALSE)
write_tsv(manifest, file.path(root,"reports","round2a_gsea_gene_set_manifest.tsv"))

summary <- read.delim(file.path(root,"reports","round2a_differential_summary.tsv"),check.names=FALSE)
core <- c("GSE76068_tumor_Escape_vs_Response","GSE76068_stroma_Escape_vs_Response","GSE73571_tumor_Resistant_vs_Sensitive",
          "GSE180687_endothelial_Resistant_vs_Sensitive","GSE64472_stroma_Cediranib_resistant_vs_sensitive",
          "GSE64472_stroma_Vandetanib_resistant_vs_sensitive","GSE26644_tumor_Bevacizumab_resistant_vs_vehicle",
          "GSE26644_stroma_Bevacizumab_resistant_vs_vehicle")

mapping_cache <- list(); rank_summaries <- list(); gsea_summaries <- list(); result_cache <- list()
for (i in seq_len(nrow(summary))) {
  id <- summary$contrast_id[i]
  stat <- read.delim(file.path(de_dir,paste0(id,"_gene_stats.tsv")),check.names=FALSE)
  species <- unique(stat$species); unit <- sub("_(Response|Escape|Resistant|Cediranib|Vandetanib|Bevacizumab).*", "", id)
  ranks <- setNames(stat$t, stat$gene_symbol); ranks <- ranks[is.finite(ranks)]
  original_n <- length(ranks); ambiguous_dropped <- 0; duplicate_human_dropped <- 0
  if (species == "Mus musculus") {
    if (is.null(mapping_cache[[unit]])) {
      om <- orthologs(unique(names(ranks)), species="mouse", human=FALSE, min_support=3, top=FALSE)
      om <- om[!is.na(om$symbol) & !is.na(om$human_symbol) & nzchar(om$symbol) & nzchar(om$human_symbol),]
      one <- ave(om$human_symbol,om$symbol,FUN=length)==1 & ave(om$symbol,om$human_symbol,FUN=length)==1
      om$one_to_one <- one
      mapping_cache[[unit]] <- om
      write_tsv(om,file.path(map_dir,paste0(unit,"_mouse_to_human_orthologs.tsv")))
    }
    om <- mapping_cache[[unit]]; usable <- om[om$one_to_one,]
    ambiguous_dropped <- length(unique(om$symbol[!om$one_to_one]))
    ix <- match(names(ranks),usable$symbol); keep <- !is.na(ix)
    ranks <- setNames(ranks[keep],usable$human_symbol[ix[keep]])
  } else names(ranks) <- toupper(names(ranks))
  if (anyDuplicated(names(ranks))) {
    ord <- order(abs(ranks),decreasing=TRUE); ranks <- ranks[ord]; duplicate_human_dropped <- sum(duplicated(names(ranks))); ranks <- ranks[!duplicated(names(ranks))]
  }
  ranks <- sort(ranks,decreasing=TRUE)
  rank_summaries[[id]] <- data.frame(contrast_id=id,species=species,original_ranked_genes=original_n,
    one_to_one_human_ranked_genes=length(ranks),ambiguous_mouse_genes_dropped=ambiguous_dropped,
    duplicate_human_symbols_dropped=duplicate_human_dropped,stringsAsFactors=FALSE)
  for (coll in names(pathways)) {
    set.seed(20260908)
    fg <- fgseaMultilevel(pathways[[coll]],ranks,minSize=15,maxSize=500,eps=0)
    fg <- as.data.frame(fg); fg$leadingEdge <- vapply(fg$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
    fg <- fg[order(fg$padj,-abs(fg$NES)),]
    fg$contrast_id <- id; fg$contrast <- stat$contrast[1]; fg$dataset <- stat$dataset[1]; fg$compartment <- stat$compartment[1]; fg$evidence_tier <- stat$evidence_tier[1]
    fg <- fg[,c("contrast_id","contrast","dataset","compartment","evidence_tier","pathway","NES","pval","padj","log2err","ES","size","leadingEdge")]
    write_tsv(fg,file.path(gsea_dir,paste0(id,"_",coll,"_full.tsv")))
    write_tsv(fg[!is.na(fg$padj) & fg$padj<0.05,],file.path(gsea_dir,paste0(id,"_",coll,".tsv")))
    result_cache[[paste(id,coll,sep="|")]] <- fg
    gsea_summaries[[paste(id,coll)]] <- data.frame(contrast_id=id,collection=coll,tested_pathways=nrow(fg),fdr_lt_0_05=sum(fg$padj<.05,na.rm=TRUE),
      positive_fdr_lt_0_05=sum(fg$padj<.05 & fg$NES>0,na.rm=TRUE),negative_fdr_lt_0_05=sum(fg$padj<.05 & fg$NES<0,na.rm=TRUE),stringsAsFactors=FALSE)
  }
}
write_tsv(do.call(rbind,rank_summaries),file.path(root,"reports","round2a_rank_mapping_summary.tsv"))
write_tsv(do.call(rbind,rank_summaries),file.path(root,"reports","ortholog_mapping_summary.tsv"))
write_tsv(do.call(rbind,gsea_summaries),file.path(root,"reports","round2a_gsea_summary.tsv"))

# Compact top-NES figures for the eight core resistance contrasts only.
for (id in core) {
  png(file.path(fig_dir,paste0(id,"_top_NES.png")),width=1800,height=1200,res=180)
  old <- par(mfrow=c(1,2),mar=c(5,15,4,2)); on.exit(par(old),add=TRUE)
  for (coll in c("hallmark","reactome")) {
    z <- result_cache[[paste(id,coll,sep="|")]]; z <- z[order(abs(z$NES),decreasing=TRUE),][1:min(12,nrow(z)),]; z <- z[order(z$NES),]
    labs <- sub("^(HALLMARK_|REACTOME_)","",z$pathway); labs <- gsub("_"," ",labs)
    bp <- barplot(z$NES,names.arg=labs,horiz=TRUE,las=1,col=ifelse(z$NES>0,"#D95F02","#1B9E77"),border=NA,
                  xlab="NES",main=paste(id,coll),cex.names=.62)
    abline(v=0,col="grey30"); points(z$NES,bp,pch=ifelse(z$padj<.05,16,1),cex=.6)
  }
  dev.off()
  for (coll in c("hallmark","reactome")) {
    z <- result_cache[[paste(id,coll,sep="|")]]; z <- z[order(abs(z$NES),decreasing=TRUE),][1:min(16,nrow(z)),]; z <- z[order(z$NES),]
    png(file.path(fig_dir,paste0(id,"_",coll,"_topNES.png")),width=1500,height=1200,res=180)
    par(mar=c(5,17,4,2)); labs <- gsub("_"," ",sub("^(HALLMARK_|REACTOME_)","",z$pathway))
    bp <- barplot(z$NES,names.arg=labs,horiz=TRUE,las=1,col=ifelse(z$NES>0,"#D95F02","#1B9E77"),border=NA,xlab="NES",main=paste(id,coll),cex.names=.68)
    abline(v=0,col="grey30"); points(z$NES,bp,pch=ifelse(z$padj<.05,16,1),cex=.65); dev.off()
  }
}

# Qualitative temporal classes are computed only within each GSE76068 compartment.
classify <- function(a,b,c) {
  clear <- function(x) is.finite(x$NES) && (abs(x$NES)>=1.2 || x$pval<.05)
  minimal <- function(x) is.finite(x$NES) && abs(x$NES)<.5 && x$pval>=.1
  if (clear(a) && clear(b) && sign(a$NES)==sign(b$NES)) return("progressive")
  if (clear(a) && clear(b) && sign(a$NES)!=sign(b$NES)) return("rebound/recovery")
  if (minimal(a) && clear(b)) return("escape-specific")
  if (clear(a) && minimal(b)) return("treatment-effect-associated")
  "mixed/indeterminate"
}
temporal <- list()
for (comp in c("tumor","stroma")) for (coll in c("hallmark","reactome")) {
  ids <- paste0("GSE76068_",comp,"_",c("Response_vs_Pretreatment","Escape_vs_Response","Escape_vs_Pretreatment"))
  aa <- result_cache[[paste(ids[1],coll,sep="|")]]; bb <- result_cache[[paste(ids[2],coll,sep="|")]]; cc <- result_cache[[paste(ids[3],coll,sep="|")]]
  pw <- Reduce(intersect,list(aa$pathway,bb$pathway,cc$pathway))
  get <- function(z,p) z[match(p,z$pathway),c("NES","pval","padj")]
  for (p in pw) { a<-get(aa,p); b<-get(bb,p); c<-get(cc,p); temporal[[paste(comp,coll,p)]] <- data.frame(
    compartment=comp,collection=coll,pathway=p,response_vs_pretreatment_NES=a$NES,response_vs_pretreatment_padj=a$padj,
    escape_vs_response_NES=b$NES,escape_vs_response_padj=b$padj,escape_vs_pretreatment_NES=c$NES,escape_vs_pretreatment_padj=c$padj,
    temporal_class=classify(a,b,c),stringsAsFactors=FALSE) }
}
temporal <- do.call(rbind,temporal); temporal <- temporal[order(temporal$compartment,temporal$collection,temporal$temporal_class),]
write_tsv(temporal,file.path(gsea_dir,"GSE76068_within_compartment_temporal_classification.tsv"))
write_tsv(temporal,file.path(gsea_dir,"GSE76068_temporal_pathway_classification.tsv"))
write_tsv(as.data.frame(table(temporal$compartment,temporal$collection,temporal$temporal_class),stringsAsFactors=FALSE),
          file.path(root,"reports","GSE76068_temporal_classification_summary.tsv"))
writeLines(capture.output(sessionInfo()),file.path(root,"reports","round2a_gsea_sessionInfo.txt"))
message("Round 2A GSEA completed for ",nrow(summary)," contrasts; temporal classification remains within GSE76068 compartments.")
