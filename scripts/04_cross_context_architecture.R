#!/usr/bin/env Rscript

# Cross-context resistance architecture.
# Reads frozen Round 2A results. It does not modify Round 1/Round 2A files.

args_all <- commandArgs(trailingOnly=FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grep("^--file=",args_all)]),mustWork=TRUE)
analysis_root_env <- Sys.getenv("VEGF_ANALYSIS_ROOT", "")
if (!nzchar(analysis_root_env)) stop("Set VEGF_ANALYSIS_ROOT to the reconstructed analysis directory before running raw-data analyses.")
root <- normalizePath(analysis_root_env, mustWork = TRUE)
suppressPackageStartupMessages({library(fgsea); library(msigdbr)})

r2a <- file.path(root,"results","round2a"); out <- file.path(root,"results","round2b"); reports <- file.path(root,"reports","round2b")
subdirs <- c("matrices","correlations","clustering","ordination","shared_programs","divergent_programs","compartment","robustness","figures","summaries")
for (d in subdirs) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
dir.create(reports,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
write_mat <- function(m,p,row_name="pathway") { z<-data.frame(pathway=rownames(m),m,check.names=FALSE); names(z)[1]<-row_name; write_tsv(z,p) }

ids <- c(
  "GSE76068_tumor_Escape_vs_Response","GSE76068_stroma_Escape_vs_Response",
  "GSE73571_tumor_Resistant_vs_Sensitive","GSE180687_endothelial_Resistant_vs_Sensitive",
  "GSE64472_stroma_Cediranib_resistant_vs_sensitive","GSE64472_stroma_Vandetanib_resistant_vs_sensitive",
  "GSE26644_tumor_Bevacizumab_resistant_vs_vehicle","GSE26644_stroma_Bevacizumab_resistant_vs_vehicle")
short <- c("76068_T","76068_S","73571_T","180687_E","64472_Ced_S","64472_Van_S","26644_T","26644_S")
metadata <- data.frame(
  contrast_id=ids, short_id=short,
  dataset=c("GSE76068","GSE76068","GSE73571","GSE180687","GSE64472","GSE64472","GSE26644","GSE26644"),
  cancer_type=c("renal cell carcinoma","renal cell carcinoma","hepatocellular carcinoma","ovarian cancer","non-small-cell lung cancer","non-small-cell lung cancer","non-small-cell lung cancer","non-small-cell lung cancer"),
  drug=c("sunitinib","sunitinib","sorafenib","B20 anti-VEGF-A antibody","cediranib","vandetanib","bevacizumab","bevacizumab"),
  drug_class=c("multikinase VEGFR-targeting TKI","multikinase VEGFR-targeting TKI","multikinase VEGFR-targeting TKI","VEGF ligand blockade","VEGFR-targeted TKI","multikinase VEGFR-targeting TKI","VEGF ligand blockade","VEGF ligand blockade"),
  target_class=c("VEGFR/PDGFR/multikinase","VEGFR/PDGFR/multikinase","RAF/VEGFR/PDGFR multikinase","VEGF-A ligand","VEGFR1/2/3 TKI","VEGFR2/EGFR/RET TKI","VEGF-A ligand","VEGF-A ligand"),
  compartment=c("tumor","stroma","tumor","endothelial","stroma","stroma","tumor","stroma"),
  organism_source=c("human PDX tumor","mouse host stroma","human xenograft tumor","mouse tumor endothelial cells","mouse host stroma","mouse host stroma","human xenograft tumor","mouse host stroma"),
  numerator=c("escape","escape","sorafenib acquired-resistant","B20 resistant","cediranib resistant","vandetanib resistant","bevacizumab resistant","bevacizumab resistant"),
  denominator=c("response","response","sorafenib sensitive","B20 sensitive","cediranib sensitive","vandetanib sensitive","PBS vehicle","PBS vehicle"),
  evidence_tier=c("Tier A","Tier A","Tier A","Tier A","Tier A","Tier A","Tier B","Tier B"),
  biological_n_numerator=c(4,4,4,4,3,2,3,3), biological_n_denominator=c(4,4,3,4,3,3,3,3),
  readiness=c("READY","READY","READY","READY_WITH_CAUTION","READY","READY_WITH_CAUTION","READY_WITH_CAUTION","READY_WITH_CAUTION"),
  caution=c("none","none","none","GSM5468046 RMA PCA/correlation deviation; primary retained","none","resistant n=2","Tier B: resistance confounded with exposure/time/tumor size","Tier B: resistance confounded with exposure/time/tumor size"),
  same_experimental_system_id=c("GSE76068_RCC_PDX_sunitinib","GSE76068_RCC_PDX_sunitinib","GSE73571_HCC_xenograft_sorafenib","GSE180687_OVCA_xenograft_B20","GSE64472_NSCLC_xenograft_stroma","GSE64472_NSCLC_xenograft_stroma","GSE26644_NSCLC_xenograft_bevacizumab","GSE26644_NSCLC_xenograft_bevacizumab"),
  stringsAsFactors=FALSE)
write_tsv(metadata,file.path(out,"summaries","core_contrast_metadata.tsv"))

# Load full Round 2A GSEA results and construct union matrices without zero filling.
gsea <- list(hallmark=list(),reactome=list())
for (coll in names(gsea)) for (id in ids) {
  f <- file.path(r2a,"gsea",paste0(id,"_",coll,"_full.tsv")); if (!file.exists(f)) stop("Missing frozen Round 2A input: ",f)
  z <- read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
  if (anyDuplicated(z$pathway)) stop("Duplicate pathway in ",f)
  gsea[[coll]][[id]] <- z
}
make_matrix <- function(coll,value) {
  paths <- sort(unique(unlist(lapply(gsea[[coll]],function(z) z$pathway))))
  m <- matrix(NA_real_,length(paths),length(ids),dimnames=list(paths,ids))
  for (id in ids) { z<-gsea[[coll]][[id]]; m[z$pathway,id] <- z[[value]] }
  m
}
mats <- list()
for (coll in c("hallmark","reactome")) {
  mats[[coll]] <- list(NES=make_matrix(coll,"NES"),pval=make_matrix(coll,"pval"),padj=make_matrix(coll,"padj"),size=make_matrix(coll,"size"))
  write_mat(mats[[coll]]$NES,file.path(out,"matrices",paste0(coll,"_core_NES_matrix.tsv")))
  write_mat(mats[[coll]]$pval,file.path(out,"matrices",paste0(coll,"_pval_matrix.tsv")))
  write_mat(mats[[coll]]$padj,file.path(out,"matrices",paste0(coll,"_padj_matrix.tsv")))
  if (coll=="hallmark") write_mat(mats[[coll]]$size,file.path(out,"matrices","hallmark_gene_set_size_matrix.tsv"))
}

# Coverage audit: per-contrast coverage plus every pairwise overlap.
coverage <- list()
for (coll in c("hallmark","reactome")) {
  union_n <- nrow(mats[[coll]]$NES)
  for (id in ids) coverage[[paste(coll,id)]] <- data.frame(record_type="contrast",collection=coll,contrast_1=id,contrast_2=NA,
    pathway_count=sum(!is.na(mats[[coll]]$NES[,id])),union_pathway_count=union_n,overlap_count=NA,coverage_fraction=sum(!is.na(mats[[coll]]$NES[,id]))/union_n,missing_fraction=mean(is.na(mats[[coll]]$NES[,id])),stringsAsFactors=FALSE)
  cmb <- combn(ids,2,simplify=FALSE)
  for (x in cmb) coverage[[paste(coll,paste(x,collapse="|"))]] <- data.frame(record_type="pair_overlap",collection=coll,contrast_1=x[1],contrast_2=x[2],pathway_count=NA,union_pathway_count=union_n,
    overlap_count=sum(!is.na(mats[[coll]]$NES[,x[1]]) & !is.na(mats[[coll]]$NES[,x[2]])),coverage_fraction=NA,missing_fraction=NA,stringsAsFactors=FALSE)
}
coverage <- do.call(rbind,coverage); write_tsv(coverage,file.path(out,"summaries","pathway_coverage_audit.tsv"))
react_cov <- coverage[coverage$collection=="reactome" & coverage$record_type=="contrast",]
if (min(react_cov$coverage_fraction)<0.85 || max(react_cov$missing_fraction)>0.15) stop("Abnormal Reactome coverage; Round 2B stopped")

# Direction gate: labels, contrast identity, and sign(logFC)==sign(t).
direction <- list()
for (i in seq_along(ids)) {
  id<-ids[i]; d<-read.delim(file.path(r2a,"differential",paste0(id,"_gene_statistics.tsv")),check.names=FALSE)
  expected_ok <- grepl(ifelse(metadata$dataset[i]=="GSE76068","Escape_vs_Response",ifelse(metadata$dataset[i]=="GSE26644","Bevacizumab_resistant_vs_vehicle",ifelse(metadata$dataset[i]=="GSE64472",ifelse(grepl("Cediranib",id),"Cediranib_resistant_vs_sensitive","Vandetanib_resistant_vs_sensitive"),"Resistant_vs_Sensitive"))),id,fixed=TRUE)
  sign_ok <- mean(sign(d$logFC)==sign(d$t),na.rm=TRUE)==1
  direction[[id]] <- data.frame(contrast_id=id,numerator=metadata$numerator[i],denominator=metadata$denominator[i],positive_NES_interpretation=paste(metadata$numerator[i],">",metadata$denominator[i]),contrast_label_ok=expected_ok,gene_logFC_t_sign_concordance=mean(sign(d$logFC)==sign(d$t),na.rm=TRUE),status=if(expected_ok&&sign_ok) "PASS" else "FAIL")
}
direction<-do.call(rbind,direction); write_tsv(direction,file.path(out,"summaries","NES_direction_check.tsv")); if(any(direction$status!="PASS")) stop("NES direction gate failed")

# Correlation and hierarchical clustering.
cor_mats <- list(); clusters <- list()
for (coll in c("hallmark","reactome")) {
  cm <- cor(mats[[coll]]$NES,use="pairwise.complete.obs",method="spearman"); cor_mats[[coll]]<-cm
  write_mat(cm,file.path(out,"correlations",paste0(coll,"_contrast_spearman.tsv")),"contrast_id")
  for (link in c("average","complete")) {
    dm <- 1-cm; dm[dm<0] <- 0; hc <- hclust(as.dist(dm),method=link)
    clusters[[paste(coll,link,sep="|")]] <- hc
    cl <- cutree(hc,k=3)
    tab <- data.frame(contrast_id=hc$labels,leaf_order=match(hc$labels,hc$labels[hc$order]),cluster_k3=cl[hc$labels],linkage=link,distance="1-Spearman",stringsAsFactors=FALSE)
    write_tsv(tab,file.path(out,"clustering",paste0(coll,"_contrast_clustering_",link,".tsv")))
    if (coll=="hallmark" && link=="average") write_tsv(tab,file.path(out,"clustering","hallmark_contrast_clustering.tsv"))
  }
}

# PCA/ordination: raw centered primary and pathway-scaled sensitivity.
for (coll in c("hallmark","reactome")) {
  x <- mats[[coll]]$NES; keep <- rowSums(is.na(x))==0 & apply(x,1,sd)>0; x<-x[keep,]
  for (scaled in c(FALSE,TRUE)) {
    pc <- prcomp(t(x),center=TRUE,scale.=scaled); tag<-if(scaled) "scaled" else "raw_centered"
    scores <- data.frame(contrast_id=rownames(pc$x),pc$x,check.names=FALSE); loads<-data.frame(pathway=rownames(pc$rotation),pc$rotation,check.names=FALSE)
    write_tsv(scores,file.path(out,"ordination",paste0(coll,"_PCA_scores_",tag,".tsv"))); write_tsv(loads,file.path(out,"ordination",paste0(coll,"_PCA_loadings_",tag,".tsv")))
    if (coll=="hallmark" && !scaled) { write_tsv(scores,file.path(out,"ordination","hallmark_PCA_scores.tsv")); write_tsv(loads,file.path(out,"ordination","hallmark_PCA_loadings.tsv")) }
  }
}

dominant_stats <- function(nes,pval,padj,cols=colnames(nes),meta=metadata) {
  m<-nes[,cols,drop=FALSE]; pv<-pval[,cols,drop=FALSE]; pa<-padj[,cols,drop=FALSE]; mt<-meta[match(cols,meta$contrast_id),]
  outx <- lapply(seq_len(nrow(m)),function(i){v<-m[i,]; avail<-is.finite(v); pos<-sum(v[avail]>0); neg<-sum(v[avail]<0); dom<-if(pos>=neg) 1 else -1
    systems<-unique(mt$same_experimental_system_id); sys_nom<-sum(vapply(systems,function(sy){ix<-mt$same_experimental_system_id==sy & avail; any(ix & sign(v)==dom & pv[i,]<.05,na.rm=TRUE)},logical(1)))
    ta<-mt$evidence_tier=="Tier A" & avail; tap<-sum(v[ta]>0); tan<-sum(v[ta]<0)
    data.frame(pathway=rownames(m)[i],n_available=sum(avail),n_positive=pos,n_negative=neg,dominant_direction=if(dom>0) "positive" else "negative",
      direction_consistency=max(pos,neg)/sum(avail),mean_NES=mean(v[avail]),median_NES=median(v[avail]),median_abs_NES=median(abs(v[avail])),max_abs_NES=max(abs(v[avail])),
      nominal_P_lt_0_05_count=sum(pv[i,]<.05,na.rm=TRUE),FDR_lt_0_05_count=sum(pa[i,]<.05,na.rm=TRUE),independent_system_nominal_recurrence=sys_nom,
      TierA_positive=tap,TierA_negative=tan,TierA_direction_consistency=max(tap,tan)/sum(ta),stringsAsFactors=FALSE)})
  do.call(rbind,outx)
}
divergence_stats <- function(nes,pval=NULL,cols=colnames(nes),meta=metadata) {
  m<-nes[,cols,drop=FALSE]; pvmat<-if(is.null(pval)) matrix(NA_real_,nrow(m),ncol(m),dimnames=dimnames(m)) else pval[,cols,drop=FALSE]; mt<-meta[match(cols,meta$contrast_id),]
  do.call(rbind,lapply(seq_len(nrow(m)),function(i){v<-m[i,]; pv<-pvmat[i,]; ok<-is.finite(v); vv<-v[ok]; pp<-pv[ok]; mm<-mt[ok,]; pos<-vv>0; neg<-vv<0; pose<-pos & pp<.05; nege<-neg & pp<.05
    ps<-length(unique(mm$same_experimental_system_id[pos])); ns<-length(unique(mm$same_experimental_system_id[neg])); pes<-length(unique(mm$same_experimental_system_id[pose])); nesys<-length(unique(mm$same_experimental_system_id[nege])); bal<-2*min(sum(pos),sum(neg))/length(vv); sdv<-sd(vv)
    data.frame(pathway=rownames(m)[i],n_available=length(vv),NES_range=max(vv)-min(vv),NES_SD=sdv,NES_IQR=IQR(vv),n_positive=sum(pos),n_negative=sum(neg),sign_balance_factor=bal,
      max_positive_NES=if(any(pos))max(vv[pos]) else NA,min_negative_NES=if(any(neg))min(vv[neg]) else NA,positive_system_count=ps,negative_system_count=ns,
      positive_nominal_count=sum(pose,na.rm=TRUE),negative_nominal_count=sum(nege,na.rm=TRUE),positive_nominal_system_count=pes,negative_nominal_system_count=nesys,
      divergence_score=sdv*bal,not_single_extreme=(sum(abs(vv)>=.75)>=3),stringsAsFactors=FALSE)}))
}

hall_shared <- dominant_stats(mats$hallmark$NES,mats$hallmark$pval,mats$hallmark$padj)
hall_shared$strong_candidate <- with(hall_shared,n_available==8 & pmax(n_positive,n_negative)>=6 & pmax(TierA_positive,TierA_negative)>=5 & median_abs_NES>=1 & independent_system_nominal_recurrence>=2)
hall_shared$moderate_candidate <- with(hall_shared,n_available==8 & pmax(n_positive,n_negative)>=6 & pmax(TierA_positive,TierA_negative)>=5 & median_abs_NES>=.75)
hall_shared$weak_directional_candidate <- with(hall_shared,n_available==8 & pmax(n_positive,n_negative)>=6 & median_abs_NES>=.5)
hall_shared$shared_tier <- ifelse(hall_shared$strong_candidate,"strong",ifelse(hall_shared$moderate_candidate,"moderate",ifelse(hall_shared$weak_directional_candidate,"weak","not_candidate")))
hall_shared<-hall_shared[order(match(hall_shared$shared_tier,c("strong","moderate","weak","not_candidate")),-hall_shared$direction_consistency,-hall_shared$median_abs_NES),]
write_tsv(hall_shared,file.path(out,"shared_programs","hallmark_shared_program_statistics.tsv")); write_tsv(hall_shared[hall_shared$shared_tier!="not_candidate",],file.path(out,"shared_programs","hallmark_shared_program_candidates.tsv"))
thresholds <- data.frame(tier=c("strong","moderate","weak"),definition=c(">=6/8 same direction; >=5/6 Tier A; median|NES|>=1; nominal recurrence in >=2 independent systems",">=6/8 same direction; >=5/6 Tier A; median|NES|>=0.75",">=6/8 same direction; median|NES|>=0.5"),candidate_count=c(sum(hall_shared$strong_candidate),sum(hall_shared$moderate_candidate),sum(hall_shared$weak_directional_candidate)),stringsAsFactors=FALSE)
write_tsv(thresholds,file.path(out,"shared_programs","hallmark_shared_threshold_sensitivity.tsv"))

hall_div <- divergence_stats(mats$hallmark$NES,mats$hallmark$pval)
hall_div$strong_preliminary <- with(hall_div,n_positive>=3 & n_negative>=3 & positive_nominal_count>=2 & negative_nominal_count>=2 & positive_nominal_system_count>=2 & negative_nominal_system_count>=2 & NES_range>=3 & not_single_extreme)
hall_div<-hall_div[order(-hall_div$strong_preliminary,-hall_div$divergence_score,-hall_div$NES_range),]

# Compartment-focused descriptive and within-study paired comparisons.
tum<-metadata$contrast_id[metadata$compartment=="tumor"]; str<-metadata$contrast_id[metadata$compartment=="stroma"]
comp <- data.frame(pathway=rownames(mats$hallmark$NES),tumor_median_NES=apply(mats$hallmark$NES[,tum,drop=FALSE],1,median,na.rm=TRUE),stroma_median_NES=apply(mats$hallmark$NES[,str,drop=FALSE],1,median,na.rm=TRUE),stringsAsFactors=FALSE)
comp$delta_tumor_minus_stroma<-comp$tumor_median_NES-comp$stroma_median_NES; comp$analysis_label<-"exploratory descriptive comparison; contexts non-independent and confounded"
comp<-comp[order(-abs(comp$delta_tumor_minus_stroma)),]; write_tsv(comp,file.path(out,"compartment","hallmark_tumor_vs_stroma_summary.tsv"))
pd <- data.frame(pathway=rownames(mats$hallmark$NES),
  delta_GSE76068=mats$hallmark$NES[,ids[1]]-mats$hallmark$NES[,ids[2]],
  delta_GSE26644=mats$hallmark$NES[,ids[7]]-mats$hallmark$NES[,ids[8]],stringsAsFactors=FALSE)
pd$direction_agreement<-ifelse(sign(pd$delta_GSE76068)==sign(pd$delta_GSE26644),ifelse(pd$delta_GSE76068>0,"tumor_higher_both","stroma_higher_both"),"discordant")
pd$mean_abs_delta<-(abs(pd$delta_GSE76068)+abs(pd$delta_GSE26644))/2; pd<-pd[order(-pd$mean_abs_delta),]
write_tsv(pd,file.path(out,"compartment","paired_tumor_stroma_deltas.tsv"))
paired_cor <- cor(pd$delta_GSE76068,pd$delta_GSE26644,method="spearman",use="complete.obs")
write_tsv(data.frame(collection="hallmark",comparison="GSE76068 delta vs GSE26644 delta",spearman=paired_cor,n_pathways=nrow(pd)),file.path(out,"compartment","paired_tumor_stroma_delta_correlation.tsv"))

# GSE64472 within-study Cediranib vs Vandetanib, Hallmark + Reactome.
cv <- list()
for (coll in c("hallmark","reactome")) { m<-mats[[coll]]$NES; z<-data.frame(collection=coll,pathway=rownames(m),Cediranib_NES=m[,ids[5]],Vandetanib_NES=m[,ids[6]],stringsAsFactors=FALSE); z$delta_Ced_minus_Van<-z$Cediranib_NES-z$Vandetanib_NES; z$direction_concordant<-sign(z$Cediranib_NES)==sign(z$Vandetanib_NES); z$analysis_label<-"exploratory / caution: Vandetanib resistant n=2"; cv[[coll]]<-z }
cv<-do.call(rbind,cv); cv<-cv[order(cv$collection,-abs(cv$delta_Ced_minus_Van)),]; write_tsv(cv,file.path(out,"compartment","GSE64472_Ced_vs_Van_pathway_comparison.tsv"))
cv_cor<-do.call(rbind,lapply(split(cv,cv$collection),function(z)data.frame(collection=z$collection[1],spearman=cor(z$Cediranib_NES,z$Vandetanib_NES,use="complete.obs",method="spearman"),n_pathways=sum(complete.cases(z[,c("Cediranib_NES","Vandetanib_NES")])))))
write_tsv(cv_cor,file.path(out,"compartment","GSE64472_Ced_vs_Van_correlation.tsv"))

# Tier-A-only architecture.
tierA <- metadata$contrast_id[metadata$evidence_tier=="Tier A"]
ta_shared<-dominant_stats(mats$hallmark$NES,mats$hallmark$pval,mats$hallmark$padj,tierA)
ta_div<-divergence_stats(mats$hallmark$NES,mats$hallmark$pval,tierA)
ta<-merge(ta_shared,ta_div,by="pathway",suffixes=c("_shared","_divergent")); ta<-ta[order(-ta$direction_consistency,-ta$median_abs_NES),]
write_tsv(ta,file.path(out,"robustness","tierA_only_hallmark_summary.tsv"))
ta_cor<-cor(mats$hallmark$NES[,tierA],method="spearman",use="pairwise.complete.obs"); write_mat(ta_cor,file.path(out,"robustness","tierA_only_hallmark_correlation.tsv"),"contrast_id")
ta_dm<-1-ta_cor; ta_dm[ta_dm<0]<-0; ta_hc<-hclust(as.dist(ta_dm),method="average"); write_tsv(data.frame(contrast_id=ta_hc$labels,leaf_order=match(ta_hc$labels,ta_hc$labels[ta_hc$order]),cluster_k3=cutree(ta_hc,3)[ta_hc$labels]),file.path(out,"robustness","tierA_only_hallmark_clustering.tsv"))

# Leave-one-contrast-out and leave-one-dataset-out recalculation for every Hallmark pathway.
loco<-list(); for (omit in ids) { keep<-setdiff(ids,omit); a<-dominant_stats(mats$hallmark$NES,mats$hallmark$pval,mats$hallmark$padj,keep); b<-divergence_stats(mats$hallmark$NES,mats$hallmark$pval,keep); z<-merge(a,b,by="pathway",suffixes=c("_shared","_divergent")); z$omitted_contrast<-omit; loco[[omit]]<-z }
loco<-do.call(rbind,loco); write_tsv(loco,file.path(out,"robustness","leave_one_contrast_out.tsv"))
lodo<-list(); for (omit in unique(metadata$dataset)) { keep<-metadata$contrast_id[metadata$dataset!=omit]; a<-dominant_stats(mats$hallmark$NES,mats$hallmark$pval,mats$hallmark$padj,keep); b<-divergence_stats(mats$hallmark$NES,mats$hallmark$pval,keep); z<-merge(a,b,by="pathway",suffixes=c("_shared","_divergent")); z$omitted_dataset<-omit; lodo[[omit]]<-z }
lodo<-do.call(rbind,lodo); write_tsv(lodo,file.path(out,"robustness","leave_one_dataset_out.tsv"))

# Attach leave-one-dataset-out robustness to divergent and shared candidates.
lodo_split <- tapply(seq_len(nrow(lodo)),lodo$pathway,function(ix) all(lodo$n_positive_divergent[ix]>=1 & lodo$n_negative_divergent[ix]>=1))
lodo_shared_cons <- tapply(lodo$direction_consistency,lodo$pathway,min,na.rm=TRUE); lodo_shared_mag<-tapply(lodo$median_abs_NES,lodo$pathway,min,na.rm=TRUE)
hall_div$lodo_split_all <- unname(lodo_split[hall_div$pathway]); hall_div$strong_candidate <- hall_div$strong_preliminary & hall_div$lodo_split_all
write_tsv(hall_div,file.path(out,"divergent_programs","hallmark_divergent_program_statistics.tsv")); write_tsv(hall_div[hall_div$strong_candidate,],file.path(out,"divergent_programs","hallmark_divergent_program_candidates.tsv"))
hall_shared$lodo_min_direction_consistency<-unname(lodo_shared_cons[hall_shared$pathway]); hall_shared$lodo_min_median_abs_NES<-unname(lodo_shared_mag[hall_shared$pathway]); hall_shared$lodo_robust<-hall_shared$lodo_min_direction_consistency>=2/3 & hall_shared$lodo_min_median_abs_NES>=.5
write_tsv(hall_shared,file.path(out,"shared_programs","hallmark_shared_program_statistics.tsv")); write_tsv(hall_shared[hall_shared$shared_tier!="not_candidate",],file.path(out,"shared_programs","hallmark_shared_program_candidates.tsv"))

# Exclude Vandetanib architecture sensitivity.
keepv<-setdiff(ids,ids[6]); evs<-dominant_stats(mats$hallmark$NES,mats$hallmark$pval,mats$hallmark$padj,keepv); evd<-divergence_stats(mats$hallmark$NES,mats$hallmark$pval,keepv); ev<-merge(evs,evd,by="pathway",suffixes=c("_shared","_divergent")); ev$primary_shared_tier<-hall_shared$shared_tier[match(ev$pathway,hall_shared$pathway)]; ev$primary_strong_divergent<-hall_div$strong_candidate[match(ev$pathway,hall_div$pathway)]
write_tsv(ev,file.path(out,"robustness","exclude_vandetanib_summary.tsv")); evc<-cor(mats$hallmark$NES[,keepv],method="spearman"); evdm<-1-evc; evdm[evdm<0]<-0; evh<-hclust(as.dist(evdm),method="average"); write_tsv(data.frame(contrast_id=evh$labels,leaf_order=match(evh$labels,evh$labels[evh$order]),cluster_k3=cutree(evh,3)[evh$labels]),file.path(out,"robustness","exclude_vandetanib_clustering.tsv"))

# GSE180687 targeted pathway-level leave-one-out sensitivity (4 resistant vs 3 sensitive).
loo <- read.delim(file.path(r2a,"qc","GSE180687_leave_GSM5468046_out_gene_stats.tsv"),check.names=FALSE)
omap <- read.delim(file.path(root,"data","processed_round2a","GSE180687_endothelial_mouse_to_human_orthologs.tsv"),check.names=FALSE)
use <- omap[omap$one_to_one,]; ix<-match(loo$gene_symbol,use$symbol); ok<-!is.na(ix)&is.finite(loo$t); ranks<-setNames(loo$t[ok],use$human_symbol[ix[ok]]); ranks<-sort(ranks[!duplicated(names(ranks))],decreasing=TRUE)
hall_db<-msigdbr(db_species="HS",species="human",collection="H"); react_db<-msigdbr(db_species="HS",species="human",collection="C2",subcollection="CP:REACTOME")
sets<-list(hallmark=split(hall_db$gene_symbol,hall_db$gs_name),reactome=split(react_db$gene_symbol,react_db$gs_name)); stabs<-list()
for (coll in names(sets)) { set.seed(20260908); fg<-as.data.frame(fgseaMultilevel(sets[[coll]],ranks,minSize=15,maxSize=500,eps=0)); primary<-gsea[[coll]][[ids[4]]]; z<-merge(primary[,c("pathway","NES","pval","padj")],fg[,c("pathway","NES","pval","padj")],by="pathway",suffixes=c("_primary","_leave_GSM5468046_out")); z$collection<-coll; z$NES_delta<-z$NES_leave_GSM5468046_out-z$NES_primary; z$direction_concordant<-sign(z$NES_primary)==sign(z$NES_leave_GSM5468046_out); z$spearman_all<-cor(z$NES_primary,z$NES_leave_GSM5468046_out,method="spearman"); stabs[[coll]]<-z }
stabs<-do.call(rbind,stabs); write_tsv(stabs,file.path(out,"robustness","GSE180687_leave_one_out_pathway_stability.tsv"))

# Reactome summary and deterministic name-family redundancy grouping.
react_shared<-dominant_stats(mats$reactome$NES,mats$reactome$pval,mats$reactome$padj); react_div<-divergence_stats(mats$reactome$NES,mats$reactome$pval)
family <- function(x) { y<-sub("^REACTOME_","",x); rules<-list(INTERFERON="INTERFERON|ANTIVIRAL",IMMUNE_CYTOKINE="IMMUNE|CYTOKINE|CHEMOKINE|INTERLEUKIN|TNF|COMPLEMENT|INFLAM",CELL_CYCLE="MITOTIC|CELL_CYCLE|E2F|G2_M|DNA_REPLICATION",RTK_MAPK="EGFR|ERBB|FGFR|MET_|MAPK|RAF|RAS",VEGF_ANGIOGENESIS="VEGF|VEGFR|ANGIO|VASCULAR",ECM_ADHESION="EXTRACELLULAR_MATRIX|COLLAGEN|INTEGRIN|CELL_JUNCTION|ADHESION",LIPID_METABOLISM="LIPID|FATTY_ACID|CHOLESTEROL|LIPOPROTEIN",RNA_TRANSLATION="RIBOSOM|TRANSLATION|RNA_|MRNA|RRNA",DNA_REPAIR="DNA_REPAIR|DOUBLE_STRAND|NUCLEOTIDE_EXCISION|HOMOLOGOUS_RECOMBINATION",METABOLISM="METABOLISM|GLYCOLYSIS|TCA|RESPIRATORY|OXIDATIVE")
    for(n in names(rules)) if(grepl(rules[[n]],y)) return(n); paste(strsplit(y,"_")[[1]][1:min(3,length(strsplit(y,"_")[[1]]))],collapse="_") }
make_themes <- function(stat,type) { if(type=="shared") keep<-stat$n_available>=6 & stat$direction_consistency>=.75 & stat$median_abs_NES>=.75 & stat$independent_system_nominal_recurrence>=2 else keep<-stat$n_available>=6 & stat$n_positive>=2 & stat$n_negative>=2 & stat$positive_nominal_count>=2 & stat$negative_nominal_count>=2 & stat$NES_range>=3
  stat<-stat[keep,]; score<-if(type=="shared") stat$direction_consistency*stat$median_abs_NES else stat$divergence_score; ord<-order(score,decreasing=TRUE); top<-stat[ord[1:min(150,nrow(stat))],]; top$theme<-vapply(top$pathway,family,character(1)); sp<-split(top,top$theme); do.call(rbind,lapply(sp,function(z){sc<-if(type=="shared")z$direction_consistency*z$median_abs_NES else z$divergence_score; rep<-z$pathway[which.max(sc)]; pattern<-paste(sprintf("%s=%+.2f",metadata$short_id,mats$reactome$NES[rep,metadata$contrast_id]),collapse=";"); data.frame(theme=z$theme[1],representative_pathway=rep,member_pathways=paste(z$pathway,collapse=";"),member_count=nrow(z),NES_pattern=pattern,contributing_contexts=paste(metadata$short_id[abs(mats$reactome$NES[rep,])>=1 & is.finite(mats$reactome$NES[rep,])],collapse=";"),theme_score=max(sc),stringsAsFactors=FALSE)})) }
write_tsv(make_themes(react_shared,"shared"),file.path(out,"shared_programs","reactome_shared_themes.tsv")); write_tsv(make_themes(react_div,"divergent"),file.path(out,"divergent_programs","reactome_divergent_themes.tsv"))

# Record probe summarization robustness without rerunning every pathway analysis.
ps<-read.delim(file.path(root,"reports","round2a_probe_collapse_sensitivity.tsv"),check.names=FALSE)
writeLines(c("# Probe summarization robustness","","Round 2B primary analysis reads the frozen highest-variance-probe Round 2A GSEA outputs.","",
  sprintf("Across the 12 Round 2A contrasts, moderated-t ranking Spearman ranged from %.3f to %.3f; top-100 direction concordance ranged from %.0f%% to %.0f%%.",min(ps$moderated_t_spearman),max(ps$moderated_t_spearman),100*min(ps$top100_sign_concordance),100*max(ps$top100_sign_concordance)),"",
  "No global pathway rerun was triggered. Targeted pathway sensitivity is reserved for a core candidate that depends on a low-stability contrast."),file.path(reports,"probe_summarization_robustness.md"))

# Compact numerical summaries for reporting.
off <- function(x) x[upper.tri(x)]
arch <- data.frame(metric=c("hallmark_pairwise_mean","hallmark_pairwise_median","hallmark_pairwise_min","hallmark_pairwise_max","reactome_pairwise_mean","reactome_pairwise_median","reactome_pairwise_min","reactome_pairwise_max","paired_tumor_stroma_delta_spearman","GSE64472_Ced_Van_hallmark_spearman","GSE180687_LOO_hallmark_spearman"),
 value=c(mean(off(cor_mats$hallmark)),median(off(cor_mats$hallmark)),min(off(cor_mats$hallmark)),max(off(cor_mats$hallmark)),mean(off(cor_mats$reactome)),median(off(cor_mats$reactome)),min(off(cor_mats$reactome)),max(off(cor_mats$reactome)),paired_cor,cv_cor$spearman[cv_cor$collection=="hallmark"],unique(stabs$spearman_all[stabs$collection=="hallmark"])))
write_tsv(arch,file.path(out,"summaries","architecture_numeric_summary.tsv"))
writeLines(capture.output(sessionInfo()),file.path(reports,"round2b_sessionInfo.txt"))
message("Round 2B architecture tables completed")
