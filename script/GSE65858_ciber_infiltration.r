library('dplyr')
library('survival')
library('survminer')
library("GSVA")
library('GSEABase')
library('estimate')
library("ComplexHeatmap")
library('circlize')
library('ConsensusClusterPlus')
library('ggthemes')

setwd("working_dir")
GSE<-read.table("input/GSE65858_exp.matrix.txt",
                sep = "\t",header = T,stringsAsFactors = F,check.names = F)

hnsc_tpm<-read.table("input/HNSC_RSEM_genes/tcga_hnsc_tpm_os.matrix.txt",
                     sep = "\t",header = T,stringsAsFactors = F,check.names = F)
hnsc_tpm<-hnsc_tpm[,-20535] ##delect last empty column

cc<-colnames(hnsc_tpm)[which(colnames(hnsc_tpm) %in% colnames(GSE))]
GSE<-GSE[,as.character(cc)]
GSE$OS_time<-(as.numeric(GSE$OS_time)/365)
rownames(GSE)<-GSE$SampleID
GSE<-as.data.frame(GSE)
GSE_for_ciber_pre<-GSE[5:ncol(GSE)]
GSE_for_ciber<-as.data.frame(t(GSE_for_ciber_pre))
colnames(GSE_for_ciber)<-as.character(colnames(GSE_for_ciber))
##
  ##
write.table(GSE_for_ciber,"input/cibersort/GSE_for_cibersort.txt",sep = "\t",quote = F)
######################
## you could run CIBERSORT at https://cibersort.stanford.edu/runcibersort.php (Permutations = 1000) using "GSE_for_cibersort.txt.txt" as input to get CIBERSORT output file 
fiber_res<-read.table("OUT/GSE_CIBERSORT_result.txt",
                      sep = "\t",header = T,stringsAsFactors = F,check.names = F)

rownames(fiber_res)<-fiber_res$`Input Sample`
GSE_ciber <- fiber_res[,2:23]
GSE_ciber_mat<-as.matrix(GSE_ciber)
GSE_ciber_mat[1:4,1:4]

samlist<-rownames(GSE_ciber_mat)
samlist<-as.data.frame(samlist)
colnames(samlist)="ann_samp"
samlist$index<-1:nrow(samlist) ###remember the order

###get the annotation, set complexheatmap annotation color parameters###
col_stromal = colorRamp2(c(-2000, 0, 2000), c("cadetblue1", "deepskyblue", "royalblue4"))
col_stromal(seq(-2000, 2000))
col_imm = colorRamp2(c(-1000, 0,2000,4000), c("ivory1", "burlywood2","coral", "coral4"))
col_imm(seq(-1000, 4000))
col_est = colorRamp2(c(-4000, 0,4000,6000), c("darkolivegreen1", "darkolivegreen3","forestgreen", "darkgreen"))
col_est(seq(-1000, 4000))

scores_gse=read.table("OUT/GSE_estimate_score.gct",skip = 2,header = T,check.names = F)
rownames(scores_gse)=scores_gse[,1]
scores_gse=t(scores_gse[,3:ncol(scores_gse)])

scores_est<-scores_gse
scores_est<-as.matrix(scores_est)
scores_est<-as.data.frame(scores_est)
scores_est$sampleID <- rownames(scores_est)

###15-gene signature score###
gene_score_15<-read.table("OUT/GSE65858_15gene_signature_scores_groupsBy3cutoffs.txt",
                          sep = "\t",header = T,stringsAsFactors = F,check.names = F)

##use the cutoff of Maxstat for 15-gene signature score##
sig_scores<-cbind(rownames(gene_score_15),gene_score_15$groupMS)
sig_scores<-as.data.frame(sig_scores)
colnames(sig_scores)<-c("sampleID","group")
####

nmf2_cluster<-read.table("OUT/gse65858_top500genesfromTCGA_nmf2_cluster_samples.txt",
                         sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(nmf2_cluster)<-c("sampleID","class","NMF_cluster","OS_time","OS_stage")

#####heatmap for hierachical clustering####
scaled_mat = t(scale(GSE_ciber_mat))

rcc6 = ConsensusClusterPlus(scaled_mat,maxK=6,reps=50,
                            pItem=0.8,pFeature=1,title="GSE65858",
                            distance='euclidean',
                            clusterAlg="pam",plot="png")

saveRDS(rcc6, file = "OUT/GSE65858_CIBERSORT_ConsensClusterPlus_maxK_6.rds")
rcc6<-readRDS(file="OUT/GSE65858_CIBERSORT_ConsensClusterPlus_maxK_6.rds")
consensusClass <- rcc6[[2]][["consensusClass"]]
consensusClass<-as.data.frame(consensusClass)
consensusClass$sampleID<-NA
consensusClass$sampleID<-rownames(consensusClass)

ann0<-merge(samlist,scores_est,by.x="ann_samp",by.y="sampleID",all.x=T)
ann1<-merge(ann0,sig_scores,by.x="ann_samp",by.y="sampleID",all.x=T)
ann2<-merge(ann1,nmf2_cluster,by.x="ann_samp",by.y="sampleID",all.x=T)
ann3<-merge(ann2,consensusClass,by.x="ann_samp",by.y="sampleID",all.x=T)
##############
ann3$ICI<-NA
ann3$ICI[which(ann3$consensusClass==1)]<-"ICI_A"
ann3$ICI[which(ann3$consensusClass==2)]<-"ICI_B"
ann3$ICI<-factor(ann3$ICI,levels = c("ICI_A", "ICI_B"))
ann3$NMF_cluster<-factor(ann3$NMF_cluster,levels = c("AS1", "AS2"))
ann3$group<-factor(ann3$group,levels = c("low", "high"))

pdf("OUT/GSE_cibersort_consensusClass_2clusters_survival.pdf",width=6,height=6)
fit_immune_os<-survfit(Surv(OS_time,OS_stage)~ICI,data = ann3)
ggsurvplot(fit_immune_os, data = ann3, pval = T,
           legend = "top",surv.median.line="hv",risk.table = T,
           legend.title = "",
           title = "Immune_subtype",
           font.legend = c(18, "plain"),
           risk.table.fontsize = 5.5,
           risk.table.y.text=FALSE,
           xlab = "Years"
           )
dev.off()

#####plot estimate score, immune score####
my_comparisons <- list(c("ICI_A","ICI_B"))
p1<-ggplot(ann3, aes(x=ICI, y=StromalScore,color=ICI)) + 
  geom_violin(trim=FALSE,show.legend = F)+ geom_boxplot(width=0.1,show.legend = F)+
  stat_compare_means(aes(group = ICI),comparisons = my_comparisons)+
  ylab("StromalScore")+labs(title = "StromalScore")+theme(plot.title = element_text(hjust = 0.5))

p2<-ggplot(ann3, aes(x=ICI, y=ImmuneScore,color=ICI)) + 
  geom_violin(trim=FALSE,show.legend = F)+ geom_boxplot(width=0.1,show.legend = F)+
  stat_compare_means(aes(group = ICI),comparisons = my_comparisons)+
  ylab("ImmuneScore")+labs(title = "ImmuneScore")+theme(plot.title = element_text(hjust = 0.5))

p3<-ggplot(ann3, aes(x=ICI, y=ESTIMATEScore,color=ICI)) + 
  geom_violin(trim=FALSE,show.legend = F)+ geom_boxplot(width=0.1,show.legend = F)+
  stat_compare_means(aes(group = ICI),comparisons = my_comparisons)+
  ylab("ESTIMATEScore")+labs(title = "ESTIMATEScore")+theme(plot.title = element_text(hjust = 0.5))

#########
pdf(file="OUT/GSE_cibersort_consensusClass_boxplot.pdf",width=6,height=3)
ggarrange(p1, p2, p3, 
          ncol = 3, nrow = 1)
dev.off()

##22 cell type plots##
library('tidyr')
library('tidyverse')
library(reshape2)
plot_info_raw<-GSE_ciber_mat
plot_info_raw<-as.data.frame(plot_info_raw)
plot_info_raw$sampleID<-NA
plot_info_raw$sampleID<-rownames(plot_info_raw)
ann4<-as.data.frame(cbind(as.character(ann3$ann_samp),as.character(ann3$ICI),ann3$index))
colnames(ann4)<-c("ann_samp","ICI","index")
plot_info_raw<-merge(ann4,plot_info_raw,by.x="ann_samp",by.y="sampleID",all.x=T)

plot_info<-as.data.frame(plot_info_raw)
plot_info<-plot_info[,-c(1,3)]
plot_info$ICI<-as.character(plot_info$ICI)

md <- melt(plot_info,id=c("ICI"),measure=-c(1))
md$ICI<-as.factor(md$ICI)

pdf(file="OUT/GSE_cibersort_consensusClass_cells_boxplot.pdf",width=7,height=5.5)
ggboxplot(md,x = "variable",y = "value",color = "black",
          fill = "ICI",
          ylab = "fraction",
          xlab = "cells",
          main = "TME Cell composition (GSE65858)")+
  stat_compare_means(label = "p.signif",
                      method = "t.test",
                      ref.group = ".all.",
                      hide.ns = T)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5,size = 13))

dev.off()

####

multi_var<-read.table("OUT/gse65858_signature_scores_groupsBy3cutoffs_multivariate_cli.txt",
                      sep = "\t",header = T,stringsAsFactors = F,check.names = F)
multi_var$Age[multi_var$Age=="<=60"]<-as.character("≤60")
multi_var<-merge(ann3,multi_var,by.x="ann_samp",by.y="SampleID",all.x=T)
multi_var$Age<-factor(multi_var$Age,levels = c('≤60',">60"))
multi_var$HPV[is.na(multi_var$HPV)]<-"N/A"
multi_var$HPV<-factor(multi_var$HPV,levels = c("HPV16 DNA+RNA+","HPV16 DNA+RNA-","Other HPVs","Negative","N/A"))
multi_var$Primary_site[is.na(multi_var$Primary_site)]<-"N/A"
multi_var$Primary_site<-factor(multi_var$Primary_site,levels = c("OPSCC","Non-OPSCC","N/A"))
multi_var$Uicc_stage<-factor(multi_var$Uicc_stage,levels = c("I-II","III","IV"))
multi_var$Smoking<-factor(multi_var$Smoking,levels = c("No","Yes"))

write.table(multi_var,"OUT/GSE65858_hnsc_signature_scores_groupsBy3cutoffs_multivariate_cli_NMFs_ICIs.txt",sep = "\t",quote = F)

pdf("OUT/GSE65858_cibersort_Heatmap_split_by_riskscore.pdf",width=11.5,height=7)
Heatmap(scaled_mat,column_km_repeats = 500,name = " ",
             top_annotation = HeatmapAnnotation(
               simple_anno_size = unit(0.35, "cm"),
               HPV=multi_var$HPV,
               Primary_site=multi_var$Primary_site,
               Uicc_stage=multi_var$Uicc_stage,
               Smoking=multi_var$Smoking,
               StromalScore=multi_var$StromalScore,
               ImmuneScore=multi_var$ImmuneScore,
               ESTIMATEScore=multi_var$ESTIMATEScore,
               NMF_cluster=multi_var$NMF_cluster,
               Risk_score=multi_var$group,
               ICI=multi_var$ICI,
               col=list(
                 HPV=c("HPV16 DNA+RNA+"="coral3","HPV16 DNA+RNA-"="yellowgreen","Other HPVs"="lightpink","Negative"="cadetblue1","N/A"="black"),
                 Primary_site=c("OPSCC"="coral3","Non-OPSCC"="yellowgreen","N/A"="black"),
                 Uicc_stage=c("I-II"="cadetblue1","III"="plum1","IV"="brown1"),
                 Smoking=c("No"="cadetblue1","Yes"="#E7B800"),
                 StromalScore=col_stromal,
                 ImmuneScore=col_imm,
                 ESTIMATEScore=col_est,
                 NMF_cluster=c("AS1"="cadetblue1","AS2"="#E7B800"),
                 Risk_score=c("low"="cadetblue1","high"="#E7B800"),
                 ICI=c("ICI_A"="cadetblue1","ICI_B"="#E7B800")
               )
             ),
             column_split = multi_var$group,
             show_column_names=F,
)
dev.off()
#
######
#compare two different ICI clusters in high risk signature score patients
pdf("OUT/GSE65858_signature_ici_nmf_survival_mix.pdf",width=5.5,height=3)
source("script/ggsurvplot_facet_risktable.R")
### ICI and signature
fit_i_h<-survfit(Surv(OS_time.x,OS_stage.x)~ICI,data = multi_var)
ggsurvplot_facet(fit_i_h,facet.by = c("Signature_score"),data = multi_var,pval=T,xlab = "Years")
fit_i_h_table <- survfit(Surv(OS_time.x,OS_stage.x) ~ ICI + Signature_score, multi_var)
ggsurvplot_facet_risktable(fit_i_h_table, multi_var,pval=T)

### signature and ICI
fit_h_i<-survfit(Surv(OS_time.x,OS_stage.x)~Signature_score,data = multi_var)
ggsurvplot_facet(fit_j_i,facet.by = c("ICI"),data = multi_var,pval=T,xlab = "Years")
fit_h_i_table <- survfit(Surv(OS_time.x,OS_stage.x) ~ Signature_score + ICI, multi_var)
ggsurvplot_facet_risktable(fit_j_i_table, multi_var,pval=T)

### ICI and NMF
fit_i_n<-survfit(Surv(OS_time.x,OS_stage.x)~ICI,data = multi_var)
ggsurvplot_facet(fit_i_n,facet.by = c("NMF_cluster"),data = multi_var,pval=T,xlab = "Years")
fit_i_n_table <- survfit(Surv(OS_time.x,OS_stage.x) ~ ICI + NMF_cluster, multi_var)
ggsurvplot_facet_risktable(fit_i_n_table, multi_var,pval=T,xlab = "Years")

### NMF and ICI
fit_n_i<-survfit(Surv(OS_time.x,OS_stage.x)~NMF_cluster,data = multi_var)
ggsurvplot_facet(fit_n_i,facet.by = c("ICI"),data = multi_var,pval=T,xlab = "Years")
fit_n_i_table <- survfit(Surv(OS_time.x,OS_stage.x) ~ ICI + NMF_cluster, multi_var)
ggsurvplot_facet_risktable(fit_n_i_table, multi_var,pval=T,xlab = "Years")

### Signature and NMF
fit_h_n<-survfit(Surv(OS_time.x,OS_stage.x)~Signature_score,data = multi_var)
ggsurvplot_facet(fit_h_n,facet.by = c("NMF_cluster"),data = multi_var,pval=T,xlab = "Years")
fit_h_n_table <- survfit(Surv(OS_time.x,OS_stage.x) ~ Signature_score + NMF_cluster, multi_var)
ggsurvplot_facet_risktable(fit_h_n_table, multi_var,pval=T,xlab = "Years")

### NMF and Signature 
fit_n_h<-survfit(Surv(OS_time.x,OS_stage.x)~NMF_cluster,data = multi_var)
ggsurvplot_facet(fit_n_h,facet.by = c("Signature_score"),data = multi_var,pval=T,xlab = "Years")
fit_n_h_table <- survfit(Surv(OS_time.x,OS_stage.x) ~ NMF_cluster + Signature_score, multi_var)
ggsurvplot_facet_risktable(fit_n_h_table, multi_var,pval=T,xlab = "Years")

fit_h_i_n <- survfit(Surv(OS_time.x,OS_stage.x) ~ Signature_score, multi_var)
ggsurvplot_facet(fit_h_i_n,facet.by = c("ICI","NMF_cluster"),data = multi_var,pval=T,xlab = "Years")

fit_i_h_n <- survfit(Surv(OS_time.x,OS_stage.x) ~ ICI, multi_var)
ggsurvplot_facet(fit_i_h_n,facet.by = c("Signature_score","NMF_cluster"),data = multi_var,pval=T,xlab = "Years")

fit_n_i_h <- survfit(Surv(OS_time.x,OS_stage.x) ~ NMF_cluster, multi_var)
ggsurvplot_facet(fit_n_i_h,facet.by = c("Signature_score","ICI"),data = multi_var,pval=T,xlab = "Years")

dev.off()

