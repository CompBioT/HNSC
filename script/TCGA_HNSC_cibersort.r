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

##before analysis, determine gene list##
setwd("working_dir")
source('input/cibersort/Cibersort.R')
GSE<-read.table("input/GSE65858_exp.matrix.txt",
                sep = "\t",header = T,stringsAsFactors = F,check.names = F)

hnsc_tpm<-read.table("input/HNSC_RSEM_genes/tcga_hnsc_tpm_os.matrix.txt",
                     sep = "\t",header = T,stringsAsFactors = F,check.names = F)
hnsc_tpm<-hnsc_tpm[,-20535] ##delect last empty column

cc<-colnames(hnsc_tpm)[which(colnames(hnsc_tpm) %in% colnames(GSE))]
cc<-as.character(cc)
hnsc_tpm<-hnsc_tpm[,as.character(cc)]
rownames(hnsc_tpm)<-hnsc_tpm$SampleID
hnsc_tpm<-as.data.frame(hnsc_tpm)
hnsc_tpm_for_ciber_pre<-hnsc_tpm[5:ncol(hnsc_tpm)]
hnsc_tpm_for_ciber_pre[,]<-sapply(hnsc_tpm_for_ciber_pre[,],as.numeric)
hnsc_tpm_for_ciber_pre<-hnsc_tpm_for_ciber_pre[,colSums(abs(hnsc_tpm_for_ciber_pre))!=0]
hnsc_tpm_for_ciber<-as.data.frame(t(hnsc_tpm_for_ciber_pre))
colnames(hnsc_tpm_for_ciber)<-as.character(colnames(hnsc_tpm_for_ciber))
write.table(hnsc_tpm_for_ciber,"input/cibersort/hnsc_tpm_for_cibersort.txt",sep = "\t",quote = F)
fiber_res <- CIBERSORT('LM22.txt','hnsc_tpm_for_cibersort.txt', perm = 1000, QN = T)
fiber_res_backup<-fiber_res
saveRDS(fiber_res_backup, file = "/OUT/TCGA_HNSC_CIBERSORT_result.rds")
fiber_res<-readRDS(file="/OUT/TCGA_HNSC_CIBERSORT_result.rds")

HNSC_ciber <- fiber_res[,1:22]
HNSC_ciber_mat<-as.matrix(HNSC_ciber)
samlist<-rownames(HNSC_ciber_mat)
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

scores_hnsc=read.table("/OUT/HNSC_estimate_score.gct",skip = 2,header = T,check.names = F)
rownames(scores_hnsc)=scores_hnsc[,1]
colnames(scores_hnsc)<-gsub("\\.", "-",colnames(scores_hnsc))
scores_hnsc=t(scores_hnsc[,3:ncol(scores_hnsc)])

scores_est<-scores_hnsc
scores_est<-as.matrix(scores_est)
scores_est<-as.data.frame(scores_est)
scores_est$sampleID<-NA
scores_est$sampleID <- gsub("\\.", "-",rownames(scores_est))

###15-gene signature score###
gene_score_15<-read.table("/OUT/hnsc_tpm_15gene_signature_scores.txt",
                          sep = "\t",header = T,stringsAsFactors = F,check.names = F)

sig_scores<-cbind(rownames(gene_score_15),gene_score_15$group)
sig_scores<-as.data.frame(sig_scores)
colnames(sig_scores)<-c("sampleID","group")

nmf2_cluster<-read.table("/OUT/tcga_hnsc_mad_top500genes_nmf2_cluster_samples.txt",
                         sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(nmf2_cluster)<-c("sampleID","class","NMF_cluster","OS_time","OS_stage")

#####heatmap for hierachical clustering####
scaled_mat = t(scale(HNSC_ciber_mat))

#####
te0<-Heatmap(scaled_mat,column_km_repeats = 500,
             column_split=3,
             show_column_names=F,
             clustering_distance_rows = "euclidean",
             clustering_method_rows = "ward.D",
             clustering_distance_columns = "euclidean",
             clustering_method_columns = "ward.D",
             #cluster_rows = F,
)

ht = draw(te0)

ann0<-merge(samlist,scores_est,by.x="ann_samp",by.y="sampleID",all.x=T)
ann0$Subtype<-NA
ann0$Subtype[column_order(ht)[[1]]]<-"C1"
ann0$Subtype[column_order(ht)[[2]]]<-"C2"
ann0$Subtype[column_order(ht)[[3]]]<-"C3"
ann1<-merge(ann0,sig_scores,by.x="ann_samp",by.y="sampleID",all.x=T)
######
rcc6 = ConsensusClusterPlus(scaled_mat,maxK=6,reps=50,
                            pItem=0.8,pFeature=1,title="TCGA_HNSC",
                            distance='euclidean',
                            clusterAlg="pam",plot="png")

consensusClass <- rcc6[[2]][["consensusClass"]]

#three consensusClasses found no significant in TCGA-HNSC cohort
#consensusClass <- rcc6[[3]][["consensusClass"]]
consensusClass<-as.data.frame(consensusClass)
consensusClass$sampleID<-NA
consensusClass$sampleID<-rownames(consensusClass)

ann2<-merge(ann1,nmf2_cluster,by.x="ann_samp",by.y="sampleID",all.x=T)
ann3<-merge(ann2,consensusClass,by.x="ann_samp",by.y="sampleID",all.x=T)
##############
setwd("/OUT")
ann3$ICI<-NA
ann3$ICI[which(ann3$consensusClass==2)]<-"ICI_A"
ann3$ICI[which(ann3$consensusClass==1)]<-"ICI_B"
ann3$ICI<-factor(ann3$ICI,levels = c("ICI_A", "ICI_B"))

pdf("HNSC_cibersort_consensusClass_2clusters_survival.pdf",width=7,height=5)
fit_immune_os<-survfit(Surv(OS_time,OS_stage)~ICI,data = ann3)
ggsurvplot(fit_immune_os, data = ann3, pval = T,legend.title = "Immune_subtype",
           legend = "top",surv.median.line="hv",risk.table = T)
dev.off()

pdf("HNSC_cibersort_consensusClass_3clusters_survival.pdf",width=8,height=6)
fit_immune_os_cens3<-survfit(Surv(OS_time,OS_stage)~consensusClass,data = ann3)
ggsurvplot(fit_immune_os_cens3, data = ann3, pval = T,legend.title = "Immune_subtype",
           legend = "top",surv.median.line="hv",risk.table = T)
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
pdf(file="HNSC_cibersort_consensusClass_boxplot.pdf",width=6,height=3)
ggarrange(p1, p2, p3, 
          ncol = 3, nrow = 1)
dev.off()

##22 cell type plots##
library('tidyr')
library('tidyverse')
library(reshape2)
plot_info_raw<-HNSC_ciber_mat
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
getwd()

#http://www.sci666.net/51296.html
pdf(file="HNSC_cibersort_consensusClass_cells_boxplot.pdf",width=8,height=5.5)
ggboxplot(md,x = "variable",y = "value",color = "black",
          fill = "ICI",
          ylab = "fraction",
          xlab = "cells",
          main = "TME Cell composition (TCGA_HNSC)")+
  stat_compare_means(label = "p.signif",
                     method = "t.test",
                     ref.group = ".all.",
                     hide.ns = T)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5,size = 10))

dev.off()

####
pdf("HNSC_cibersort_Heatmap_split_by_riskscore.pdf",width=11,height=7)
Heatmap(scaled_mat,column_km_repeats = 500,name = " ",
        top_annotation = HeatmapAnnotation(
          StromalScore=ann3$StromalScore,
          ImmuneScore=ann3$ImmuneScore,
          ESTIMATEScore=ann3$ESTIMATEScore,
          NMF_cluster=ann3$NMF_cluster,
          Risk_score=ann3$group,
          ICI=ann3$ICI,
          col=list(
            StromalScore=col_stromal,
            ImmuneScore=col_imm,
            ESTIMATEScore=col_est,
            NMF_cluster=c("AS1"="cadetblue1","AS2"="#E7B800"),
            Risk_score=c("low"="cadetblue1","high"="#E7B800"),
            ICI=c("ICI_B"="cadetblue1","ICI_A"="#E7B800")
          )
        ),
        column_split = ann3$group,
        show_column_names=F,
)
dev.off()

######