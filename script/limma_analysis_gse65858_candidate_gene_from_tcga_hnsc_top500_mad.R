library('survival')
library('survminer')
#library("NMF")
library("limma")
library("ggrepel")
library('dplyr')

install.packages('NMF')
library("NMF")

setwd("working_dir")
############ load GSE65858 intersect data
GSE<-readRDS(file = "input/gse65858_hnsc_microarray.intersect.cli.rds")

######test for tcga top500 genes##
tcga_top500_t<-read.table("OUT/tcga_hnsc_mad_top500genes.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
tcga_top500<-as.character(tcga_top500_t$x)
#
rna_gse<-readRDS(file = "OUT/gse65858_hnsc_microarray.intersect.withoutcli.rds")
gse_t500<-rna_gse
gse_t500_t<-gse_t500[which(rownames(gse_t500)%in%as.character(tcga_top500)),]
saveRDS(gse_t500_t, file = "OUT/GSE65858_rna_tp500fromtcga.rds")

if(all(gse_t500_t>=0))
  print("checked") else
    (gse_t500_t<-(gse_t500_t*2)) & (gse_t500_t[gse_t500_t<0]<-0)

gse_t500_t<-as.matrix(gse_t500_t)

gse_t500_res2<-nmf(gse_t500_t, 2, nrun=100, seed=123456)
pdf(file="OUT/gse_t50065858_nmf_nrun100.allsample.top500genesfromtcga.consensusmap.split2.pdf",width=6,height=5.5)
consensusmap(gse_t500_res2, labCol=NA, labRow=NA)
dev.off()

gse_estim_tp500.r <- nmf(gse_t500_t, 2:6, nrun=100, seed=123456)
# Save an object to a file
saveRDS(gse_estim_tp500.r, file = "OUT/GSE65858_nmf_nrun100_estimate_2t06_tp500_fromtcga.rds")
# Restore the object
gse_rrr<-readRDS(file = "OUT/GSE65858_nmf_nrun100_estimate_2t06_tp500_fromtcga.rds")

pdf(file="OUT/GSE65858_nmf_nrun100_estimate_2t06_tp500_from_tcga.pdf",width=6,height=5)
plot(gse_rrr)
dev.off()

pdf(file="OUT/GSE65858_nmf_nrun100_estimate_2t06_tp500_consensusmap.pdf",width=13,height=8)
consensusmap(gse_rrr, labCol=NA, labRow=NA)
dev.off()

#saveRDS(object = gse_t500_res2, file = "OUT/gse65858_hnsc_microarray.intersect.gse_t500_nmf_res2.rds")
gse_t500_res2<-readRDS(file = "OUT/gse65858_hnsc_microarray.intersect.gse_t500_nmf_res2.rds")

gse_t500_classifier2 <- as.matrix(apply(coef(gse_t500_res2),2,which.max))
gse_t500_claf2<-as.data.frame(cbind(as.character(rownames(gse_t500_classifier2)),as.numeric(as.character(gse_t500_classifier2))))
colnames(gse_t500_claf2)<-c("SampleID","class")
gse_t500_claf2$group<-NA
gse_t500_claf2$group[which(gse_t500_claf2$class==1)]<-"AS1"
gse_t500_claf2$group[which(gse_t500_claf2$class==2)]<-"AS2"
gse_t500_claf2$group <- factor(gse_t500_claf2$group, levels = c("AS1","AS2"))

gse_t500_claf2<-merge(gse_t500_claf2,GSE[,1:3],by.x = "SampleID",by.y = "SampleID",all.x = T,sort = F)
#creat classifier input
gse_t500_claf2_classifier_input<-merge(gse_t500_claf2,GSE[,-c(2:3)],by.x="SampleID",by.y="SampleID",all.x=T, sort =F)
gse_t500_claf2_classifier_input$label<-NA
gse_t500_claf2_classifier_input$label[which(gse_t500_claf2_classifier_input$group=="AS1")]<-as.numeric(1)
gse_t500_claf2_classifier_input$label[which(gse_t500_claf2_classifier_input$group=="AS2")]<-as.numeric(0)
gse_t500_claf2_classifier_input<-as.data.frame(gse_t500_claf2_classifier_input)

gse_t500_claf2_classifier_input_save<-select(gse_t500_claf2_classifier_input,c("SampleID","class","group","OS_time","OS_stage"))

write.table(gse_t500_claf2_classifier_input_save,"OUT/tcga_hnsc_mad_top500genesfromTCGA_nmf2_cluster_samples.txt",sep = "\t")

gse_t500_claf2_classifier_input_all_genes<-select(gse_t500_claf2_classifier_input,-c("class","group","OS_time","OS_stage"))
#write.csv(x=gse_t500_claf2_classifier_input_all_genes, file="OUT/gse_t50065858_nmf_nrun100.allsample.top500genesfromtcga_all_genes_classifier_input.csv")
gse_t500_claf2_classifier_input_top500_genes<-select(gse_t500_claf2_classifier_input_all_genes,c("SampleID","label",tcga_top500))
#write.csv(x=gse_t500_claf2_classifier_input_top500_genes, file="OUT/gse_t50065858_nmf_nrun100.allsample.top500genesfromtcga_top500_genes_classifier_input.csv")
#

gse_t500_fit_claf2<-survfit(Surv(OS_time,OS_stage)~group,data = gse_t500_claf2)
pdf(file="OUT/gse_t50065858_nmf_nrun100.allsample.top500genesfromtcga.ggsurvplot.split2.pdf",width=6,height=6)
ggsurvplot(gse_t500_fit_claf2, data = gse_t500_claf2, pval = T,legend.title="gse_t500_65858 (n=270)",legend.labs = levels(gse_t500_claf2$group),
           legend = "top",surv.median.line="hv",risk.table = T,palette = c("#2E9FDF","#E7B800"))
dev.off()

###differential expression use limma###
exprSet<-gse_t500_t
  
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)
gsm_samp<-colnames(exprSet)
gsm_samp<-as.data.frame(gsm_samp)
ann_gse<-merge(gsm_samp,gse_t500_claf2,by.x = "gsm_samp",by.y = "SampleID",all.x = T,sort=F)
group<-as.factor(as.character(ann_gse$group))

design <- model.matrix(~0+group)
colnames(design)=levels(factor(group))
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts(paste0(unique(group),collapse = "-"),levels = design)
fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput <- topTable(fit2, coef=1, n=Inf)
nrDEG <- na.omit(tempOutput) 

nrDEG_t<-as.data.frame(nrDEG)
colnames(nrDEG_t)<-c("log2FoldChange","AveExpr","t","P.Value","padj","B")
#Volcano plot
alpha <- 0.01

gse_deg <- mutate(nrDEG_t, sig=ifelse(nrDEG_t$padj<alpha, "moderate", "Not Sig"))
gse_deg$sig[which(gse_deg$log2FoldChange > 1 & gse_deg$padj<0.01)]<-"up"
gse_deg$sig[which(gse_deg$log2FoldChange< -1 & gse_deg$padj<0.01)]<-"down"
vol_input <- cbind(gene=rownames(nrDEG), gse_deg)
vol_input<-na.omit(vol_input)

p<-ggplot(vol_input, aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=sig))+
  scale_color_manual(values=c("blue","black","gray","red"))+
  geom_hline(yintercept=-log10(alpha),color="brown")+
  geom_vline(xintercept=c(1,-1),color="brown")
pdf(file="OUT/limma_analysis_gse65858_candidate_gene_from_tcga_top500_mad.pdf",width=9,height=8)
p+geom_text_repel(data=vol_input[which(vol_input$padj<0.01 & abs(vol_input$log2FoldChange) >1),],
                  aes(label=gene))+scale_y_continuous(breaks = c(0,2,25,50,60),limits = c(0,60))+xlim(c(-4,2.5))
dev.off()

range(-log10(vol_input$padj))
#

