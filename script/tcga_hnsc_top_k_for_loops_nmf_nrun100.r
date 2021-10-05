library('survival')
library('randomForestSRC')
library('Hmisc')
library('dplyr')
library('CancerSubtypes')
library('doParallel')
library('caret')
library('survminer')
library('timeROC')
library('data.table')
library('dplyr')
library("DESeq2")
library('gProfileR')
library("affy")
library("ggrepel")

install.packages('NMF')
library("NMF")

setwd("working_dir")
##before analysis, determine gene list##
GSE<-read.table("/input/GSE65858_exp.matrix.txt",
                sep = "\t",header = T,stringsAsFactors = F,check.names = F)

#part1, read tpm#

hnsc_tpm<-read.table("/input/HNSC_RSEM_genes/tcga_hnsc_tpm_os.matrix.txt",
                     sep = "\t",header = T,stringsAsFactors = F,check.names = F)

hnsc_tpm<-hnsc_tpm[,-20535]
cc<-colnames(hnsc_tpm)[which(colnames(hnsc_tpm) %in% colnames(GSE))]
hnsc_tpm<-hnsc_tpm[,as.character(cc)]
##for backup
hnsc_tpm_cp<-hnsc_tpm

rownames(hnsc_tpm)<-hnsc_tpm$SampleID
hnsc_tpm$OS_time<-as.numeric(hnsc_tpm$OS_time)
hnsc_tpm$OS_time<-(hnsc_tpm$OS_time/365)
hnsc_tpm$OS_stage<-as.numeric(hnsc_tpm$OS_stage)

hnsc_tpm<-as.data.frame(hnsc_tpm)
#saveRDS(object = hnsc_tpm, file = "OUT/tcga_hnsc_tpm.cli.rds")
hnsc_tpm <- readRDS(file = "OUT/tcga_hnsc_tpm.cli.rds")
#
dataEXP<-hnsc_tpm[4:ncol(hnsc_tpm)]
dataEXP[,]<-sapply(dataEXP[,],as.numeric)
dataEXP<-dataEXP[,colSums(abs(dataEXP))!=0]

rna<-as.data.frame(t(dataEXP))
rna<-as.matrix(rna)
rna_cp<-rna
#saveRDS(rna_cp, file = "OUT/TCGA_HNSC_rna_cp.rds")

###RNA loop###
CI=data.frame()
pdf(file="OUT/TCGA_HNSC_nmf_nrun100.allsample.topkgenes.ggsurvplot.split2.pdf",width=5.5,height=6)
for(k in seq(100,5000,100)){
  rna_topgene<-FSbyMAD(rna_cp,cut.type = "topk",value = k)
  if(all(rna_topgene>=0))
    print("checked") else
      (rna_topgene<-(rna_topgene*2)) & (rna_topgene[rna_topgene<0]<-0)
  
  res_k <- nmf(rna_topgene, 2, nrun=100, seed=123456)
  classifier_k <- as.matrix(apply(coef(res_k),2,which.max))
  claf_k<-as.data.frame(cbind(as.character(rownames(classifier_k)),as.numeric(as.character(classifier_k))))
  colnames(claf_k)<-c("SampleID","class")
  claf_k$group<-NA
  claf_k$group[which(claf_k$class==1)]<-"AS1"
  claf_k$group[which(claf_k$class==2)]<-"AS2"
  claf_k<-merge(claf_k,hnsc_tpm[,1:3],by.x = "SampleID",by.y = "SampleID",all.x = T,sort = F)
  fit_clafk<-survfit(Surv(OS_time,OS_stage)~group,data = claf_k)
  pval_k<-surv_pvalue(fit_clafk, claf_k)$pval.txt
  
  cindex_k_c<-summary(coxph(Surv(OS_time,OS_stage)~group,data = claf_k))$concordance[1]
  cindex_k_sec<-summary(coxph(Surv(OS_time,OS_stage)~group,data = claf_k))$concordance[2]
  tit<-paste0("c-index=",cindex_k_c,";" ,"k=",";",k,"p-value=",pval_k)
  svp<-ggsurvplot(fit_clafk, data = claf_k, pval = T,legend.title="tit",legend.labs = levels(claf_k$group),
                  legend = "top",surv.median.line="hv",risk.table = F)
  

  info<-cbind(cindex_k_c,cindex_k_sec,k,pval_k)
  CI<-rbind(CI,info)
  print(svp)
}
dev.off()
CI
write.table(CI,"OUT/kgenes_cindex.txt",sep = "\t",row.names = F,append = F,quote=F)

###### NMF clustering k=2#
rna_tp500<-FSbyMAD(rna_cp,cut.type = "topk",value = 500)

tcga_top500<-as.character(rownames(rna_tp500))

if(all(rna_tp500>=0))
  print("checked") else
    (rna_tp500<-(rna_tp500*2)) & (rna_tp500[rna_tp500<0]<-0)

saveRDS(rna_tp500, file = "OUT/TCGA_HNSC_rna_tp500.rds")
rna_tp500<-readRDS(file="OUT/TCGA_HNSC_rna_tp500.rds")

estim_tp500.r <- nmf(rna_tp500, 2:6, nrun=100, seed=123456)
# Save an object to a file
saveRDS(estim_tp500.r, file = "OUT/TCGA_HNSC_nmf_nrun100_estimate_2t06_tp500.rds")
# Restore the object
rrr<-readRDS(file = "OUT/TCGA_HNSC_nmf_nrun100_estimate_2t06_tp500.rds")

pdf(file="OUT/TCGA_HNSC_nmf_nrun100_estimate_2t06_tp500.pdf",width=6,height=5)
plot(rrr)
dev.off()

pdf(file="OUT/TCGA_HNSC_nmf_nrun100_estimate_2t06_tp500_consensusmap.pdf",width=13,height=8)
consensusmap(rrr, labCol=NA, labRow=NA)
dev.off()
###set final top k = 500###
res2 <- nmf(rna_tp500, 2, nrun=100, seed=123456)
res2_cp<-res2

saveRDS(res2, file = "OUT/TCGA_HNSC_nmf_nrun100.top500genes.res2.rds")

pdf(file="OUT/TCGA_HNSC_nmf_nrun100.top500genes.estim.consensusmap.split2.pdf",width=6,height=5.5)
consensusmap(res2, labCol=NA, labRow=NA)
dev.off()

#evaluate the survival differences for top 500 MAD genes in TCGA HNSC RnaSeq data
res2<-readRDS(file = "OUT/TCGA_HNSC_nmf_nrun100.top500genes.res2.rds")
hnsc_tpm <- readRDS(file = "OUT/tcga_hnsc_tpm.cli.rds")
classifier2 <- as.matrix(apply(coef(res2),2,which.max))
claf2<-as.data.frame(cbind(as.character(rownames(classifier2)),as.numeric(as.character(classifier2))))
colnames(claf2)<-c("SampleID","class")
claf2$group<-NA
claf2$group[which(claf2$class==1)]<-"AS1"
claf2$group[which(claf2$class==2)]<-"AS2"
claf2$group <- factor(claf2$group, levels = c("AS1","AS2"))

claf2<-merge(claf2,hnsc_tpm[,1:3],by.x = "SampleID",by.y = "SampleID",all.x = T,sort = F)

claf2_classifier_input<-merge(claf2,hnsc_tpm[,-c(2:3)],by.x="SampleID",by.y="SampleID",all.x=T, sort =F)
claf2_classifier_input$label<-NA
claf2_classifier_input$label[which(claf2_classifier_input$group=="AS1")]<-as.numeric(1)
claf2_classifier_input$label[which(claf2_classifier_input$group=="AS2")]<-as.numeric(0)
claf2_classifier_input<-as.data.frame(claf2_classifier_input)

claf2_classifier_input_t<-select(claf2_classifier_input,-c("class","group","OS_time","OS_stage"))

#write.csv(x=claf2_classifier_input_t, file="OUT/tcga_hnsc_mad_top500genes_nmf2_cluster_samples_all_genes_classifier_input.csv")

claf2_classifier_input_top500_genes<-select(claf2_classifier_input_t,c("SampleID","label",all_of(tcga_top500)))
#write.csv(x=claf2_classifier_input_top500_genes, file="OUT/tcga_hnsc_mad_top500genes_nmf_nrun100.allsample.top500genesfromtcga_top500_genes_classifier_input.csv")

#write.table(claf2,"OUT/tcga_hnsc_mad_top500genes_nmf2_cluster_samples.txt",sep = "\t",row.names = F,append = F,quote=F)

fit_claf2<-survfit(Surv(OS_time,OS_stage)~group,data = claf2)
pdf(file="OUT/TCGA_HNSC_nmf_nrun100.allsample.top500genes.ggsurvplot.split2.pdf",width=6,height=6)
ggsurvplot(fit_claf2, data = claf2, pval = T,legend.title="TCGA HNSC (n=519)",
           legend = "top",surv.median.line="hv",risk.table = T,palette = c("#2E9FDF","#E7B800"))
dev.off()

#sava TCGA HNSC top 500 MAD genes
tcga_hnsc_mad_top500genes<-as.character(rownames(rna_tp500))
#write.table(tcga_hnsc_mad_top500genes,"OUT/tcga_hnsc_mad_top500genes.txt",sep = "\t",row.names = F,append = F,quote=F)


##############
##differential expression volcano plot##
##############
epsilon<-1
#tcga_top500<-read.table("OUT/tcga_hnsc_mad_top500genes.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
#tcga_top500<-as.character(tcga_top500$x)

raw_count<-read.table("input/HNSC_RSEM_genes/HNSC_RSEM_raw_count.txt",
                      sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
raw_count<-raw_count[,-567]

raw_count<-select(raw_count,-c("TCGA-CQ-A4CA-01A-11R-A24Z-07"))
tuid<-colnames(raw_count)[which(substr(colnames(raw_count),14,15)=="01")]

##TCGA-CQ-A4CA
raw_count<-select(raw_count,c(all_of(tuid)))
colnames(raw_count)<-substr(colnames(raw_count),1,12)

raw_count<-raw_count[which(rownames(raw_count)%in%tcga_top500),]
###
deg_input<-as.data.frame(raw_count)
deg_input[1:4,1:4]
rownames(deg_input)<-rownames(deg_input)

deg_samp<-colnames(deg_input)
deg_samp<-as.data.frame(deg_samp)

ann<-merge(deg_samp,claf2,by.x = "deg_samp",by.y = "SampleID",all.x = T,sort=F)
rownames(ann)<-ann$deg_samp
ann<-as.data.frame(ann)
ann<-select(ann,c("group"))
names(ann)[names(ann) == "group"] <- "Cluster"

ann$color<-NA
ann$color[which(ann$Cluster=="AS1")]<-"green"
ann$color[which(ann$Cluster=="AS2")]<-"orange"
ann$Cluster<-as.factor(ann$Cluster)
ann$Cluster<-factor(ann$Cluster, levels = c("AS1","AS2"))
ann$Cluster

head(ann)

getint <- function(x) {
  return(as.integer(x))
}
deg_input_int<-apply(deg_input, 2, getint)
deg_input_int<-as.data.frame(deg_input_int)
rownames(deg_input_int)<-rownames(deg_input)
##
#Differential analysis with DESeq2

deg_input_int[1:4,1:4]

#use the DESeqDateSetFromMatrix to create a DESeqDataSet object
dds0<-DESeqDataSetFromMatrix(countData = deg_input_int,colData = ann, design = ~ Cluster)
dds0$Cluster<-relevel(dds0$Cluster,"AS1")
#dds$condition <- relevel(dds$condition, "2")

print(dds0)
is(dds0)
isS4(dds0)
#what does it contain?
#the list of slot names
slotNames(dds0)
#cds is a countDataset
estimSf<-function(cds){
  cts<-counts(cds)
  #compute the geometric mean
  geomMean<-function(x) prod(x)^(1/length(x))
  #compute the geometric mean over the line
  gm.mean<-apply(cts, 1, geomMean)
  
  #Zero values are set to NA (avoid subsequentcdsdivision by 0)
  gm.mean[gm.mean == 0] <- NA
  #Divide each line by its corresponding geometric mean
  #sweep(x,MARGIN,STATS,FUN="-",check.margin=T, ...)
  #MARGIN: 1 or 2 (line or columns)
  #STATS: a vector of length nrow(x) or ncol(x),depending on MARGIN
  #FUN: the function to be applied
  cts<-sweep(cts, 1, gm.mean, FUN="/")
  #Compute the median over the columns
  med<-apply(cts,2,median,na.rm=T)
  return(med)
}
dds.norm<-estimateSizeFactors(dds0)
head(sizeFactors(dds.norm))
head(estimSf(dds0))
all(round(estimSf(dds0),6)==round(sizeFactors(dds.norm),6))

head(counts(dds.norm,normalized=T))[1:4,1:4]

#checking the normalization
par(mfrow=c(2,2),cex.lab=0.7)

mycol<-ann$color
pdf(file="OUT/tp1.pdf",width=8,height=100)
boxplot(log2(counts(dds.norm)+epsilon),  col=mycol, cex.axis=0.7, 
        las=1, xlab="log2(counts+1)", horizontal=TRUE, main="Raw counts")

boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=mycol, cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts")
dev.off()
plotDensity(log2(counts(dds.norm)+epsilon), col=mycol, 
            xlab="log2(counts+1)", cex.lab=0.7, panel.first=grid()) 

plotDensity(log2(counts(dds.norm)+epsilon),lty=1,col=ann$color,lwd=2)

plotDensity(log2(counts(dds.norm, normalized=TRUE)+epsilon), col=mycol, 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
# Restore default parameters
par(mfrow=c(1,1), cex.lab=1)

# Computing mean and variance
norm.counts <- counts(dds.norm, normalized=TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)

# sum(mean.counts==0) # Number of completely undetected genes
norm.counts.stats <- data.frame(
  min=apply(norm.counts, 2, min),
  mean=apply(norm.counts, 2, mean),
  median=apply(norm.counts, 2, median),
  max=apply(norm.counts, 2, max),
  zeros=apply(norm.counts==0, 2, sum),
  percent.zeros=100*apply(norm.counts==0, 2, sum)/nrow(norm.counts),
  perc05=apply(norm.counts, 2, quantile, 0.05),
  perc10=apply(norm.counts, 2, quantile, 0.10),
  perc90=apply(norm.counts, 2, quantile, 0.90),
  perc95=apply(norm.counts, 2, quantile, 0.95)
)

head(norm.counts.stats)
# Mean and variance relationship
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
#abline(a=8,b=1.4, col="brown")

# performing estimation of dispersion parameter
dds.disp<-estimateDispersions(dds.norm)
#A diagnostic plot which shows the mean of normalized counts (x axis) and dispersion estimate for each genes
plotDispEsts(dds.disp)

#Performing differential expression call
alpha <- 0.0001
waldTestResult <- nbinomWaldTest(dds.disp)

resultDESeq2 <- results(waldTestResult, alpha=alpha, pAdjustMethod="BH")
# What is the object returned by nbinomTest()
class(resultDESeq2)
is(resultDESeq2) # a data.frame
slotNames(resultDESeq2)  ## can view by resultDESeq2@
head(resultDESeq2)

# The column names of the data.frame
# Note the column padj 
# contains FDR values (computed Benjamini–Hochberg procedure)
colnames(resultDESeq2)
# Order the table by decreasing p-valuer
resultDESeq2 <- resultDESeq2[order(resultDESeq2$padj),]
head(resultDESeq2)
#ttt<-as.data.frame(resultDESeq2)
#
#Draw an histogram of the p-values
h1<-hist(resultDESeq2$pvalue,breaks=20,col="grey",main="DESeq2 p-value distribution",
         xlab="Nominal P-value",ylab="Number of genes")

#plot an horizontal bar with the estimation of m0 according to Storey-Tibshirani method (2003)
m<-sum(h1$counts) #Number of genes with a p-value
m0 <- 2*(sum(as.vector(na.omit(resultDESeq2$pvalue)) >0.5)) # Estimation of the number of genes under H0
m1 <- m - m0 # Estimation of the number of genes under m1
nbins <- length(h1$counts)
abline(v=m0/nbins, col="blue")

#Volcano plot
alpha <- 0.01 # Threshold on the p-value
# par(mfrow=c(1,2))
# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
resultDESeq2$sig <- -log10(resultDESeq2$padj)
sum(is.infinite(resultDESeq2$sig))
#change infinite p value
resultDESeq2[is.infinite(resultDESeq2$sig),"sig"] <- 350
#View(resultDESeq2[is.na(resultDESeq2$pvalue),])
# Select genes with a defined p-value (DESeq2 assigns NA to some genes)
genes.to.plot <- !is.na(resultDESeq2$pvalue)
range(resultDESeq2[genes.to.plot, "log2FoldChange"])

#Volcano plot of adjusteed p-values
cols<-densCols(resultDESeq2$log2FoldChange,resultDESeq2$sig)
cols[resultDESeq2$pvalue ==0] <- "purple"
resultDESeq2$pch <- 19
resultDESeq2$pch[resultDESeq2$pvalue ==0] <- 6
plot(resultDESeq2$log2FoldChange, 
     resultDESeq2$sig, 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=resultDESeq2$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
gn.selected <- abs(resultDESeq2$log2FoldChange) > 1 & resultDESeq2$padj < alpha 

text(resultDESeq2$log2FoldChange[gn.selected],
     -log10(resultDESeq2$padj)[gn.selected],
     lab=rownames(resultDESeq2)[gn.selected], cex=0.6,pos=2)


tttt<-as.data.frame(resultDESeq2)
mutateddf <- mutate(tttt, sig=ifelse(tttt$padj<0.01, "moderate", "Not Sig"))

mutateddf$sig[which(mutateddf$log2FoldChange>1 & mutateddf$padj<0.01)]<-"up"
mutateddf$sig[which(mutateddf$log2FoldChange< -1 & mutateddf$padj<0.01)]<-"down"

input <- cbind(gene=rownames(tttt), mutateddf)
input<-na.omit(input)

p<-ggplot(input, aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=sig))+
  scale_color_manual(values=c("blue","black","gray","red"))+
  geom_hline(yintercept=-log10(alpha),color="brown")+
  geom_vline(xintercept=c(1,-1),color="brown")
p+geom_text_repel(data=input[which(input$padj<0.01 & abs(input$log2FoldChange) >1 ),],
                    aes(label=gene))+scale_y_continuous(breaks = c(0,2,25,50,75,100),limits = c(0,100))

#
pdf(file="OUT/TCGA_HNSC_nmf_nrun100.allsample.top500genes.DESEQ2.pdf",width=9,height=8)
p+geom_text_repel(data=input[which(input$padj<0.01 & abs(input$log2FoldChange) >1 ),],
                  aes(label=gene))
dev.off()

write.table(input,"OUT/TCGA_AS1_AS2_differential_expression_result.txt",sep = "\t",row.names = F,append = F,quote=F)
