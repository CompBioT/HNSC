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
library("pheatmap")
library("circlize")
library("glmnet")

setwd("working_dir")
##before analysis, determine gene list##
GSE<-read.table("/input/GSE65858_exp.matrix.txt",
                sep = "\t",header = T,stringsAsFactors = F,check.names = F)

#part1, read tpm#
hnsc_tpm<-read.table("/input/HNSC_RSEM_genes/tcga_hnsc_tpm_os.matrix.txt",
                         sep = "\t",header = T,stringsAsFactors = F,check.names = F)

hnsc_tpm[1:4,1:4]

hnsc_tpm<-hnsc_tpm[,-20535]
cc<-colnames(hnsc_tpm)[which(colnames(hnsc_tpm) %in% colnames(GSE))]
cc<-as.character(cc)
hnsc_tpm<-hnsc_tpm[,as.character(cc)]

rownames(hnsc_tpm)<-hnsc_tpm$SampleID
hnsc_tpm$OS_time<-as.numeric(hnsc_tpm$OS_time)
hnsc_tpm$OS_time<-(hnsc_tpm$OS_time/365)
hnsc_tpm$OS_stage<-as.numeric(hnsc_tpm$OS_stage)

hnsc_tpm<-as.data.frame(hnsc_tpm)
dataEXP<-hnsc_tpm[4:ncol(hnsc_tpm)]

#top mad 8000 genes
dataEXP[,]<-sapply(dataEXP[,],as.numeric)
dataEXP<-dataEXP[,colSums(abs(dataEXP))!=0]
rna<-as.data.frame(t(dataEXP))
rna<-FSbyMAD(rna,cut.type = "topk",8000)

dataEXP<-as.data.frame(t(rna))
dataEXP<-scale(dataEXP,center=T)
data<-cbind(hnsc_tpm[,2:3],dataEXP)

##########
          ###############

# get survival related candidate genes########
coxR=data.frame()
coxf<-function(x){
  fmla1 <- as.formula(Surv(OS_time,OS_stage)~data[,x])
  mycox <- coxph(fmla1,data=data)
}

##part 2, select candidate prognosis related genes##
cores=detectCores()
cl<-makeCluster(cores[1]-1)
registerDoParallel(cl)

coxR<-foreach(index=3:ncol(data),.combine = rbind)%dopar%{
  library('survival')
  a=colnames(data)[index]
  mycox=coxf(a)
  coxResult=summary(mycox)
  cbind(id=a,HR=coxResult$coefficients[,"exp(coef)"],
        P=coxResult$coefficients[,"Pr(>|z|)"])
}

stopCluster(cl)
coxR<-as.data.frame(coxR)
write.table(coxR,"OUT/hnsc.tpm.coxR.txt",sep = "\t")
coxR<-read.table("OUT/hnsc.tpm.coxR.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
sigGene<-as.character(coxR[as.numeric(as.character(coxR$P))<0.01,]$id)
sigGene<-sigGene[!is.na(sigGene)]
write.table(sigGene,"OUT/hnsc.tpm.sigGene.txt",sep = "\t")
sigGene<-read.table("OUT/hnsc.tpm.sigGene.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
sigGene<-as.character(sigGene$x)

data1 <- cbind(data[,1:2],data[,as.character(sigGene)])

patients=rownames(hnsc_tpm)
##bootstrap for 1000 times to narrow down prognosis candidate genes##
outTab=data.frame()
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
outTab<-foreach(index=3:ncol(data1),.combine = rbind) %dopar% {
  library('survival')
  gene=colnames(data1)[index]
  Mboot <- replicate(1000, expr = {
    indices <- sample(patients, size=nrow(data1)*0.8, replace = F)
    data<-data1[indices,]
    fmla1 <- as.formula(Surv(data[,"OS_time"],data[,"OS_stage"])~data[,gene])
    mycox <- coxph(fmla1,data=data)
    coxResult = summary(mycox)
    P=coxResult$coefficients[,"Pr(>|z|)"]
  })
  times=length(Mboot[which(Mboot<0.05)])   #Mboot<0.05 fix
  t<-cbind(gene=gene,times=as.numeric(times))
}
stopCluster(cl)

##save genes which reach p<0.05 over 900 times##
outTab<-as.data.frame(outTab)
#write.table(outTab,"OUT/hgnc.tpm.outTab.timesp0.05.txt",sep = "\t") #fix at 20210226
#
outTab<-read.table("OUT/hgnc.tpm.outTab.timesp0.05.txt",
           sep = "\t",header = T,stringsAsFactors = F,check.names = F)

bootGene<-outTab[as.numeric(as.character(outTab$times))>900,]

#random forest to select top variable genes#
data2 <- cbind(data1[,1:2],data1[,as.character(bootGene$gene)])
#####
data4 <- data2
res.rsf<-rfsrc(Surv(OS_time,OS_stage)~., data4, nodesize = 20, proximity=T, tree.err = T, 
               forest = T, ntree = 1000, splitrule = "logrank", importance = TRUE,seed=123456)

plot(res.rsf)
#############
res.trc<-c()
res.trcoob<-c()
res.testc<-c()

topvars <-vector(mode="character",length=100) 
for (j in 1:100) { 
  print(paste("trying for",j,"times"))
  vars<-var.select(object=res.rsf,
                   cause =1,
                   method = "md", 
                   conservative = c("high"), 
                   ntree = 1000,
                   nodesize = 20, splitrule = "logrank", nsplit = 3, xvar.wt = NULL,
                   refit = T, fast = T,
                   na.action = c("na.impute"), 
                   always.use = NULL, nrep = 10,
                   prefit =  list(action = T, ntree = 1000,
                                  nodesize = 20, nsplit = 3),
                   verbose = TRUE)
  
  trc<-rcorr.cens(-vars$rfsrc.refit.obj$predicted, 
                  Surv(data4$OS_time, data4$OS_stage))["C Index"]
  trcoob<-rcorr.cens(-vars$rfsrc.refit.obj$predicted.oob, 
                     Surv(data4$OS_time, data4$OS_stage))["C Index"]
  
  res.trc<-rbind(res.trc, trc)
  res.trcoob<-rbind(res.trcoob, trcoob)
  topvars<-rbind(topvars,vars$topvars)
}

##
topvars_cp<-topvars
topvars_cp2<-topvars
write.table(topvars_cp,"OUT/hgnc.tpm.timesp0.05.randfost100times.topvars.txt",sep = "\t") #fix at 20210227
result<-data.frame(res.trc,res.trcoob,row.names = 1:nrow(res.trc))
colnames(result)<-c("res.trc.cindex","res.trcoob.cindex")
write.table(result,"OUT/hgnc.tpm.timesp0.05.randforst100times.trcoob.txt",sep = "\t")
bestresult<-result[result$res.trcoob.cindex==max(result$res.trcoob.cindex),] 
write.table(bestresult,"OUT/hgnc.tpm.timesp0.05.randforst100times.bestresult.txt",sep = "\t")

topvars_cp<-as.matrix(topvars_cp)[-1,] 
rownames(topvars_cp)<-c(1:dim(topvars_cp)[1])
bestvars<-unique(topvars_cp[rownames(bestresult),])
rd_vars<-bestvars
write.table(rd_vars,"OUT/hgnc.tpm.timesp0.05.randforst100times.best_variable.txt",sep = "\t")
############### 
               #######

#######lasso regression
rd_vars<-read.table("OUT/hgnc.tpm.timesp0.05.randforst100times.best_variable.txt",
           sep = "\t",header = T,stringsAsFactors = F,check.names = F)
rd_vars<-as.character(rd_vars$x)
#
data1 <- cbind(data[,1:2],data[,as.character(rd_vars)])
data_lso <- cbind(data1[,1:2],data1[,as.character(rd_vars)])
set.seed(123456)
x_1se<-data_lso[,-c(1:2)]
y_1se<-data_lso[,c(1:2)]
fit_1se<-Surv(y_1se$OS_time,y_1se$OS_stage)
gml_1se<-glmnet(as.matrix(x_1se),fit_1se,family = "cox",alpha = 1)
cvfit_1se<-cv.glmnet(as.matrix(x_1se),fit_1se,nfolds = 10,family="cox")

pdf("OUT/hgnc.tpm.timesp0.05.randforst100times.lasso_1se_15genes.pdf",width = 5.5,height = 6,paper = 'special')
plot(gml_1se,label=T,xvar="lambda")
plot(cvfit_1se)
print(cvfit_1se)

coef_1se<-coef(cvfit_1se,s="lambda.1se")

coef_1se<-as.data.frame(as.matrix(coef_1se))
colnames(coef_1se)<-"coef_1se"
genes_1se<-rownames(coef_1se)

coef_1se_nonzero<-coef_1se[which(coef_1se$coef_1se!=0),]

coef_1se[which(coef_1se$coef_1se!=0),]
coef_1se_gene<-genes_1se[which(coef_1se$coef_1se!=0)]
cbind(coef_1se_gene,round(coef_1se_nonzero, digits=4))

write.table(coef_1se_gene,"OUT/hgnc.tpm.timesp0.05.randfost100times_N_lasso_1se_15genes.txt",sep = "\t")

coef_1se_gene_train <- cbind(data[,1:2],data1[,as.character(coef_1se_gene)])

coef_1se_res.coxts<-coxph(Surv(OS_time,OS_stage)~ZNF266+SLC2A3+POLR1D+DDX3Y+CBWD3+
                            SC5DL+PGK1+CAMK2N1+AIG1+POMP+SEC11A+C9orf123+EFNB2+STC2+HPRT1,data=data1)


write.table(cbind(coef_1se_gene,coef_1se_res.coxts$coefficients),"OUT/hgnc.tpm.timesp0.05.randfost100times_N_lasso_1se_15genescoefficients.txt",sep = "\t")

coef_1se_TR<-data1[,as.character(coef_1se_gene)]
coef_1se_Score<-rowSums(sweep(coef_1se_TR,MARGIN = 2,coef_1se_nonzero,`*`))
Score<-coef_1se_Score
coef_1se_TR2<-as.data.frame(cbind(data1[,1:2],coef_1se_TR,Score))
##plot for genes that selected by lambda.1se
plot(coef_1se_res.coxts$coefficients,xlab="gene index",ylab="Cox coefficients")
coef_1se_TR2_roc<-timeROC(T=coef_1se_TR2$OS_time,delta=coef_1se_TR2$OS_stage,marker=coef_1se_TR2$Score,cause=1,weighting='marginal',ROC=T,iid=T,times=c(1,3,5))
plotAUCcurve(coef_1se_TR2_roc,conf.int = T,col = "tomato")

col<-c("#0073C2FF","firebrick1","orange","green","pink")
plot(coef_1se_TR2_roc,time=1,title=F,lwd=1.5,col=col[1])
plot(coef_1se_TR2_roc,time=3,title=F,lwd=1.5,col=col[2],add = T)
plot(coef_1se_TR2_roc,time=5,title=F,lwd=1.5,col=col[3],add = T)
id<-c(paste0("1-year AUC = ",round(coef_1se_TR2_roc$AUC[1],3)),
      paste0("3-year AUC = ",round(coef_1se_TR2_roc$AUC[2],3)),
      paste0("5-year AUC = ",round(coef_1se_TR2_roc$AUC[3],3))
)
legend("bottomright",id,fill=col[1:3],
       bty="0",cex=1,border = NA)
abline(0,1,lty=2,lwd=0.5)

dev.off()
quantile(coef_1se_TR2$Score)

pdf("OUT/hgnc.tpm.timesp0.05.randforst100times.lasso_1se_15genes.survival.pdf",width = 6,height = 6,paper = 'special')
coef_1se_TR2$group<-NA
coef_1se_TR2$group[which(coef_1se_TR2$Score<=median(coef_1se_TR2$Score))]<-"low"
coef_1se_TR2$group[which(coef_1se_TR2$Score>median(coef_1se_TR2$Score))]<-"high"
coef_1se_TR2$group <-factor(coef_1se_TR2$group,levels = c("low", "high"))

write.table(coef_1se_TR2,"OUT/hnsc_tpm_15gene_signature_scores.txt",sep = "\t")

coef_1se_TR2_fit_test<-survfit(Surv(OS_time,OS_stage)~group,data = coef_1se_TR2)
ggsurvplot(coef_1se_TR2_fit_test, data = coef_1se_TR2, pval = T,legend.title="Risk Score(TCGA-HNSC, n=519)",legend.labs = levels(coef_1se_TR2$group),
           legend = "top",surv.median.line="hv",risk.table = T,palette = c("#2E9FDF","#E7B800"))

coef_1se_TR2$groupQ<-NA
coef_1se_TR2$groupQ[which(coef_1se_TR2$Score<=quantile(coef_1se_TR2$Score)[4])]<-"Q1-Q3"
coef_1se_TR2$groupQ[which(coef_1se_TR2$Score>quantile(coef_1se_TR2$Score)[4])]<-"Q4"
coef_1se_TR2$groupQ <-factor(coef_1se_TR2$groupQ,levels = c("Q1-Q3", "Q4"))
coef_1se_TR2_fit_testQ<-survfit(Surv(OS_time,OS_stage)~groupQ,data = coef_1se_TR2)
ggsurvplot(coef_1se_TR2_fit_testQ, data = coef_1se_TR2, pval = T,legend.title="Risk Score(TCGA-HNSC, n=519)",legend.labs = levels(coef_1se_TR2$groupQ),
           legend = "top",surv.median.line="hv",risk.table = T,palette = c("#2E9FDF","#E7B800"))

coxph(Surv(OS_time, OS_stage) ~ groupQ, data = coef_1se_TR2) %>% 
  gtsummary::tbl_regression(exp = TRUE)
dev.off()

####plot risk score figure####

rk_dat<-select(coef_1se_TR2,c("OS_time","OS_stage","Score"))
rk_dat_sort<-rk_dat[order(rk_dat$Score),]
rk_dat_sort$index<-1:length(rk_dat$Score)
rk_dat_sort$Survival_state<-NA
rk_dat_sort$Survival_state[which(rk_dat_sort$OS_stage==0)]<-"Alive"
rk_dat_sort$Survival_state[which(rk_dat_sort$OS_stage==1)]<-"Death"
rk_dat_sort$Survival_state<-as.factor(rk_dat_sort$Survival_state)

spot<-max(rk_dat_sort$index[which(rk_dat_sort$Score<median(rk_dat_sort$Score))])

pdf("OUT/hgnc.tpm.timesp0.05.randforst100times.lasso_1se_riskscore_15genes.pdf",width = 7,height = 3,paper = 'special')
ggplot(rk_dat_sort,aes(x=index, y=Score,fill=Survival_state))+
  geom_bar(stat="identity") +xlab("index")+ylab("Risk Score")+
  geom_segment(aes(x=0,xend=spot,y=median(coef_1se_TR2$Score),yend=median(coef_1se_TR2$Score)),size=0.15,linetype=2)+
  geom_segment(aes(x = spot , y = -Inf, xend = spot, yend = median(coef_1se_TR2$Score)),size=0.15,linetype=2)+
  scale_x_continuous(breaks=seq(0,550,100))+scale_fill_manual(values=c("yellowgreen", "orangered2"))

ggplot(rk_dat_sort, aes(x=index, y=OS_time, color=Survival_state)) + 
  geom_point()+scale_color_manual(values=c("yellowgreen", "orangered2"))
dev.off()

#use TPM value for heatmap view#
nes<-select(coef_1se_TR2,c("OS_time","OS_stage","Score","group"))
heat_input<-cbind(nes,hnsc_tpm[,as.character(coef_1se_gene)])
heat_input$Risk_Group<-heat_input$group
heat_input_order<-with(heat_input, heat_input[order(group),])
ann_input<-select(heat_input_order,c("Risk_Group"))
ann_input<-as.data.frame(ann_input)
ann_input$Risk_Group<-as.factor(ann_input$Risk_Group)

heat_input_order<-select(heat_input_order,coef_1se_gene)
#check sample ID
all(rownames(heat_input_order)==rownames(ann_input))
heat_input_order<-t(heat_input_order)

library('gplots')
library('RColorBrewer')
library('devtools')

pdf("OUT/hgnc.tpm.timesp0.05.randforst100times.lasso_1se_15genes_heatmap.pdf",width = 7,height = 4,paper = 'special')
ggplot(rk_dat_sort, aes(x=index, y=OS_time, color=Survival_state)) + 
  geom_point()+scale_color_manual(values=c("yellowgreen", "orangered2"))

plot_color = c('green','orange')[ann_input$Risk_Group]
heatmap.2(heat_input_order,               #Input must be matrix
          trace="none",    
          scale="row",  
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)),
          density.info=c("none"),
          ColSideColors = plot_color,
          dendrogram = "none",
          Colv = F,
          breaks = seq(-2,2,0.2),
          labCol = F,
          key = T,
          keysize = 0.8
)

plot.new()
par(lend = 1)
legend("topright",
       legend = c("low", "high"),
       col = c("green", "orange"),
       lty= 1,
       lwd = 10
)

dev.off()


save.image(file='OUT/tcga_hnsc_tpm_random_survival_forest_image.RData')

