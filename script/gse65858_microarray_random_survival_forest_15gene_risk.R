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
library("gtsummary")

#part1, read gse65858 microarray data#
setwd("/setwd("working_dir")")
col<-c("#0073C2FF","firebrick1","orange","green","pink")
GSE<-read.table("input/GSE65858_exp.matrix.txt",
                          sep = "\t",header = T,stringsAsFactors = F,check.names = F)
GSE<-GSE[,-22006]
GSE<-GSE[,cc]
GSE_raw<-GSE

ARR<-GSE[4:ncol(GSE)]
ARR<-scale(ARR)
GSE<-cbind(GSE[,1:3],ARR)
GSE$OS_time<-(as.numeric(GSE$OS_time)/365)

###by random forest,times0.05###
rfls_gse_gene_test <- cbind(GSE[,2:3],GSE[,as.character(coef_1se_gene)])
rfls_gse_TS<-GSE[,as.character(coef_1se_gene)]  #rfls_gse_res.coxts$coefficients

gse_coef<-rowSums(sweep(data.matrix(rfls_gse_TS),MARGIN = 2,coef_1se_res.coxts$coefficients,`*`))
gse_coef_1se<-rowSums(sweep(data.matrix(rfls_gse_TS),MARGIN = 2,coef_1se_nonzero,`*`))

rfls_gse_TS2<-as.data.frame(cbind(rfls_gse_gene_test[,1:2],rfls_gse_TS,gse_coef,gse_coef_1se))

#rfls_gse_TS2$Score<-rfls_gse_TS2$gse_coef
rfls_gse_TS2$Score<-rfls_gse_TS2$gse_coef_1se

rownames(rfls_gse_TS2)<-as.character(GSE$SampleID)
plot(coef_1se_res.coxts$coefficients)
###plot final result###
pdf("OUT/gse65858.timesp0.05.randforst100times.lasso_1se_15genes.pdf",width = 5.5,height = 6,paper = 'special')
rfls_gse_ts2_roc<-timeROC(T=rfls_gse_TS2$OS_time,delta=rfls_gse_TS2$OS_stage,marker=rfls_gse_TS2$Score,cause=1,weighting='marginal',ROC=T,iid=T,times=c(1,3,5))
plotAUCcurve(rfls_gse_ts2_roc,conf.int = T,col = "tomato")
plot(rfls_gse_ts2_roc,time=1,title=F,lwd=1.5,col=col[1])
plot(rfls_gse_ts2_roc,time=3,title=F,lwd=1.5,col=col[2],add = T)
plot(rfls_gse_ts2_roc,time=5,title=F,lwd=1.5,col=col[3],add = T)
id<-c(paste0("1-year AUC = ",round(rfls_gse_ts2_roc$AUC[1],3)),
      paste0("3-year AUC = ",round(rfls_gse_ts2_roc$AUC[2],3)),
      paste0("5-year AUC = ",round(rfls_gse_ts2_roc$AUC[3],3))
)
legend("bottomright",id,fill=col[1:3],
       bty="0",cex=1,border = NA)

abline(0,1,lty=2,lwd=0.5)
summary(rfls_gse_TS2$Score)
dev.off()

pdf("OUT/gse65858.timesp0.05.randforst100times.lasso_1se_15genes.survival.pdf",width = 6,height = 6,paper = 'special')
rfls_gse_TS2$group<-NA
rfls_gse_TS2$group[which(rfls_gse_TS2$Score<=median(rfls_gse_TS2$Score))]<-"low"
rfls_gse_TS2$group[which(rfls_gse_TS2$Score>median(rfls_gse_TS2$Score))]<-"high"
rfls_gse_TS2$group <-factor(rfls_gse_TS2$group,levels = c("low", "high"))
rfls_gse_TS2_fit_test<-survfit(Surv(OS_time,OS_stage)~group,data = rfls_gse_TS2)
ggsurvplot(rfls_gse_TS2_fit_test, data = rfls_gse_TS2, pval = T,legend.title="Risk Score(GSE65858, n=270)",legend.labs = levels(rfls_gse_TS2$group),
           legend = "top",surv.median.line="hv",risk.table = T,palette = c("#2E9FDF","#E7B800"))

#rfls_gse_TS2$group <-factor(rfls_gse_TS2$group,levels = c("low", "high"))
coxph(Surv(OS_time, OS_stage) ~ group, data = rfls_gse_TS2) %>% 
  gtsummary::tbl_regression(exp = TRUE)

rfls_gse_TS2$groupQ<-NA
rfls_gse_TS2$groupQ[which(rfls_gse_TS2$Score<=quantile(rfls_gse_TS2$Score)[4])]<-"Q1-Q3"
rfls_gse_TS2$groupQ[which(rfls_gse_TS2$Score>quantile(rfls_gse_TS2$Score)[4])]<-"Q4"
rfls_gse_TS2$groupQ <-factor(rfls_gse_TS2$groupQ,levels = c("Q1-Q3", "Q4"))

rfls_gse_TS2_fit_testQ<-survfit(Surv(OS_time,OS_stage)~groupQ,data = rfls_gse_TS2)
ggsurvplot(rfls_gse_TS2_fit_testQ, data = rfls_gse_TS2, pval = T,legend.title="Risk Score(GSE65858, n=270)",legend.labs = levels(rfls_gse_TS2$groupQ),
           legend = "top",surv.median.line="hv",risk.table = T,palette = c("#2E9FDF","#E7B800"))

coxph(Surv(OS_time, OS_stage) ~ groupQ, data = rfls_gse_TS2) %>% 
  gtsummary::tbl_regression(exp = TRUE)

dev.off()

####plot risk score figure####
rk_geo<-select(rfls_gse_TS2,c("OS_time","OS_stage","Score"))
rk_geo_sort<-rk_geo[order(rk_geo$Score),]

rk_geo_sort$index<-1:length(rk_geo$Score)
rk_geo_sort$Survival_state<-NA
rk_geo_sort$Survival_state[which(rk_geo_sort$OS_stage==0)]<-"Alive"
rk_geo_sort$Survival_state[which(rk_geo_sort$OS_stage==1)]<-"Death"
rk_geo_sort$Survival_state<-as.factor(rk_geo_sort$Survival_state)

spot_gse<-max(rk_geo_sort$index[which(rk_geo_sort$Score<median(rk_geo_sort$Score))])

pdf("OUT/gse65858.timesp0.05.randforst100times.lasso_1se_15genes_riskscore.pdf",width = 7,height = 3,paper = 'special')
ggplot(rk_geo_sort,aes(x=index, y=Score,fill=Survival_state))+
  geom_bar(stat="identity") +xlab("index")+ylab("Risk Score")+
  geom_segment(aes(x=0,xend=spot_gse,y=median(rk_geo_sort$Score),yend=median(rk_geo_sort$Score)),size=0.15,linetype=2)+
  geom_segment(aes(x = spot_gse , y = -Inf, xend = spot_gse, yend = median(rk_geo_sort$Score)),size=0.15,linetype=2)+
  scale_x_continuous(breaks=seq(0,550,100))+scale_fill_manual(values=c("yellowgreen", "orangered2"))

ggplot(rk_geo_sort, aes(x=index, y=OS_time, color=Survival_state)) + 
  geom_point()+scale_color_manual(values=c("yellowgreen", "orangered2"))

dev.off()
nes_gse<-select(rfls_gse_TS2,c("OS_time","OS_stage","Score","group"))
write.table(nes_gse,"OUT/GSE65858_15gene_signature_scores.txt",sep = "\t")
heat_gse_input<-cbind(nes_gse,GSE_raw[,as.character(coef_1se_gene)])
heat_gse_input$Risk_Group<-heat_gse_input$group
heat_gse_input_order<-with(heat_gse_input, heat_gse_input[order(group),])
ann_input<-select(heat_gse_input_order,c("Risk_Group"))
ann_input<-as.data.frame(ann_input)
ann_input$Risk_Group<-as.factor(ann_input$Risk_Group)

heat_gse_input_order<-select(heat_gse_input_order,coef_1se_gene)
#check sample ID
#all(rownames(heat_gse_input_order)==rownames(ann_input))
heat_gse_input_order<-t(heat_gse_input_order)

pdf("OUT/gse65858.timesp0.05.randforst100times.lasso_1se_15genes_pheatmap.pdf",width = 7,height = 4,paper = 'special')
plot_color = c('green','orange')[ann_input$Risk_Group]
heatmap.2(heat_gse_input_order,               #Input必须是matrix
          trace="none",    # trace可以给每个色块中添加一条线，与行平行或者与列平行。其与色块中心的距离代表了这个值被显示的比例。
          scale="row",  # scale在这里。
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)),
          density.info=c("none"),
          ColSideColors = plot_color,
          dendrogram = "none",
          Colv = F,
          breaks = seq(-2,2,0.2),
          labCol = F,
          key = T,
          #margins=c(6,12),
          keysize = 0.8
)

dev.off()

save.image(file='OUT/gse65858_microarray_random_survival_forest_image.RData')
