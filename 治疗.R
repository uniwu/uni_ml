#----------新辅助治疗预测价值------------
##-----------------GSE16446----------
library(AnnoProbe)
library(GEOquery)
library(affy)
library(Biobase)
library(ggExtra)
library(ggpubr)
library(limma)
gset <- getGEO('GSE16446', destdir="D:/RStudio/Ful数据/F3/测试集/GSE16446",
               AnnotGPL = T,     
               getGPL = T) 
GPL<-fData(gset[[1]])
pdata<-pData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$Gene_Symbol <- data.frame(sapply(gpl$Gene_Symbol, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
#代入模型计算评分
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))
GSE16446<-exprset_unique[comsam,]
GSE16446<-t(GSE16446)
comsam<-intersect(rownames(GSE16446os),rownames(GSE16446))
GSE16446<-GSE16446[comsam,]
GSE16446<-scale(GSE16446)
GSE16446<-cbind(GSE16446os,GSE16446)
GSE16446$OS<-as.numeric(GSE16446$OS)
GSE16446$OS.time<-as.numeric(GSE16446$OS.time)
rsf_v<- predict(fit,newdata = GSE16446,proximity=T)
score_16446<- data.frame(GSE16446[,c(1,2)],Score=rsf_v$predicted)
pcr <- pdata[, 61, drop = FALSE]
pcr<-pcr[comsam,,drop = FALSE]
score_16446<-cbind(score_16446,pcr)
#pam50
BiocManager::install("genefu") # 包
BiocManager::install("breastCancerTRANSBIG") # 数据集
library(genefu)
library(breastCancerTRANSBIG)
## 内置测试数据
data(pam50.robust)
dd <- get(data(transbig))
ddata <- t(exprs(dd))
dannot <- featureData(dd)@data
annot_data <- dannot %>% dplyr::select(EntrezGene.ID,probe,Gene.symbol)
#数据处理
exprset_unique$Gene.symbol<-rownames(exprset_unique)
exp1 <- dplyr::inner_join(annot_data,exprset_unique)
exp2 <- exp1[,-c(1:3)]
exp3 <- t(exp2)
colnames(exp3) <- exp1$probe
exp_anno <- exp1[,c(1:3)]
exp_pam50 <- molecular.subtyping(sbt.model = "pam50",data=exp3,annot=exp_anno,do.mapping=TRUE)
pam<-as.data.frame(exp_pam50$subtype)
pam<-pam[comsam,,drop = FALSE]
score_16446<-cbind(score_16446,pam)
colnames(score_16446)[4]<-"pcr"
colnames(score_16446)[5]<-"pam"
save(score_16446,file = 'score16446.RData')
#寻找最佳截断
score_16446<-score_16446[order(score_16446$pam), ]
# 58     12      6     20     11 
basal16446<-score_16446[(1:58),]
library(maxstat)
basal16446$pcr<-as.factor(basal16446$pcr)
ms <- maxstat.test(pcr~Score,
                   data = basal16446,
                   pmethod="HL")# 106.9019 
ms
# 添加一个新列以标记Score高低
basal16446<- basal16446 %>%
  group_by(pam) %>%
  mutate(Group = ifelse(Score >  106.9019 , "high", "low"))
# 计算每个组的pCR率
library(dplyr)
score_16446$pcr <- as.numeric(score_16446$pcr)


ms
pcr_rates <- score_16446 %>%
  group_by(pam, Group) %>%
  summarise(pCR_rate = mean(pcr, na.rm = TRUE) * 100, .groups = 'drop')
pcr_rates <- basal16446 %>%
  group_by(pam, Group) %>%
  summarise(pCR_rate = mean(pcr, na.rm = TRUE) * 100, .groups = 'drop')
#plot





p <- ggplot(basal16446, aes(x = factor(pcr), y = Score, color = factor(pcr))) +
  geom_boxplot(aes(fill = factor(pcr)), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_16446$Score),max(score_16446$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
p


##---------GSE22226------------
gset <- getGEO('GSE22226', destdir="D:/RStudio/Ful数据/训练6/F6/gse22226",
               AnnotGPL = T,     
               getGPL = T) 
GPL<-fData(gset[[1]])
pdata<-pData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$Gene_Symbol <- data.frame(sapply(gpl$Gene_Symbol, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))
#代入模型计算评分

GSE22226<-exprset_unique[comsam,]
GSE4779<-t(GSE4779)
GSE4779<-scale(GSE4779)




rsf_v<- predict(fit,newdata = trainrsf,proximity=T)
score_4779<- data.frame(trainrsf[,c(1,2)],Score=rsf_v$predicted)
pcr <- pdata[, 61, drop = FALSE]
pcr<-pcr[comsam,,drop = FALSE]
score_4779<-cbind(score_4779,pcr)
#--------辅助治疗-------
##---------20685---------
ad<-pdata[,48,drop=F]
colnames(ad)<-'adjuvant_chemo'
ad$adjuvant_chemo<- ifelse(ad$adjuvant_chemo== "yes", 1, ifelse(ad$adjuvant_chemo== "no", 0,ad$adjuvant_chemo))
# 删除无效信息
ad <- ad[ad$adjuvant_chemo != "unknown", ,drop=F]
comsam<-intersect(rownames(score_20685),rownames(ad))
score_20685<-score_20685[comsam,]
ad<-ad[comsam,,drop=F]
score_20685<-cbind(ad,score_20685)

##------37751----------
therapy<-pdata[,c(46,49,50),drop=F]
comsam<-intersect(rownames(therapy),rownames(score_37751))
therapy<-therapy[comsam,,drop=F]
therapy$`chemotherapy:ch1`<- ifelse(therapy$`chemotherapy:ch1`== "Yes", 1, ifelse(therapy$`chemotherapy:ch1`== "No", 0,therapy$`chemotherapy:ch1`))
therapy$`hormone-therapy:ch1`<- ifelse(therapy$`hormone-therapy:ch1`== "Yes", 1, ifelse(therapy$`hormone-therapy:ch1`== "No", 0,therapy$`hormone-therapy:ch1`))
therapy$`neoadjuvant therapy:ch1`<- ifelse(therapy$`neoadjuvant therapy:ch1`== "Yes", 1, ifelse(therapy$`neoadjuvant therapy:ch1`== "No", 0,therapy$`neoadjuvant therapy:ch1`))
therapy <- therapy[!(apply(therapy, 1, function(row) any(row == '--'))), ]
#所有列转为数值型
therapy[] <- lapply(therapy, as.numeric)
therapy$number_of_treatments <- rowSums(therapy)
# 标记接受0、1、2或3种治疗的病人
therapy$category <- ifelse(therapy$number_of_treatments == 0, 0,
                      ifelse(therapy$number_of_treatments == 1, 1,
                             ifelse(therapy$number_of_treatments == 2, 2, 3)))
score_37751a<-score_37751[comsam,,drop=F]
score_37751a<-cbind(score_37751a,therapy$category)
colnames(score_37751a)[ncol(score_37751a)] <- "category"
#plot
score_37751a$category <- as.factor(score_37751a$category)
score_37751a$category[score_37751a$category == 1] <- 0
pdf("boxplot_metastasis20685.pdf", width = 6, height = 5.5)
ggplot(score_37751a, aes(x = category,y =Score, color =category)) +
  geom_boxplot(aes(fill =category), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2",'#00A087B2' )) +  # 设置点的颜色
  scale_fill_manual(values = c("#8491B4B2","#E64B35B2",'#00A087B2')) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_37751a$Score),max(score_37751a$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()




#--------复发--------
##----GSE9893-------
gset <- getGEO('GSE9893', destdir="D:/RStudio/Ful数据/训练6/F6/GSE9893",
               AnnotGPL = F,     
               getGPL = T) 
GPL<-fData(gset[[1]])
pdata<-pData(gset[[1]])
gpl<-GPL[,c(1,13)]
gpl$Gene_Symbol <- data.frame(sapply(gpl$Gene_Symbol, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$Gene_Symbol))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$Gene_Symbol)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))
##------GSE69031--------
gset <- getGEO('GSE69031', destdir="D:/RStudio/Ful数据/训练6/F6/GSE69031",
               AnnotGPL = T,     
               getGPL = T) 
GPL<-fData(gset[[1]])
pdata<-pData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<- data.frame(sapply(gpl$`Gene symbol`, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
#部分样本没有表达值，删去NA列
exprset <- exprset[, colSums(is.na(exprset)) == 0]
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
#代入模型计算评分
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))
##------GSE20685---------
setwd("D:/RStudio/Ful数据/训练6/F6")
gset <- getGEO('GSE20685', destdir="D:/RStudio/Ful数据/F3/测试集/GSE20685",
               AnnotGPL = T,     
               getGPL = T) 
pdata<-pData(gset[[1]])

rr<-pdata[,58,drop = FALSE]
score_20685<-cbind(score_20685,rr)
colnames(score_20685)[ncol(score_20685)]<-'reginal_relapse'
score_20685<- score_20685[score_20685$reginal_relapse != "NA", ]
score_20685$reginal_relapse<-as.factor(score_20685$reginal_relapse)
save(score_20685,file="score_20685.RData")
ms <- maxstat.test(reginal_relapse~Score,
                   data = score_20685,
                   pmethod="HL")

ms#113.7603 
score_20685$group <- ifelse(score_20685$Score >=113.7603, "high", "low")
p <- ggplot(score_20685, aes(x = reginal_relapse, y =Score, color = reginal_relapse)) +
  geom_boxplot(aes(fill =reginal_relapse), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_20685$Score),max(score_20685$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
p



#--------ki67---------
##-----88770------
#根据score列将score_88770分为high和low两组
library(maxstat)
ms <- maxstat.test(ki67~Score,
                   data = score_88770,
                   pmethod="HL")

ms#109.6904 
score_88770$group <- ifelse(score_88770$Score >=109.6904, "high", "low")
score_88770$`ki67:ch1`<-as.numeric(score_88770$`ki67:ch1`)
score_88770<- na.omit(score_88770)
colnames(score_88770)[1]<-'ki67'
##----km-----
survifit <- survfit(Surv(OS.time,OS)~group,score_88770)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~group,score_88770)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
pdf("KM_88770_ki67.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data =score_88770, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("high", "low"),
  xlab = "Time (years)",
  break.time.by =5,
  palette = c("#E64B35B2", "#8491B4B2"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
##----boxplot----
pdf("boxplot_88770_ki67.pdf", width = 6, height = 5.5)
ggplot(score_88770, aes(x = group, y =ki67, color = group)) +
  geom_boxplot(aes(fill = group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values =  c("#E64B35B2", "#8491B4B2"))  +  # 设置点的颜色
  scale_fill_manual(values =  c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_88770$ki67),max(score_88770$ki67))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

#-------转移----------
##--------20685--------
library(survival)
survifit <- survfit(Surv(OS.time,OS)~metastasis,score_20685)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~metastasis,score_20685)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
###-------------------K-M曲线------------------
pdf("KM_metastasis20685.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data = score_20685, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "metastasis",
  legend.labs = c("0", "1"),
  xlab = "Time (years)",
  break.time.by =3.5,
  palette = c( "#8491B4AA","#E64B35AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
score_20685$metastasis<-as.factor(score_20685$metastasis)
###-----boxplot--------
pdf("boxplot_metastasis20685.pdf", width = 6, height = 5.5)
ggplot(score_20685, aes(x = metastasis, y =Score, color =metastasis)) +
  geom_boxplot(aes(fill =metastasis), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2" )) +  # 设置点的颜色
  scale_fill_manual(values = c("#8491B4B2","#E64B35B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_20685$Score),max(score_20685$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

##-----48390-----
gset <- getGEO('GSE48390', destdir="D:/RStudio/Ful数据/F3/测试集/GSE48390",
               AnnotGPL = T,     
               getGPL = T) 
pdata<-pData(gset[[1]])
rr<-pdata[,40,drop=FALSE]
colnames(rr)<-'event'
rr$event <- ifelse(rr$event == "disease-free", 0, 1)
score_48390<-cbind(score_48390,rr)
score_48390$event<-as.factor(score_48390$event)
survifit <- survfit(Surv(OS.time,OS)~event,score_48390)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~event,score_48390)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
###-------------------K-M曲线------------------
pdf("KM_event48390.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data = score_48390, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "metastasis or recurrence",
  legend.labs = c("0", "1"),
  xlab = "Time (years)",
  break.time.by =1.5,
  palette = c( "#8491B4AA","#E64B35AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()

###-----boxplot--------
pdf("boxplot_metastasis48390.pdf", width = 6, height = 5.5)
ggplot(score_48390, aes(x = event, y =Score, color =event)) +
  geom_boxplot(aes(fill =event), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#8491B4B2","#E64B35B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_48390$Score),max(score_48390$Score))) +
  stat_compare_means(method = "t.test")  # 添加两组之间的p值
dev.off()
##------58812--------
gset <- getGEO('GSE58812', destdir="D:/RStudio/Ful数据/F3/测试集/GSE58812",
               AnnotGPL = T,     
               getGPL = T) 
pdata<-pData(gset[[1]])
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<- data.frame(sapply(gpl$`Gene symbol`, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
#代入模型计算评分
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))
gse58812<-exprset_unique[comsam,,drop=F]
gse58812<-t(gse58812)
exprSet <- log2(gse58812+1)
exprSet <-scale(exprSet )
GSE58812os$OS<-as.numeric(GSE58812os$OS)
gse58812<-cbind(GSE58812os,exprSet)
rsf_v<- predict(fit,newdata = gse58812,proximity=T)
score_58812<- data.frame(gse58812[,c(1,2)],Score=rsf_v$predicted)
meta<-pdata[,44,drop=F]
colnames(meta)<-'metastasis'
score_58812<-cbind(meta,score_58812)
score_58812$metastasis<-as.factor(score_58812$metastasis)
survifit <- survfit(Surv(OS.time,OS)~metastasis,score_58812)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~metastasis,score_58812)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
###-------------------K-M曲线------------------
pdf("KM_metastasis58812.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data = score_58812, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "metastasis",
  legend.labs = c("0", "1"),
  xlab = "Time (years)",
  break.time.by =1.5,
  palette = c( "#8491B4AA","#E64B35AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
score_58812$metastasis<-as.factor(score_58812$metastasis)
###-----boxplot--------
pdf("boxplot_metastasis58812.pdf", width = 6, height = 5.5)
ggplot(score_58812, aes(x = metastasis, y =Score, color =metastasis)) +
  geom_boxplot(aes(fill =metastasis), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2" )) +  # 设置点的颜色
  scale_fill_manual(values = c("#8491B4B2","#E64B35B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_58812$Score),max(score_58812$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()
##---------GSE10886------
gset <- getGEO('GSE10886', destdir="D:/RStudio/Ful数据/训练6/F6/GSE10086",
               AnnotGPL = T,     
               getGPL = T) 
pdata<-pData(gset[[1]])
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<- data.frame(sapply(gpl$`Gene symbol`, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
#代入模型计算评分
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))



#-----lymph node-----
##----37751-----
gset <- getGEO('GSE37751', destdir="D:/RStudio/Ful数据/训练6/F6/GSE37751",
               AnnotGPL = T,     
               getGPL = T) 
pdata<-pData(gset[[1]])
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<- data.frame(sapply(gpl$`Gene symbol`, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
#代入模型计算评分
comsam<-intersect(rownames(exprset_unique),rownames(importance_gene))
gse37751<-exprset_unique[comsam,]
GSE37751os$OS<-as.numeric(GSE37751os$OS)
gse37751<-t(gse37751)
gse37751<-scale(gse37751)
comsam<-intersect(rownames(gse37751),rownames(GSE37751os))
gse37751<-gse37751[comsam,]
gse37751<-cbind(GSE37751os,gse37751)
gse37751$OS.time<-gse37751$OS.time/365
rsf_v<- predict(fit,newdata = gse37751,proximity=T)
score_37751<- data.frame(gse37751[,c(1,2)],Score=rsf_v$predicted)
node<-pdata[,51,drop=F]
colnames(node)<-'lymph_node'
node<-node[comsam,]
score_37751<-cbind(node,score_37751)
score_37751$node<-as.factor(score_37751$node)
library(ggplot2)
library(ggalluvial)
save(score_37751,file="score_37751.RData")
cut <- surv_cutpoint(score_37751,'OS.time','OS','Score')
summary(cut)#106.3041
score_37751$group <- ifelse(score_37751$Score >=106.3041, "high", "low")
###---plot---
# 设置图形主题
plot_theme <- theme_minimal() +
  theme(legend.position = "none",  # 隐藏图例
        text = element_text(family = "sans",
                            face = "bold",  # 设置字体为粗体
                            size = 14),  # 设置字体大小
        panel.grid = element_blank())  # 去除背景网格线
# 创建数据可视化图
score_37751$frequency <- with(score_37751, ave(seq_along(node), node, group, FUN = length))
gg <- ggplot(score_37751, aes(y = frequency, axis1 =node, axis2=group)) +
  # geom_alluvium 用于绘制中间的曲线，包括颜色填充、曲线类型和线宽设置
  geom_alluvium(aes(fill = node), curve_type = "spline", size = 2, width = 0.2) +
  # 设置ID列框的格式和大小
  geom_stratum(aes(fill = node),  width = 0.2) +
  # 设置另外两列框格式和大小
  geom_stratum(aes(fill = group), width = 0.2) +
  # 使用颜色板 "Set3" 进行填充，也可自定义颜色版
  scale_fill_brewer(palette = "Set3") +
  # 设置框上的标签
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3,
            angle = 0) +
  # 自定义X和Y轴上的标签
  scale_x_discrete(labels = c("node", "group")) +
  scale_y_continuous(expand = c(0, 0)) +
  # 定义X和Y轴的名称
  labs(x = "Class", y = "frequency")
# 输出图形到屏幕
print(gg + plot_theme)

##--42568-----
score_42568$`lymph node status:ch1`<-as.factor(score_42568$`lymph node status:ch1`)
colnames(score_42568)[2]<-'node'
p <- ggplot(score_42568, aes(x = node, y =Score, color = node)) +
  geom_boxplot(aes(fill =node), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_42568$Score),max(score_42568$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
p


