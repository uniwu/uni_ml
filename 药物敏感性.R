
install.packages("oncoPredict")
BiocManager::install('GenomicFeatures')
BiocManager::install('TCGAbiolinks')
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(oncoPredict)

library(tidyverse)
library(cowplot)
library(rstatix)
dir <- "D:/RStudio/Ful数据/耐药/DataFiles/DataFiles/DataFiles/Training Data"
#训练集表达谱数据
GDSC2_Expr = readRDS(file = file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
#IC50数据
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)

#读取自己的表达矩阵,TPM格式
setwd('D:/RStudio/Ful数据/训练6/F7')
# 按照type列排序
score_tcga <- score_tcga %>%
  arrange(type)
exp<-t(tpms)
exp <- exp[rownames(score_tcga), ]
table(score_tcga$type)#534
group <- rownames(exp) %>% as.data.frame()
group$group <- c(rep('high',534),rep('low',534))
group$Score<-score_tcga$Score
colnames(group)[1] <- 'ID'
#log化
dataexp2 <- log(exp+1)
dataexp2<-t(dataexp2)
oncoPredict::calcPhenotype(trainingExprData = GDSC2_Expr,
                           trainingPtype = GDSC2_Res,
                           testExprData = dataexp2,
                           batchCorrect = 'eb',  #   "eb" for ComBat  
                           powerTransformPhenotype = TRUE,
                           removeLowVaryingGenes = 0.2,
                           minNumSamples = 10, 
                           printOutput = TRUE, 
                           removeLowVaringGenesFrom = 'rawData' )

#读取结果以及合并组别
resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', header = T ,stringsAsFactors = F ,check.names = F)
names(resultPtype)[1] <- "ID"
drugs <- group %>%
  inner_join(resultPtype,by='ID')
rownames(importance_gene)<-importance_gene$gene
comsam<-intersect(rownames(importance_gene),rownames(dataexp2))
dataexp2<-dataexp2[comsam,,drop=F]
exp3<-t(dataexp2)
rownames(drugs)<-drugs$ID
write(drugs,'drugs.xlsx')
#计算相关性
##----score----
results <- lapply(ic50_data, function(column) {
  cor.test(drugs$Score, column)
})
correlation_results <- data.frame(
  Column = names(ic50_data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
correlation_results$x<-'Score'
#筛选p
filteredcor<- correlation_results [correlation_results $P_Value < 0.05, ]
#按照相关系数的绝对值排序
sortedcor<- filteredcor[order(abs(filteredcor$Correlation), decreasing = TRUE), ]
save(sortedcor,file='Scorecor.RData')
#筛选p
exp3cor_1<- exp3cor [exp3cor $P_Value < 0.05, ]
#筛选相关性前50
exp3cor_1<- exp3cor_1[order(abs(exp3cor_1$Correlation), decreasing = TRUE), ]
exp3cor_2<- exp3cor_1[1:50,]
library(tidyverse)
sortedcor<- filteredcor[order(abs(filteredcor$Correlation), decreasing = TRUE), ]
score50<-sortedcor[1:50,]
rownames(score50) <- gsub("\\.cor", "", rownames(score50))
comsam<-intersect(rownames(score50),colnames(ic50_data))
ic50data<-ic50_data[,comsam]

##----------LAMA3--------
results <- lapply(ic50data, function(column) {
  cor.test(exp3$LAMA3, column)
})
LAMA3 <- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
LAMA3$x<-'LAMA3'
LAMA3$drug<-rownames(LAMA3)
##----------morn3------
results <- lapply(ic50data, function(column) {
  cor.test(exp3$MORN3, column)
})
MORN3<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
MORN3$x<-'MORN3'
MORN3$drug<-rownames(MORN3)
##-----il27ra-------
results <- lapply(ic50data, function(column) {
  cor.test(exp3$IL27RA, column)
})
IL27RA<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
IL27RA$x<-'IL27RA'
IL27RA$drug<-rownames(IL27RA)
##----IGSF10----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$IGSF10, column)
})
IGSF10<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
IGSF10$x<-'IGSF10'
IGSF10$drug<-rownames(IGSF10)
##---PXK----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$PXK, column)
})
PXK<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
PXK$x<-'PXK'
PXK$drug<-rownames(PXK)
##----CIR1----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$CIR1, column)
})
CIR1<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
CIR1$x<-'CIR1'
CIR1$drug<-rownames(CIR1)
##----PRAME----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$PRAME, column)
})
PRAME<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
PRAME$x<-'PRAME'
PRAME$drug<-rownames(PRAME)
##---SIX2---
results <- lapply(ic50data, function(column) {
  cor.test(exp3$SIX2, column)
})
SIX2<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
SIX2$x<-'SIX2'
SIX2$drug<-rownames(SIX2)
##----SPDL1----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$SPDL1, column)
})
SPDL1<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
SPDL1$x<-'SPDL1'
SPDL1$drug<-rownames(SPDL1)
##---ATG12----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$ATG12, column)
})
ATG12<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
ATG12$x<-'ATG12'
ATG12$drug<-rownames(ATG12)
##----AATF---
results <- lapply(ic50data, function(column) {
  cor.test(exp3$AATF, column)
})
AATF<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
AATF$x<-'AATF'
AATF$drug<-rownames(AATF)
##-----CD59----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$CD59, column)
})
CD59<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
CD59$x<-'CD59'
CD59$drug<-rownames(CD59)
##---IQUB----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$IQUB, column)
})
IQUB<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
IQUB$x<-'IQUB'
IQUB$drug<-rownames(IQUB)
##----TOMM34---
results <- lapply(ic50data, function(column) {
  cor.test(exp3$TOMM34, column)
})
TOMM34<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
TOMM34$x<-'TOMM34'
TOMM34$drug<-rownames(TOMM34)
##----ZNF148---
results <- lapply(ic50data, function(column) {
  cor.test(exp3$ZNF148, column)
})
ZNF148<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
ZNF148$x<-'ZNF148'
ZNF148$drug<-rownames(ZNF148)
##----VPS37A----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$VPS37A, column)
})
VPS37A<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
VPS37A$x<-'VPS37A'
VPS37A$drug<-rownames(VPS37A)
##------C12orf57-----
results <- lapply(ic50data, function(column) {
  cor.test(exp3$C12orf57, column)
})
C12orf57<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
C12orf57$x<-'C12orf57'
C12orf57$drug<-rownames(C12orf57)
##-----SPATA7---
results <- lapply(ic50data, function(column) {
  cor.test(exp3$SPATA7, column)
})
SPATA7<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
SPATA7$x<-'SPATA7'
SPATA7$drug<-rownames(SPATA7)
##-----EP400NL---
results <- lapply(ic50data, function(column) {
  cor.test(exp3$EP400NL, column)
})
EP400NL<- data.frame(
  Column = names(ic50data),
  Correlation = sapply(results, function(x) x$estimate),
  P_Value = sapply(results, function(x) x$p.value)
)
EP400NL$x<-'EP400NL'
EP400NL$drug<-rownames(EP400NL)
exp3cor<- rbind(LAMA3[, c(2:5)], MORN3[, c(2:5)], IL27RA[, c(2:5)], 
                IGSF10[, c(2:5)], PXK[, c(2:5)], CIR1[, c(2:5)], 
                PRAME[, c(2:5)], SIX2[, c(2:5)], SPDL1[, c(2:5)], 
                ATG12[, c(2:5)], AATF[, c(2:5)], CD59[, c(2:5)], 
                IQUB[, c(2:5)], TOMM34[, c(2:5)], ZNF148[, c(2:5)], 
                VPS37A[, c(2:5)], C12orf57[, c(2:5)], SPATA7[, c(2:5)], 
                EP400NL[, c(2:5)])
score50$x<-'Score'
score50$drug<-rownames(score50)
exp3cor_3<-rbind(exp3cor, score50[,c(2:5)])
exp3cor_3$drug<-gsub("\\.cor", "", exp3cor_3$drug)
write.xlsx(exp3cor_3, "exp3cor_3.xlsx")

#根据正负相关划分p值
exp3cor_3$fdr_group <- cut(exp3cor_3$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
exp3cor_3$fdr_group <- factor(exp3cor_3$fdr_group, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
exp3cor_3$fdr_group_2 <- cut(exp3cor_3$P_Value, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05"))
exp3cor_3$fdr_group_2 <- factor(exp3cor_3$fdr_group_2, levels = c(labels = c("<0.0001","<0.001", "<0.01", "<0.05", ">0.05")))
x <- c('AATF', 'ATG12', 'C12orf57', 'CD59', 'CIR1', 
       'EP400NL', 'IGSF10', 'IL27RA', 'IQUB', 'LAMA3', 
       'MORN3', 'PRAME', 'PXK', 'Score', 'SIX2', 
       'SPATA7', 'SPDL1', 'TOMM34', 'VPS37A', 'ZNF148')
##-plot------
ggplot2::ggplot() +
  geom_hline(yintercept = x, color = "#E8E8E8") +
  geom_point(data = subset(exp3cor_3, Correlation > 0), shape = 19, stroke = 0, # 正相关
             aes(x = drug, y = x, size = abs(Correlation), color = fdr_group)) +
  geom_point(data = subset(exp3cor_3, Correlation < 0), shape = 21, stroke = 0.1, # 负相关
             aes(x = drug, y = x, size = abs(Correlation), fill = fdr_group_2)) +
  scale_size_continuous(limits = c(0, 0.6), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
  scale_fill_manual(values = c("#8491B4", "#E8E8E8", "#B0C4DE", "#778899", "#ADD8E6")) + # 自定义填充颜色
  scale_color_manual(values = c("#E64B35", "#E8E8E8", "#FFA07A", "#FF6347", "#FA8072")) + # 自定义边框颜色
  cowplot::theme_cowplot() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(color = "Positive\ncorrelation\nP_Value", size = "Correaltion",
       fill = "Negative\ncorrelation\nP_Value", x = "", y = "",
       title = "Correlation of drug sensitive") +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 3),
         color = guide_legend(override.aes = list(size = 4), order = 2),
         size = guide_legend(override.aes = list(size = c(2:7), fill = "white"), order = 1))


#------------boxplot------------
##------------Camptothecin--------
pdf("boxplot_Camptothecin_ic50.pdf", width = 6, height = 5.5)
ggplot(drugs, aes(x = group, y =Camptothecin_1003, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(drugs$Camptothecin_1003),max(drugs$Camptothecin_1003))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

##----------Epirubicin------------
pdf("boxplot_Epirubicin_ic50.pdf", width = 6, height = 5.5)
ggplot(drugs, aes(x = group, y =Epirubicin_1511, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(drugs$Epirubicin_1511),max(drugs$Epirubicin_1511))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

##------Mitoxantrone----------
pdf("boxplot_Mitoxantroneic50.pdf", width = 6, height = 5.5)
ggplot(drugs, aes(x = group, y =Mitoxantrone_1810, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(drugs$Mitoxantrone_1810),max(drugs$Mitoxantrone_1810))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

##-------Palbociclib------
pdf("boxplot_Palbociclib_ic50.pdf", width = 6, height = 5.5)
ggplot(drugs, aes(x = group, y =Palbociclib_1054, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(drugs$Palbociclib_1054),max(drugs$Palbociclib_1054))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

##-------Trametinib--------
pdf("boxplot_Trametinib_1372_ic50.pdf", width = 6, height = 5.5)
ggplot(drugs, aes(x = group, y =Trametinib_1372, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(drugs$Trametinib_1372),max(drugs$Trametinib_1372))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()

#--------散点图------
library(ggplot2)
library(ggpubr)
library(ggExtra)
##------------Camptothecin--------
x=as.numeric(drugs$Score)
y=as.numeric(drugs$Camptothecin_1003)
pdf("scatter_Camptothecin_ic50.pdf", width = 6, height = 5.5)
p1 <- ggplot(drugs, aes(x, y)) + 
  xlab('Score') + ylab('Camptothecin_1003') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
dev.off()
##----------Epirubicin------------
y=as.numeric(drugs$Epirubicin_1511)
pdf("scatter_Epirubicin_ic50.pdf", width = 6, height = 5.5)
p1 <- ggplot(drugs, aes(x, y)) + 
  xlab('Score') + ylab('Epirubicin_1511') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
dev.off()
##------Mitoxantrone----------
y=as.numeric(drugs$Mitoxantrone_1810)
pdf("scatter_Mitoxantroneic50.pdf", width = 6, height = 5.5)
p1 <- ggplot(drugs, aes(x, y)) + 
  xlab('Score') + ylab('Mitoxantrone_1810') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
dev.off()
##-------Palbociclib------
y=as.numeric(drugs$Palbociclib_1054)
pdf("scatter_Palbociclib_ic50.pdf", width = 6, height = 5.5)
p1 <- ggplot(drugs, aes(x, y)) + 
  xlab('Score') + ylab('Palbociclib_1054') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
dev.off()
##-------Trametinib-------
y=as.numeric(drugs$Trametinib_1372)
pdf("scatter_Trametinib_1372_ic50.pdf", width = 6, height = 5.5)
p1 <- ggplot(drugs, aes(x, y)) + 
  xlab('Score') + ylab('Trametinib_1372') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
dev.off()


#---------------TIDE------------
#未log化的tpm数据进行标准化
exprSet <- log2(tpms+1)
mydat <- t(apply(exprSet, 1, function(x)x-(mean(x))))
write.table(mydat,file = 'TIDE1.txt',sep = '\t') 
patient_ICBresponse <- read.csv("TIDEoutpu1.csv",header = TRUE)
rownames(patient_ICBresponse)<-patient_ICBresponse$Patient
ICBresponse <- patient_ICBresponse
##统计免疫响应患者数
table(ICBresponse$Responder=="False")
#无响应729 响应339
library(stringr)
ICBresponse$Responder <- ifelse(str_detect(ICBresponse$Responder,"False"),"NR","R")
ICBresponse1 <- dplyr::pull(ICBresponse, 3) ##将data.frame中的一列转换为向量
names(ICBresponse1) <- ICBresponse$Patient
ICBresponse1[1:10]
save(ICBresponse,file = "ICBresponse.Rdata")
colnames(ICBresponse)[1]<-'ID'
merged<-merge(ICBresponse,drugs,by="ID")
##-----------发散条形图----------
## barplot
dat_plot <- data.frame(id = ICBresponse$ID,
                       t = ICBresponse$TIDE,
                       threshold=ICBresponse$Responder)
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
library(ggprism)
ggplot(data = dat_plot, aes(x = id, y = t, fill = threshold)) +
  geom_col() +
  scale_fill_manual(values = c('R' = '#FF1493', 'NR' = '#40E0D0')) +
  geom_hline(yintercept = c(-2, 2), color = 'white', linewidth = 0.5, linetype = 'dashed') + # 使用 linewidth 替换 size
  xlab('') + 
  ylab('TIDE score') + 
  guides(fill = guide_legend(key.width = 3, key.height = 5, nrow = 2, ncol = 1, byrow = TRUE)) +
  theme_prism(border = TRUE) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position.inside = c(.95, .95) # 使用 legend.position.inside 替换 legend.position
  )

ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col(position = position_dodge2()) + 
  # scale_x_continuous(breaks = seq(0,50,100))+ # X 轴每隔 50 个单位显示一个刻度
  scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0')) +
  geom_hline(yintercept = c(-2,2),color = 'white',linewidth = 0.5,lty='dashed') +
  xlab('') + 
  ylab('TIDE score') + 
  guides(fill=guide_legend(key.width = 3, key.height = 5, nrow = 2, ncol = 1, byrow = TRUE))+ #显示图例
  theme_prism(border = T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position.inside = c(.95, .95),#图例位置
    legend.justification = c("right", "top")#固定右上角
  )
##-----------堆积条形图---------------

group <- c("high","high","low","low")
Response <- c("NR","R","NR","R")
table(merged$group[merged$Responder=="R"])
#high  low  193  146 
table(merged$group[merged$Responder=="NR"])
#high  low  341  388 
Percent <- c(0.57,0.43,0.47,0.53)
num <- c(57,43,47,53)
data <- data.frame(group,Response,Percent,num)
##添加统计数值(卡方检验)
NR <- c(57,43)
R <- c(47,53)
dat <- data.frame(NR,R)
rownames(dat) <- c("high","low")
chisq.test(dat)
#X-squared = 1.6226, df = 1, p-value = 0.2027
ggplot( data, aes( x = group, weight = Percent, fill = Response))+
  geom_bar( position = "stack")+xlab("X-squared = 1.6226, df = 1, p-value = 0.2027")+
  scale_fill_manual(values = c("R" = "#E64B35B2", "NR" = "#8491B4B2"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5))

##--------箱线图---------
pdf("boxplot_Camptothecin_ic50.pdf", width = 6, height = 5.5)
ggplot(merged, aes(x = group, y =TIDE, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(merged$TIDE),max(merged$TIDE))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()
x=as.numeric(merged$TIDE)
y=as.numeric(merged$Score)
p1 <- ggplot(merged, aes(x, y)) + 
  xlab('TIDE') + ylab('Score') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
###-----------MSI-----------
merged2<-merge(group,ICBresponse,by="ID")
pdf("highlowMSI.pdf", width = 6, height = 5.5)
ggviolin(merged2, "group", "MSI.Expr.Sig" ,
         fill = "group", #小提琴内部颜色对应的数据列
         color = "group", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T,
         palette = c("#8491B4B2","#E64B35B2"),  #小提琴自定义颜色
         legend = "right", #图例添加在图的右侧
         legend.title = " ",#图例的标题
         font.legend = c(10, "plain", "black"), #图例字体的大小/样式/颜色
         font.y = 10,  #y轴标题字体大小
         font.tickslab = c(10,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  labs(x = NULL, #设置x轴标题
       y = "MSI", #设置y轴标题
       title = NULL)+
  stat_compare_means(method = "t.test")
dev.off()
###---------Dysfunction-------
pdf("highlowDysfunction.pdf", width = 6, height = 5.5)
ggviolin(merged2, "group", "Dysfunction" ,
         fill = "group", #小提琴内部颜色对应的数据列
         color = "group", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T,
         palette = c("#8491B4B2","#E64B35B2"),  #小提琴自定义颜色
         legend = "right", #图例添加在图的右侧
         legend.title = " ",#图例的标题
         font.legend = c(10, "plain", "black"), #图例字体的大小/样式/颜色
         font.y = 10,  #y轴标题字体大小
         font.tickslab = c(10,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  labs(x = NULL, #设置x轴标题
       y = "Dysfunction", #设置y轴标题
       title = NULL)+
  stat_compare_means(method = "t.test")
dev.off()

###----------Exclusion----------
pdf("highlowExclusion.pdf", width = 6, height = 5.5)
ggviolin(merged2, "group", "Exclusion" ,
         fill = "group", #小提琴内部颜色对应的数据列
         color = "group", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T,
         palette = c("#8491B4B2","#E64B35B2"),  #小提琴自定义颜色
         legend = "right", #图例添加在图的右侧
         legend.title = " ",#图例的标题
         font.legend = c(10, "plain", "black"), #图例字体的大小/样式/颜色
         font.y = 10,  #y轴标题字体大小
         font.tickslab = c(10,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  labs(x = NULL, #设置x轴标题
       y = "Exclusion", #设置y轴标题
       title = NULL)+
  stat_compare_means(method = "t.test")
dev.off()


#-----------Easier包预测免疫治疗响应,无TMB-------
BiocManager::install("easier")
BiocManager::install("easier", dependencies = TRUE)
# 下载开发版本
BiocManager::install("olapuentesantana/easier")
library("easier")
##------------计算免疫反应的特征---------
#compute_scores_immune_response对多个不同的免疫特征进行打分
immune_response_scores <- compute_scores_immune_response(RNA_tpm =tpms)
save(immune_response_scores,file = "immune_response_scores.RData")
##----------对TME定量描述的计算-----------
###----------细胞成分---------
cell_fractions <- compute_cell_fractions(RNA_tpm = tpms)
save(cell_fractions,file = "cell_fractions.RData")
###-------通路活性方面-----------
pathway_activities <- compute_pathway_activity(RNA_counts = counts, 
                                               remove_sig_genes_immune_response = TRUE)
save(pathway_activities,file = "pathway_activities.RData")
###---------------转录因子活性----------------
tf_activities <- compute_TF_activity(RNA_tpm = tpms)
save(tf_activities,file = "tf_activities.RData")
###----配体-受体权重-------
lrpair_weights <- compute_LR_pairs(RNA_tpm = tpms,cancer_type = "BRCA")
save(lrpair_weights,file = "lrpair_weights.RData")
###--------------量化细胞间的关联强度------------
ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, 
                                  cancer_type = "pancan")
save(ccpair_scores,file = "ccpair_scores.RData")
###----获取患者对免疫反应的预测--------
predictions <- predict_immune_response(pathways = pathway_activities,
                                       immunecells = cell_fractions,
                                       tfs = tf_activities,
                                       lrpairs = lrpair_weights,
                                       ccpairs = ccpair_scores,
                                       cancer_type = 'BRCA', 
                                       verbose = TRUE)
save(predictions,file = "predictions.RData")
easier_derived_scores <- retrieve_easier_score(predictions_immune_response = predictions,
                                               weight_penalty = 0.5)
easier_derived_scores$ID<-rownames(easier_derived_scores)
easier_derived_scores<-merge(easier_derived_scores,group,by="ID")
###---------------plot---------------
pdf("无TMBeasierscore.pdf", width = 6, height = 5.5)
ggplot(easier_derived_scores, aes(x = group, y =easier_score, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2","#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2","#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(easier_derived_scores$easier_score),max(easier_derived_scores$easier_score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()
x=as.numeric(easier_derived_scores$Score)
y=as.numeric(easier_derived_scores$easier_score)
pdf("scatter_easierscore.pdf", width = 6, height = 5.5)
p1 <- ggplot(easier_derived_scores, aes(x, y)) + 
  xlab('Score') + ylab('easier_score') +
  geom_point(color = "#8491B4B2") + # 自定义散点颜色
  geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = x, y = y))
ggMarginal(p1, type = "density", xparams = list(fill = "#E64B35B2"), yparams = list(fill = "#8491B4B2"))
dev.off()

#-----------Easier包预测免疫治疗响应,有TMB-------
##-----------肿瘤突变负荷TMB------------
library("easier")
##------------计算免疫反应的特征---------
#compute_scores_immune_response对多个不同的免疫特征进行打分
immune_response_scores1 <- immune_response_scores[comsam,]
##----------对TME定量描述的计算-----------
###----------细胞成分---------
cell_fractions1 <- cell_fractions[comsam,]
###-------通路活性方面-----------
pathway_activities1 <- pathway_activities[comsam,]

###---------------转录因子活性----------------
tf_activities1 <- tf_activities[comsam,]
###----配体-受体权重-------
lrpair_weights1 <- lrpair_weights[comsam,]
###--------------量化细胞间的关联强度------------
ccpair_scores1 <- ccpair_scores[comsam,]
###----获取患者对免疫反应的预测--------
predictions1 <- predict_immune_response(pathways = pathway_activities1,
                                       immunecells = cell_fractions1,
                                       tfs = tf_activities1,
                                       lrpairs = lrpair_weights1,
                                       ccpairs = ccpair_scores1,
                                       cancer_type = 'BRCA', 
                                       verbose = TRUE)
save(predictions1,file = "predictions1.RData")
#TMB
library(TCGAmutations)
BRCA <- TCGAmutations::tcga_load(study = "BRCA")#source = "MC3"
maf = tmb(maf = BRCA,
          captureSize = 50,
          logScale = TRUE)   #若logScale=F，则不会出现log
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
maf$sample <- gsub("-",".",maf$sample)
maf$sample <- substr(maf$sample, 1, nchar(maf$sample) - 4)
#重复行取第一次出现
maf<- maf[!duplicated(maf$sample), ]
maf<-as.data.frame(maf)
tpms1<-t(tpms)
tpms1<-as.data.frame(tpms1)
tpms1$sample<-rownames(tpms1)
maf1<-merge(maf,tpms1,by="sample")
TMB<-maf1$total_perMB_log
TMB<-na.omit(TMB)
easier_derived_scores1 <- retrieve_easier_score(predictions_immune_response = predictions1,
                                               TMB_values = TMB,
                                               easier_with_TMB = c("weighted_average", 
                                                                   "penalized_score"),
                                               weight_penalty = 0.5)

###------plot------
#高低评分的TMB差异
pdf("TMB.pdf", width = 6, height = 5.5)
group$sample<-group$ID
merged1<-merge(maf1,group,by="sample")
#剔除离群值
#可视化离群值
boxplot(merged1$total_perMB) #画箱式图
q1 <- quantile(merged1$total_perMB, 0.25)
q3 <- quantile(merged1$total_perMB, 0.75)
iqr <- q3 - q1
# 计算概述统计量
summary(merged1$total_perMB)
# 使用IQR方法检测离群值
low <- q1 - 1.5 * iqr
high <- q3 + 1.5 * iqr
# 从数据集中删除离群值
merged1 <- merged1[merged1$total_perMB >= low & merged1$total_perMB <= high, ]
pdf("highlowTMB.pdf", width = 6, height = 5.5)
ggviolin(merged1, "group", "total_perMB" ,
         fill = "group", #小提琴内部颜色对应的数据列
         color = "group", #小提琴边框颜色对应的数据列  有时候这里会用默认即black，
         trim = T,
         palette = c("#8491B4B2","#E64B35B2"),  #小提琴自定义颜色
         legend = "right", #图例添加在图的右侧
         legend.title = " ",#图例的标题
         font.legend = c(10, "plain", "black"), #图例字体的大小/样式/颜色
         font.y = 10,  #y轴标题字体大小
         font.tickslab = c(10,"plain","black"), #x轴 y轴刻度大小/样式/颜色
         add = "boxplot",  #叠加箱线图
         add.params = list(
           fill = "white", #设置箱线图内部颜色
           color = "black",  #设置箱线图边框
           width = 0.2,   #箱线图的宽度
           linetype = 1)) +
  labs(x = NULL, #设置x轴标题
       y = "TMB", #设置y轴标题
       title = NULL)+
  stat_compare_means(method = "t.test")
dev.off()




#--------免疫相关热图------
immune_response_scores$ID<-rownames(immune_response_scores)
imm<-merge(immune_response_scores,group,by="ID")
library(corrplot)  
cor <- cor(imm)  # 计算相关性矩阵
res <- cor.mtest(imm, conf.level = 0.95)
p <- res$p
mycol <- colorRampPalette(c("#8491B4B2","white" ,"#E64B35B2"),alpha = TRUE)  # 自定义颜色向量mycol，用于热图颜色映射
pdf("immheatmap.pdf", width = 6, height = 5.5)
# 绘制上半部分的相关性矩阵，使用饼图表示相关性，同时添加显著性水平标记：
corrplot(
  cor,                    # 相关性矩阵
  method = c('pie'),      # 绘图方法，使用饼图表示相关性
  type = c('upper'),      # 绘制的部分，这里绘制矩阵的上半部分
  col = mycol(100),       # 颜色设置
  outline = 'grey',       # 边框颜色
  order = c('AOE'),       # 相关性矩阵元素的排列顺序
  diag = TRUE,            # 包括对角线上的相关性
  tl.cex = 0.4,             # 刻度标签文本大小
  tl.col = 'black',       # 刻度标签文本颜色
  tl.pos = 'd',           # 刻度标签位置为下方
  p.mat = p,              # 显著性水平矩阵
  sig.level = c(0.001, 0.01, 0.05),  # 显著性水平阈值
  insig = "label_sig",    # 显著性不满足时的标注样式，这里使用标签 "label_sig"
  pch.cex = 0.8,          # 显著性标记的大小
  pch.col = 'black'       # 显著性标记的颜色
)
# 绘制下半部分的相关性矩阵，不显著的相关性使用"X"标记：
corrplot(
  cor,                    # 相关性矩阵
  add = TRUE,             # 将这次的绘图结果添加到之前的图中
  method = c('number'),   # 绘图方法，使用数字表示相关性
  type = c('lower'),      # 绘制的部分，这里绘制矩阵的下半部分
  col = mycol(100),       # 颜色设置
  order = c('AOE'),       # 相关性矩阵元素的排列顺序
  diag = FALSE,           # 不包括对角线上的相关性
  number.cex = 0.6,       # 相关性系数数字标签大小
  tl.pos = 'n',           # 不显示刻度标签
  cl.pos = 'n',           # 不显示颜色条
  p.mat = p,              # 显著性水平矩阵
  insig = "pch"           # 不显著的相关性使用 "X" 标记
)
dev.off()


#drug1
drug1<-t(drugs)
drug1<-as.data.frame(drug1)
drug1$drug<-rownames(drug1)
#Tamoxifen_1199
ggplot(drugs, aes(x = group, y = Tamoxifen_1199, color = group)) +
  geom_boxplot(aes(fill = group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  
  scale_color_manual(values = c("#E64B35","#83B8D7")) +  
  scale_fill_manual(values = c("#E64B35","#83B8D7")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("high", "low")), 
              map_signif_level = TRUE, test='t.test',
              textsize = 4, vjust = 0.5,
              color = "black") +
  labs(x = "Tamoxifen", y = "Estimated ic50")

#Palbociclib_1054
ggplot(drugs, aes(x = group, y =Palbociclib_1054, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35","#83B8D7")) +  
  scale_fill_manual(values = c("#E64B35","#83B8D7")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("high", "low")), 
              map_signif_level = TRUE, test='t.test',
              textsize = 4, vjust = 0.5,
              color = "black") +
  labs(x = "Palbociclib", y = "Estimated ic50")

#Olaparib_1017
ggplot(drugs, aes(x = group, y =Olaparib_1017, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35","#83B8D7")) +  
  scale_fill_manual(values = c("#E64B35","#83B8D7")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("high", "low")), 
              map_signif_level = TRUE, test='t.test',
              textsize = 4, vjust = 0.5,
              color = "black") +
  labs(x = "Olaparib", y = "Estimated ic50")

#Talazoparib_1259
ggplot(drugs, aes(x = group, y =Talazoparib_1259, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35","#83B8D7")) +  
  scale_fill_manual(values = c("#E64B35","#83B8D7")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("high", "low")), 
              map_signif_level = TRUE, test='t.test',
              textsize = 4, vjust = 0.5,
              color = "black") +
  labs(x = "Talazoparib", y = "Estimated ic50")

#	Alpelisib_1560
ggplot(drugs, aes(x = group, y =Alpelisib_1560, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35","#83B8D7")) +  
  scale_fill_manual(values = c("#E64B35","#83B8D7")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("high", "low")), 
              map_signif_level = TRUE, test='t.test',
              textsize = 4, vjust = 0.5,
              color = "black") +
  labs(x = "Alpelisib", y = "Estimated ic50")


#Gefitinib_1010
ggplot(drugs, aes(x = group, y =Gefitinib_1010, color =group)) +
  geom_boxplot(aes(fill =group), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35","#83B8D7")) +  
  scale_fill_manual(values = c("#E64B35","#83B8D7")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("high", "low")), 
              map_signif_level = TRUE, test='t.test',
              textsize = 4, vjust = 0.5,
              color = "black") +
  labs(x = "Gefitinib", y = "Estimated ic50")
