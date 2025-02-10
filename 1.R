setwd('D:/ful数据/训练6')
#-----------------------------合并TCGA和METABRIC数据------------------------------
#删除均值为0的行
# 仅保留在一半以上样本里表达的基因
TCGAtumor = TCGAtumor [apply(TCGAtumor , 1, function(x) sum(x > 0) > 0.5*ncol(TCGAtumor )), ]
library(sva)
batch<-c(rep('TCGA',1069),rep('METABRIC',1979))#批次信息
comsam<-intersect(rownames(TCGAtumor),rownames(meta))#16681
meta<-meta[comsam,]
tcga<-TCGAtumor[comsam,]
tm<-cbind(tcga,meta)
tm_numeric <- apply(tm, 2, as.numeric)
adjustedtm <- ComBat(dat = tm_numeric, batch = batch)
rownames(adjustedtm)<-rownames(tm)
save(adjustedtm,file = '去批次的tcga+meta.RData')



#--------------------------免疫浸润--------------------------------
##---------IOBR包计算各种免疫浸润细胞类型评分---------------------

#-------------------------------immunedeconv包---------------------------
library(immunedeconv)
cancer="BRCA"#TCGA癌症缩写
brca_vector <- rep(cancer, ncol(adjustedtm))
res_estimate <- deconvolute(adjustedtm, method="estimate",tumor = T,indications = brca_vector)
##------------------timer
res_timer <- deconvolute(adjustedtm, method="timer",tumor = T,indications = brca_vector)
write.table(res_timer,"timer.txt",sep="\t",row.names = F)
res_timer$cell_type=paste0(res_timer$cell_type,"_TIMER")
##-----------abis
anyNA(adjustedtm)
any(is.nan(adjustedtm))
any(is.infinite(adjustedtm))
#用行的中位数替换NA
adjustedtm <- t(apply(adjustedtm, 1, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)))
save(adjustedtm,file = '替换掉NA的tcga+metabric合并数据.RData')
res_abis <- deconvolute(adjustedtm, method="abis",tumor = T,indications = brca_vector)
write.table(res_abis,"abis.txt",sep="\t",row.names = F)
res_abis$cell_type=paste0(res_abis$cell_type,"_ABIS")
##------consensus_tme
res_consensus_tme <- deconvolute(adjustedtm, method="consensus_tme",tumor = T,indications = brca_vector)
write.table(res_consensus_tme,"consensus_tme.txt",sep="\t",row.names = F)
res_consensus_tme$cell_type=paste0(res_consensus_tme$cell_type,"_ConsensusTME")

##-------xcell
library(xCell)
res_xcell <- xCellAnalysis(adjustedtm,rnaseq=F)
write.table(res_xcell,"xcell.txt",sep="\t",row.names = F)
rownames(res_xcell)=paste0(rownames(res_xcell),"_xCell")
#Calling gsva(expr=., gset.idx.list=., method=., ...) is deprecated; use a method-specific parameter object (see '?gsva'). 
#In MulticoreParam(progressbar = verbose, workers = parallel.sz,
##------epic
res_epic <- deconvolute(adjustedtm, method="epic",tumor = T,indications = brca_vector)
write.table(res_epic,"epic.txt",sep="\t",row.names = F)
res_epic$cell_type=paste0(res_epic$cell_type,"_EPIC") 
#Warning message:In (function (bulk, reference = NULL, mRNA_cell = NULL, mRNA_cell_sub = NULL
     
##------quantiseq
res_quantiseq <- deconvolute(adjustedtm, method="quantiseq",tumor = T,indications = brca_vector) 
write.table(res_quantiseq,"quantiseq.txt",sep="\t",row.names = F)
res_quantiseq$cell_type=paste0(res_quantiseq$cell_type,"_quanTIseq")
#Signature genes found in data set: 129/138 (93.48%)

##---------cibersort
#采用cibersort方法时，比较麻烦，
#需要自己准备CIBERSORT.R和LM22.txt两个文件
#所以直接用CIBERSORT R包来做
library(devtools)
devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
#读取包自带的LM22文件（免疫细胞特征基因文件）
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#perm置换次数=1000，QN分位数归一化=TRUE
res_ciber <- cibersort(sig_matrix, adjustedtm,perm = 0,QN = T)
#perm：表示置换次数，数字越大运行时间越长，一般文章都设置为1000；
#QN：如果为芯片数据这里设为“T”；如果为测序数据设为“F”
res_ciber =res_ciber [,1:22] #最后3列不需要
write.table(res_ciber,"cibersort.txt",sep="\t",row.names = F)
#改成与其他算法格式相同的数据
res_ciber<-t(res_ciber)
rownames(res_ciber)=paste0(rownames(res_ciber),"_CIBERSORT")
res_ciber=data.frame(cell_type=rownames(res_ciber),res_ciber)
#非免疫细胞数据剔除
res_consensus_tme=res_consensus_tme[1:(nrow(res_consensus_tme)-1),]
res_xcell=res_xcell[1:(nrow(res_xcell)-3),]
file_names <- c("res_xcell.RData", "res_abis.RData", "res_ciber.RData", "res_consensus_tme.RData", "res_epic.RData", "res_quantiseq.RData", "res_timer.RData")
# 使用循环批量保存数据
for (i in seq_along(file_names)) {
  data_name <- paste0("res_", c("xcell", "abis", "ciber", "consensus_tme", "epic", "quantiseq", "timer")[i])
  save(list = data_name, file = file_names[i])
}
##----------------------免疫浸润淋巴细胞的筛选 ------------------------------------
###epic
rowepic<-c(1,3,4,7)
lepic<-res_epic[rowepic,]
###quantiseq
rowquan<-c(1,6,7,8,9)
lquan<-res_quantiseq[rowquan,]
###cibersort
rowciber<-c(1,2,4,5,6,7,8,9,10,11,12)
lciber<-res_ciber[rowciber,]
###xcell
# 定义要筛选的行名
res_xcell<-as.data.frame(res_xcell)
selected_row_names <- c("CD4+ T-cells_xCell", "CD4+ memory T-cells_xCell", "CD4+ naive T-cells_xCell", 
                        "CD4+ Tcm_xCell", "CD4+ Tem_xCell", "CD8+ T-cells_xCell", "CD8+ naive T-cells_xCell", 
                        "CD8+ Tem_xCell", "CD8+ Tcm_xCell", "Tregs_xCell", "Th1 cells_xCell", "Th2 cells_xCell", 
                        "B-cells_xCell", "pro B-cells_xCell", "naive B-cells_xCell", "Memory B-cells_xCell", 
                        "Class-switched memory B-cells_xCell", "NK cells_xCell", "NKT_xCell")
lxcell<- res_xcell[row.names(res_xcell) %in% selected_row_names, ]
###abis
rowabis<-c(2,3,4,5,6,7,8,9,11,14,17)
labis<-res_abis[rowabis,]
###timer
rowtimer<-c(1,2,3)
ltimer<-res_timer[rowtimer,]
###consensus_tme
rowconsensus<-c(1,7,9,10,11,12)
lconsensus<-res_consensus_tme[rowconsensus,]
save(labis,file = '淋巴abis.RData')
save(lconsensus,file = '淋巴consensus.RData')
save(lciber,file = '淋巴ciber.RData')
save(lepic,file = '淋巴epic.RData')
save(lquan,file = '淋巴quan.RData')
save(ltimer,file = '淋巴timer.RData')
save(lxcell,file = '淋巴xcell.RData')

##------------------相关系数计算----------------------------------------
###abis
rownames(labis)<-labis$cell_type
labis<-labis[,-1]
col_sums <- colSums(labis)
labis<- rbind(labis, col_sums)
rownames(labis)[12]<-'score'
#从表达矩阵中提取marker基因
marker<- subset(adjustedtm, rownames(adjustedtm) %in% c("CD3D", "CD4", "CD8A", "FOXP3", "PDCD1"))
marker<-marker[,(1:1068)]
###判断数据是否符合正态分布,这里只用了EPIC评分查看
library(tidyverse)
library(mice)
m<-t(marker)
le<-t(lepic)
le<-as.data.frame(le[,5])
check<-cbind(m,le)
# 提取基因表达值和评分列
gene_expression <- check[, 1:5] # 基因表达值的前五列
score <- check[, ncol(check)]   # 最后一列为评分
# 设置图形布局
par(mfrow = c(2, 1))
for (i in 1:5) {
  hist(gene_expression[, i], main = paste("Gene Expression", i), xlab = "Expression", ylab = "Frequency", col = "skyblue")
}
lines(density(na.omit(gene_expression)))
hist(score, xlab = "Score", ylab = "Frequency", col = "skyblue")
ks.test(scale(score), 'pnorm', alternative = 'two.sided') # 一定要加scale,小样本更好
#D = 0.098044, p-value < 2.2e-16,数据不符合正态分布，用斯皮尔曼相关系数计算


# 提取score行作为一个向量# 提取scscoreore行作为一个向量
# 使用转置t()是因为cor函数需要向量作为输入
score_vector <- t(labis["score", ])  
score_vector<-score_vector[(1:1068),]
# 初始化一个向量来存储相关性值
correlations <- numeric(nrow(marker))
# 对每个基因计算与score的相关性
for (i in 1:nrow(marker)) {
  gene_expression <- marker[i, ]    #默认皮尔森，斯皮尔曼要加method = "spearman", 
  correlations[i] <- cor(score_vector, gene_expression, method = "spearman", use = "complete.obs")
}
# 将相关性值与基因名称结合
correlation_df <- data.frame(Gene = rownames(marker), Correlation = correlations)
###lciber
lciber<-lciber[,-1]
col_sums <- colSums(lciber)
lciber<- rbind(lciber, col_sums)
rownames(lciber)[12]<-'score'
score_vector <- t(lciber["score", ])  #运行上面的循环即可
correlation_df2<- data.frame(Gene = rownames(marker), Correlation = correlations)
###lconsensus
rownames(lconsensus)<-lconsensus$cell_type;lconsensus<-lconsensus[,-1]
col_sums <- colSums(lconsensus)
lconsensus<- rbind(lconsensus, col_sums)
rownames(lconsensus)[7]<-'score'
score_vector <- t(lconsensus["score", ])  #运行上面的循环即可
correlation_df3<- data.frame(Gene = rownames(marker), Correlation = correlations)
###epic
rownames(lepic)<-lepic$cell_type;lepic<-lepic[,-1]
col_sums <- colSums(lepic)
lepic<- rbind(lepic, col_sums)
rownames(lepic)[5]<-'score'
score_vector <- t(lepic["score", ])  #运行上面的循环即可
correlation_df4<- data.frame(Gene = rownames(marker), Correlation = correlations)
###quan
rownames(lquan)<-lquan$cell_type;lquan<-lquan[,-1]
col_sums <- colSums(lquan)
lquan<- rbind(lquan, col_sums)
rownames(lquan)[6]<-'score'
score_vector <- t(lquan["score", ])  #运行上面的循环即可
correlation_df5<- data.frame(Gene = rownames(marker), Correlation = correlations)
###timer
rownames(ltimer)<-ltimer$cell_type;ltimer<-ltimer[,-1]
col_sums <- colSums(ltimer)
ltimer<- rbind(ltimer, col_sums)
rownames(ltimer)[4]<-'score'
score_vector <- t(ltimer["score", ])  #运行上面的循环即可
correlation_df6<- data.frame(Gene = rownames(marker), Correlation = correlations)
##xcell
col_sums <- colSums(lxcell)
lxcell<- rbind(lxcell, col_sums)
rownames(lxcell)[20]<-'score'
score_vector <- t(lxcell["score", ])  #运行上面的循环即可
correlation_df7<- data.frame(Gene = rownames(marker), Correlation = correlations)
file_names <- c("correlation_df.xlsx", "correlation_df2.xlsx", "correlation_df3.xlsx", 
                "correlation_df4.xlsx", "correlation_df5.xlsx", "correlation_df6.xlsx", 
                "correlation_df7.xlsx")
data_frames <- list(correlation_df, correlation_df2, correlation_df3, 
                    correlation_df4, correlation_df5, correlation_df6, 
                    correlation_df7)
# 批量导出数据框为xlsx文件
for (i in 1:length(data_frames)) {
  write.xlsx(data_frames[[i]], file = file_names[i])
}
##------------------相关系数计算（2）----------------------------------------------
colnames(lepic)<-gsub("-", ".", colnames(lepic))
colnames(lquan)<-gsub("-", ".", colnames(lquan))
colnames(lxcell)<-gsub("-", ".", colnames(lxcell))
score_labis <- tail(labis, 1)
score_lciber <- tail(lciber, 1)
score_lconsensus <- tail(lconsensus, 1)
score_lepic <- tail(lepic, 1)
score_lquan <- tail(lquan, 1)
score_ltimer <- tail(ltimer, 1)
score_lxcell <- tail(lxcell, 1)
score<-rbind(score_labis,score_lciber,score_lconsensus,score_lepic,score_lquan,score_ltimer,score_lxcell)
colnames(marker)<-gsub("-", ".", colnames(marker))
library(psych)  #函数要求行数一致
corr_matrix <- corr.test(t(marker), t(score), method = 'spearman', adjust = 'BH')
r <- corr_matrix$r    #相关系数矩阵
p <- corr_matrix$p    #p 值矩阵
###----------------------------plot----------------------------------------------
library(ggplot2)
library(reshape2)
correlationall<-as.data.frame(correlationall)
rownames(correlationall)<-correlationall$Gene
# 行名为基因名，列名为不同算法的名称
# 将数据框转换为长格式
long_data <- melt(correlationall, id.vars = "Gene", variable.name = 'Algorithm', value.name = 'Correlation')
# 计算相关系数值的范围
cor_range <- range(long_data$Correlation)
heatmap_plot <- ggplot(long_data, aes(x = Algorithm, y = Gene, fill = Correlation)) +
  geom_tile(color = "white") + # 添加瓷砖图层，并设置瓷砖边框为白色
  scale_fill_gradient2(high = "#B2182B", low = "#2166AC", mid = "white", midpoint = 0.5, 
                       limit = cor_range, space = "Lab", name = "Correlation") + # 修改颜色映射
  geom_text(aes(label = sprintf("%.3f", Correlation)), color = "black", size = 3) + # 添加相关系数的标签，保留小数点后三位
  theme_classic() + # 使用经典主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # 旋转x轴标签以避免重叠，并增大字体大小
        axis.text.y = element_text(size = 8), # 增大y轴标签字体大小
        axis.title = element_text(size = 10), # 增大轴标题字体大小
        legend.title = element_text(size = 10), # 增大图例标题字体大小
        legend.text = element_text(size = 8), # 增大图例文本字体大小
        panel.grid.major = element_blank(), # 移除主要网格线
        panel.grid.minor = element_blank(), # 移除次要网格线
        panel.border = element_blank(), # 移除面板边框
        panel.background = element_blank(), # 移除面板背景
        plot.title = element_text(size = 12, hjust = 0.5)) + # 增大图标题字体大小，并居中显示
  labs(title = "Correlation Heatmap", # 添加图标题
       x = "Algorithm", # 添加x轴标签
       y = "Gene") + # 添加y轴标签
  guides(fill = guide_colorbar(title.position = "top", title.vjust = 0.5)) # 修改颜色条图例标题位置
# 打印绘制的热图
print(heatmap_plot)

#-------------------一致性聚类（根据consensus评分）---------------------------------
library(ConsensusClusterPlus)
#ConsensusClusterPlus包提供了三种无监督聚类方式，分别是Hierarchical Clustering (层次聚类 - hc)
#Partitioning Around Medoids (基于中心点的划分 - pam)、K-means (K均值聚类 - km)
#distance="pearson",clusterAlg="hc")
#distance="euclidean",clusterAlg="km")
#distance="manhattan",clusterAlg="pam")
lconsensus1<-lconsensus [-7,]
boxplot(lconsensus1[,1:20])
lconsensus2 <- sweep(lconsensus1,1,apply(lconsensus1,1,median))#中位数归一化
boxplot(lconsensus2 [,1:20])
lconsensus2<- as.matrix(lconsensus2)
#多尝试几种聚类方式找出最优
#本数据：km>hc
#行为细胞类型，列为样本
ccres <- ConsensusClusterPlus(lconsensus2,
                              maxK=9,
                              reps=100,
                              pItem=0.8,
                              pFeature=1,
                              tmyPal = c("white","#C75D30"),
                              title='ConsensusCluster/',
                              clusterAlg="km",
                              distance="euclidean",
                              seed=23,
                              plot="png"
)
ccres <- ConsensusClusterPlus(lconsensus2,
                              maxK=9,
                              reps=100,
                              pItem=0.8,
                              pFeature=1,
                              tmyPal = c("white","#C75D30"),
                              title='ConsensusCluster/',
                              clusterAlg="pam",
                              distance="spearman",
                              seed=23,
                              plot="png"
)#最后选这个，色块更干净
iclres <- calcICL(ccres,title="ics of ssgsea res")
## PAC = Proportion of ambiguous clustering 模糊聚类比例
#通过比较不同K值下的PAC值，可以帮助确定最优的聚类数
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = ccres[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
## [1] 2
#提取结果
sample_subtypes <-as.data.frame( ccres[[2]][["consensusClass"]])
sample_subtypes$sample <- rownames(sample_subtypes)
colnames(sample_subtypes)[1]<-'cluster'
table(sample_subtypes)
save(sample_subtypes,file = '一致性聚类分组结果.RData')
#sample_subtypes
#1    2 
#1572 1476 
##------------------------------plot----------------------------------------------
###-------------------------------免疫浸润淋巴细胞差异---------------------
tlconsensus<-as_tibble(t(lconsensus))
rownames(tlconsensus)<-colnames(lconsensus) #行为样本列为细胞类型，矩阵为tibble格式
tlconsensus$ID<-rownames(tlconsensus)
tlconsensus<-tlconsensus[,-7]
tlconsensus %>%
  mutate(sample_subtypes = factor(sample_subtypes)) %>%
  pivot_longer(-c(ID, sample_subtypes), names_to = "cell_type", values_to = "value") %>%
  ggplot(aes(cell_type, value, fill = sample_subtypes)) +
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) + # 添加自定义颜色
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
###--------------所有免疫浸润细胞差异------------------------------
library(tidyr)
library(dplyr)
res_consensus_tme<-as.data.frame(t(res_consensus_tme))
res_consensus_tme$ID<-rownames(res_consensus_tme)
res_consensus_tme<-as_tibble(res_consensus_tme)

###-------------------------------样本信息热图-----------------------------------
#导入TCGA和META临床信息
clinical.indexed$submitter_id<-gsub("-", ".", clinical.indexed$submitter_id)
trainpam50$ID<-gsub("-", ".", trainpam50$ID)
colnames(trainpam50)[2]<-'sample'
tcgaage<-clinical.indexed[,c(2,40)]
metaage<-metaphe[,c(1,13)]
metaage$PATIENT_ID<-gsub("-", ".", metaage$PATIENT_ID)
metaage$AGE_AT_DIAGNOSIS<-round(metaage$AGE_AT_DIAGNOSIS)
colnames(metaage)[1]<-'sample'
colnames(tcgaage)[1]<-'sample'
colnames(metaage)[2]<-'Age'
colnames(tcgaage)[2]<-'Age'
age<-rbind(tcgaage,metaage)
Train_surv$sample<-rownames(Train_surv)
Train_surv$sample<-gsub("-", ".", Train_surv$sample)
cli<-merge(Train_surv,age,by='sample')
cli<-merge(cli,trainpam50,by='sample')
cli<-cli[,-3]
cli<-merge(cli,sample_subtypes,by='sample')
cli$cluster<-as.factor(cli$cluster)
cli$OS<-as.factor(cli$OS)
save(cli,file = '整合分组与临床信息.RData')
library(ComplexHeatmap)
columnAnno <- HeatmapAnnotation(status = cli$OS,
                                age = cli$Age,
                                pam500=cli$pam50,
                                cluster = cli$cluster
                                #,na_col = "white"
)
scaled_ssgsea <- scale(lconsensus)#标准化
scaled_ssgsea[scaled_ssgsea>2] <- 2#大于2设置为2
scaled_ssgsea[scaled_ssgsea< -2] <- -2#小于-2设置为-2
rownames(scaled_ssgsea)<- sub("_.*", "", rownames(scaled_ssgsea))
ComplexHeatmap::Heatmap(scaled_ssgsea, na_col = "white",show_column_names = F,
                        row_names_side = "left",name = "fraction",
                        column_order = c(colnames(lconsensus)[c(grep("1",cli$cluster),grep("2",cli$cluster))]),
                        column_split = cli$cluster, column_title = NULL,
                        cluster_columns = F,
                        top_annotation = columnAnno
)

#-------------------------------cluster分析--------------------------------------
##estimate
library(immunedeconv)
cancer="BRCA"#TCGA癌症缩写
brca_vector <- rep(cancer, ncol(adjustedtm))
res_estimate <- deconvolute(adjustedtm, method="estimate",tumor = T,indications = brca_vector)
write.table(res_estimate,"estimate.txt",sep="\t",row.names = F)
res_estimate<-as.data.frame(t(res_estimate))
rownames(res_estimate)<-gsub('-','.',rownames(res_estimate))
colnames(res_estimate) <- as.character(res_estimate[1,])
res_estimate <- res_estimate[-1,]
save(res_estimate,file = 'res_estimate.RData')
estimate<-cbind(res_estimate,sample_subtypes)
estimate<-estimate[,c(2,5)]
colnames(estimate)[1]<-'immune_score'
library(ggpubr) # 继承ggplot语法
library(patchwork) # 拼图包
library(ggsci) #配色包
estimate$cluster <- as.factor(estimate$cluster)
estimate$immune_score<-as.numeric(estimate$immune_score)#一定要改为数值型！
p <- ggplot(estimate, aes(x = cluster, y = immune_score, color = cluster)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-2000,4000)) +  # 设置可视范围为 -1000 到 3000
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
p
#不同亚型的免疫评分是否仍然是2>1
estimate$sample<-rownames(estimate)
espam<-merge(trainpam50,estimate,by='sample')
save(espam,file='计算免疫评分矩阵.RData')
p <- ggplot(espam, aes(x = cluster, y = immune_score, color = cluster)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(-2000,4000)) +  # 设置可视范围为 -1000 到 3000
  stat_compare_means(method = "wilcox.test") +  # 添加两组之间的p值
  facet_grid(cols = vars(pam50))  # 根据 pam50 列拆分子图
p
##-------------------------------------CYT评分---------------------------------------
#计算GZMA和PRF1表达的几何平均数
CYTmarker<- subset(adjustedtm, rownames(adjustedtm) %in% c("GZMA","PRF1"))
CYTmarker<-as.data.frame(t(CYTmarker))
library(psych) #加载包
for (i in 1:nrow(CYTmarker)) {
  CYTmarker$CYT[i] <- geometric.mean(c(CYTmarker[i,"GZMA"]+0.01,CYTmarker[i,"PRF1"]+0.01))
}  #计算每个样本的CYT
CYTmarker$sample<-rownames(CYTmarker)
CYTmarker$sample<-gsub('-','.',CYTmarker$sample)
CYTmarker<-merge(CYTmarker,trainpam50,by='sample')
CYTmarker<-merge(CYTmarker,sample_subtypes,by='sample')
save(CYTmarker,file = '计算CYT评分矩阵.RData')
p <- ggplot(CYTmarker, aes(x = cluster, y = CYT, color = cluster)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
p
#亚型CYT评分
p <- ggplot(CYTmarker, aes(x = cluster, y = CYT, color = cluster)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +  # 设置可视范围为 -1000 到 3000
  stat_compare_means(method = "wilcox.test") +  # 添加两组之间的p值
  facet_grid(cols = vars(pam50))  # 根据 pam50 列拆分子图
p


##----------------耗竭T细胞评分-----------------------
eT<-subset(adjustedtm, rownames(adjustedtm) %in% c("HAVCR2","LAG3","PDCD1","TIGIT"))
#计算marker的几何平均数
eT<-as.data.frame(t(eT))
for (i in 1:nrow(eT)) {
  eT$score[i] <- geometric.mean(c(eT[i,"HAVCR2"]+0.01,eT[i,"LAG3"]+0.01,eT[i,"PDCD1"]+0.01,eT[i,"TIGIT"]+0.01))
}  #计算每个样本的耗竭T细胞评分
eT<-as.data.frame(t(eT))
eT$sample<-rownames(eT)
eT$sample<-gsub('-','.',eT$sample)
eT<-merge(eT,sample_subtypes,by='sample')
save(eT,file = '计算耗竭T细胞评分矩阵.RData')
library(ggplot2)
eT$cluster<-as.factor(eT$cluster)
#发现2型耗竭T细胞评分显著高于1
p <- ggplot(eT, aes(x = cluster, y = score, color = cluster)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
p
#亚型耗竭T细胞评分
p <- ggplot(eT, aes(x = cluster, y = score, color = cluster)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35B2", "#8491B4B2")) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35B2", "#8491B4B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +  # 设置可视范围为 -1000 到 3000
  stat_compare_means(method = "wilcox.test") +  # 添加两组之间的p值
  facet_grid(cols = vars(pam50))  # 根据 pam50 列拆分子图