library(tidyverse)
library(impute) # 用于KNN填补药敏数据
library(pRRophetic) # 用于药敏预测
library(SimDesign) # 用于禁止药敏预测过程输出的信息
library(ggplot2) # 绘图
library(cowplot) #拼图
#分出肿瘤和正常样本
normsam <- colnames(TCGA_gset[,which(substr(colnames(TCGA_gset),14,15) == "11")])
nortpm <- as.matrix(TCGA_gset[,normsam])    
rowsam<-rownames(tpms)
nortpm<-nortpm[rowsam,]
#fpkm转tpm
expMatrix <- nortpm
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
nortpm <- apply(expMatrix,2,fpkmToTpm)
nortpm[1:3,]
colSums(nortpm)
colnames(nortpm)<-gsub("-",".",colnames(nortpm))





#读入评分数据group
# 标准化，把Score处理到0-1之间
rownames(group)<-group$ID
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
npps <- range01(group$Score)

# 创建样本信息
Sinfo <- data.frame(PPS = npps,
                    row.names = rownames(group),
                    stringsAsFactors = F)
# 把pps保存到文件
write.csv(Sinfo, "output_PPS.csv", quote = F)

library(ISOpureR)
runpure <- T # 如果想运行就把这个改为T
if(runpure) {
  set.seed(123)
  # Run ISOpureR Step 1 - Cancer Profile Estimation
  ISOpureS1model <- ISOpure.step1.CPE(tpms, nortpm)
  # For reproducible results, set the random seed
  set.seed(456);
  # Run ISOpureR Step 2 - Patient Profile Estimation
  ISOpureS2model <- ISOpure.step2.PPE(tpms,nortpm,ISOpureS1model)
  pure.tpms <- ISOpureS2model$cc_cancerprofiles
}

if(!runpure) {
  pure.tpms <- tpms
}






# 制作CTRP ctrpauc矩阵，保存到CTRP_ctrpauc.txt文件
ctrpauc <- read.table("C:/Users/xinyiwu/Desktop/Ful数据/训练6/F7/CTRP_ctrpauc_raw.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) 
# Supplementary Data Set 3
ctrpauc$comb <- paste(ctrpauc$index_cpd,ctrpauc$index_ccl,sep = "-")
ctrpauc <- apply(ctrpauc[,"area_under_sensitivity_curve",drop = F], 2, function(x) tapply(x, INDEX=factor(ctrpauc$comb), FUN=max, na.rm=TRUE)) # 重复项取最大ctrpauc
ctrpauc <- as.data.frame(ctrpauc)
ctrpauc$index_cpd <- sapply(strsplit(rownames(ctrpauc),"-",fixed = T),"[",1)
ctrpauc$index_ccl <- sapply(strsplit(rownames(ctrpauc),"-",fixed = T),"[",2)
ctrpauc <- reshape(ctrpauc,
               direction = "wide",
               timevar = "index_cpd",
               idvar = "index_ccl")
colnames(ctrpauc) <- gsub("area_under_sensitivity_curve.","",colnames(ctrpauc),fixed = T)
ccl_anno <- read_excel("C:/Users/xinyiwu/Desktop/ccl_anno.xlsx")
# Supplementary Data Set 1
cpd_anno <- read_excel("C:/Users/xinyiwu/Desktop/cpd.anno.xlsx")
# Supplementary Data Set 2
# 保存到文件    
write.table(ctrpauc,"CTRP_ctrpauc.txt",sep = "\t",row.names = F,col.names = T,quote = F)
# 加载药敏ctrpauc矩阵并进行数据处理
ctrp.ctrpauc <- read.table("CTRP_ctrpauc.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

## a. 移除缺失值大于20%的药物
ctrp.ctrpauc <- ctrp.ctrpauc[,apply(ctrp.ctrpauc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.ctrpauc)]
prism.ctrpauc <- prism.ctrpauc[,apply(prism.ctrpauc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.ctrpauc)]
## b. 移除CTRP数据里源自haematopoietic_and_lymphoid_tissue的细胞系
rmccl <- paste0("CCL",na.omit(ccl_anno[which(ccl_anno$ccle_primary_site == "haematopoietic_and_lymphoid_tissue"),"index_ccl"]))
rownames(ctrp.ctrpauc) <- paste0("CCL",rownames(ctrp.ctrpauc))
ctrp.ctrpauc <- ctrp.ctrpauc[setdiff(rownames(ctrp.ctrpauc),rmccl),]
## c. KNN填补缺失值
ctrp.ctrpauc.knn <- impute.knn(as.matrix(ctrp.ctrpauc))$data
prism.ctrpauc.knn <- impute.knn(as.matrix(prism.ctrpauc))$data
## d. 数据量级修正
ctrp.ctrpauc.knn <- ctrp.ctrpauc.knn/ceiling(max(ctrp.ctrpauc.knn)) # 参考Expression Levels of Therapeutic Targets as Indicators of Sensitivity to Targeted Therapeutics (2019, Molecular Cancer Therapeutics)
prism.ctrpauc.knn <- prism.ctrpauc.knn/ceiling(max(prism.ctrpauc.knn))
# 加载CCLE细胞系的表达谱，作为训练集
#下载链接https://depmap.org/portal/download/all/
ccl.expr <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 加载基因注释文件，用于基因ID转换
Ginfo <- read.table("Cell_lines_annotations_20181226.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)


library(stringr)
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Hs.eg.db)
#library(org.Hs.eg.db)人类
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
# 查看org.Hs.eg.db包支持的键类型
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
gene <- bitr(rownames(ccl.expr),fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
#data$Probe输入的基因  fromType：当前ID类型  toType：转换成什么ID   OrgDb：注释数据库
ccl.expr$ENSEMBL<-rownames(ccl.expr)
ccl.expr<-merge(ccl.expr,gene,by='ENSEMBL')
#去重
duplicated_symbols <- duplicated(ccl.expr$SYMBOL)
ccl.expr <- ccl.expr[!duplicated_symbols, ]
rownames(ccl.expr)<-ccl.expr$SYMBOL
ccl.expr <- ccl.expr[, -c(1, ncol(ccl.expr))]

keepgene <- apply(ccl.expr, 1, mad) > 0.5 # 保留表达值有效的基因
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) # 重置细胞系名
trainPtype <- as.data.frame(ctrp.ctrpauc.knn)
#ccl_anno中的细胞系名称有重复 只取第一次出现
ccl_anno <- ccl_anno[!duplicated(ccl_anno$cell_line_name), ]
ccl.name <- ccl.miss <- c() # 替换细胞系名
for (i in rownames(trainPtype)) {
  if(!is.element(gsub("CCL","",i),ccl_anno$index_ccl)) {
    cat(i,"\n")
    ccl.miss <- c(ccl.miss, i) # 没有匹配到的细胞系
    ccl.name <- c(ccl.name, i) # 插入未匹配的细胞系    
  } else {
    ccl.name <- c(ccl.name,  ccl_anno[which(ccl_anno$index_ccl == gsub("CCL","",i)),"cell_line_name"]) # 插入匹配的细胞系
  }
}
cpd.name <- cpd.miss <- c() # 替换药物名
for (i in colnames(trainPtype)) {
  if(!is.element(i,cpd_anno$index_cpd)) {
    cat(i,"\n")
    cpd.miss <- c(cpd.miss, i) # 没有匹配到的药物
    cpd.name <- c(cpd.name, i) # 插入未匹配的药物
  } else {
    cpd.name <- c(cpd.name,  cpd_anno[which(cpd_anno$index_cpd == i),"compound_name"]) # 插入匹配的药物
  }
}
rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] # 去除未匹配的细胞系
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] # 去除未匹配的药物
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) # 提取有表达且有药敏的细胞系
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]
# 测试集
keepgene <- apply(tpms, 1, mad) > 0.5 # 纯化的测试集取表达稳定的基因    
testExpr <- log2(tpms[keepgene,] + 1) # 表达谱对数化
# 取训练集和测试集共有的基因
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])



#----------ctrp预测药物敏感性------------
outTab <- NULL
# 循环很慢，请耐心
for (i in 1:ncol(trainPtype)) {
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于CTRP的AUC可能有0值，因此加一个较小的数值防止报错
  # 岭回归预测药物敏感性
  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
top.pps <- Sinfo[Sinfo$PPS >= quantile(Sinfo$PPS,probs = seq(0,1,0.1))[10],,drop=FALSE] # 定义上十分位的样本
bot.pps <- Sinfo[Sinfo$PPS <= quantile(Sinfo$PPS,probs = seq(0,1,0.1))[2],,drop=FALSE] # 定义下十分位的样本
ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  d <- rownames(ctrp.pred.auc)[i]
  a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值    
  b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}
candidate.ctrp <- ctrp.log2fc[ctrp.log2fc > 0.015] # 这里调整了阈值，控制结果数目
ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  d <- rownames(ctrp.pred.auc)[i]
  a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)])
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  ctrp.cor <- c(ctrp.cor,r)
  ctrp.cor.p <- c(ctrp.cor.p,p)
}
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.2]  # 这里我调整了阈值，控制结果数目    
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))


#------prism预测----------
prism.ctrpaucall <-  read.csv("D:/下载/secondary-screen-dose-response-curve-parameters.csv")
# 数据来自https://depmap.org/portal/download/ Drug sensitivity ctrpauc (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.ctrpauc[,1:5] # 前5列为细胞系注释信息
prism.ctrpauc <- prism.ctrpaucall[,c(3,9,12)]
#长数据转为宽数据
#直接转换会出现重复数据
library(tidyr)
extract_first_value <- function(x) {
  if (length(x) > 0) {
    return(x[[1]])
  } else {
    return(NA)
  }
}
# 使用 pivot_wider，并传递自定义函数给 values_fn 参数
prismauc <- pivot_wider(prism.ctrpauc,
                        names_from = "name", 
                        values_from = "auc",
                        values_fn = list(auc = extract_first_value))

#删除细胞名称列的遗漏行
prismauc <- prismauc[complete.cases(prismauc$ccle_name), ]
prismauc<-as.data.frame(prismauc)
rownames(prismauc)<-prismauc$ccle_name
prismauc<-prismauc[,-1]
## a. 移除缺失值大于20%的药物
prismauc <- prismauc[,apply(prismauc,2,function(x) sum(is.na(x))) < 0.2*nrow(prismauc)]
## c. KNN填补缺失值
prism.auc.knn <- impute.knn(as.matrix(prismauc))$data
## d. 数据量级修正
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))
trainprism <- as.data.frame(prism.auc.knn)
#表达矩阵的列名细胞系要与敏感性的行名细胞系一致
trainExprp <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExprp) <- sapply(strsplit(colnames(trainExprp),"_",fixed = T),"[",1) # 重置细胞系名
# 提取有表达且有药敏的细胞系
trainprism$ID<-rownames(trainprism)
trainprism$ID<-sub("_.*", "", trainprism$ID)
trainprism <- trainprism[!duplicated(trainprism$ID), ]
rownames(trainprism)<-trainprism$ID
trainprism<-trainprism[,-1292]
comccl1 <- intersect(rownames(trainprism),colnames(trainExprp)) 
trainExprp <- trainExprp[,comccl1]
trainprism <- trainprism[comccl1,]
trainExprp<-as.matrix(trainExprp)

##-----------预测----------
outTab1 <- NULL
# 循环很慢，请耐心
for (i in 1:ncol(trainprism)) {
  d <- colnames(trainprism)[i]
  tmp <- log2(as.vector(trainprism[,d]) + 0.00001) 
  # 岭回归预测药物敏感性
  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExprp,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
  outTab1 <- rbind.data.frame(outTab1,ptypeOut)
}
dimnames(outTab1) <- list(colnames(trainprism),colnames(testExpr))
prism.pred.auc <- outTab1

prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
  d <- rownames(prism.pred.auc)[i]
  a <- mean(as.numeric(prism.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
  b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  prism.log2fc <- c(prism.log2fc,log2fc)
}
candidate.prism <- prism.log2fc[prism.log2fc > 0.08] # 这里我调整了阈值，控制结果数目

prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
  d <- rownames(prism.pred.auc)[i]
  a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)])
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  prism.cor <- c(prism.cor,r)
  prism.cor.p <- c(prism.cor.p,p)
}
candidate.prism2 <- prism.cor[prism.cor < -0.2] 
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))


#------plot----------

#绘制相关性棒棒糖图
# 设置颜色
darkblue <- "#0772B9"
lightblue <- "#48C8EF"
cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) +
  theme(legend.position = "bottom",
        axis.line.y = element_blank())
cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)
p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +    
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) +
  theme(legend.position = "bottom",
        axis.line.y = element_blank())
ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)])
  b <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High Risk","Low Risk"),c(nrow(top.pps),nrow(bot.pps))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
#绘制显著差异箱线图
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) +
  geom_boxplot(aes(col = group),outlier.shape = NA) +
  geom_text(aes(drug, y=max(auc)),
            label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) +
  scale_fill_manual(values = c(darkblue, lightblue)) +
  scale_color_manual(values = c(darkblue, lightblue)) +
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank())
dat <- ggplot_build(p3)$data[[1]]
p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)    
prism.boxdata <- NULL
for (d in prism.candidate) {
  a <- as.numeric(prism.pred.auc[d,rownames(top.pps)])
  b <- as.numeric(prism.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                    data.frame(drug = d,
                                               auc = c(a,b),
                                               p = p,
                                               s = s,
                                               group = rep(c("High Risk","Low Risk"),c(nrow(top.pps),nrow(bot.pps))),
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)
p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) +
  geom_boxplot(aes(col = group),outlier.shape = NA) +
  geom_text(aes(drug, y=max(auc)),
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) +
  scale_fill_manual(values = c(darkblue, lightblue)) +
  scale_color_manual(values = c(darkblue, lightblue)) +
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]
p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
#进行拼图
plot_grid(p1, p3, p2, p4, labels=c("A", "", "B", ""),
          ncol=2, 
          rel_widths = c(2, 2)) #左右两列的宽度比例
#保存图片    
ggsave(filename = "drug.pdf",width = 8,height = 8)
