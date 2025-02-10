setwd('D:/ful数据/训练6/F2')
#-------------------------WGCNA分析-------------------------------------------------
BiocManager::install("preprocessCore")
library(WGCNA)
# 判断一下数据质量
# 通过goodSamplesGenes函数生成一个质量检查结果对象gsg
gsg <- goodSamplesGenes(t(adjustedtm), verbose = 3)#列为基因，行为样本
# `verbose` 参数被设置为 `3`，用于控制函数 `goodSamplesGenes` 的详细输出程度。
# `verbose = 0`：不产生任何输出，只返回结果，通常用于静默模式。
# `verbose = 1`：产生基本的信息输出，以提供一些关于函数执行进度的信息。
# `verbose = 2`：产生更详细的输出，可能包括一些中间步骤的信息。
# `verbose = 3`：产生最详细的输出，通常包括每个步骤的详细信息，用于调试和详细分析。
# 检查质量检查结果对象中的allOK属性，如果为FALSE，表示质量不好
gsg$allOK #[1] TRUE
##-----------------------------------构建样本的系统聚类树---------------------------
sampleTree <- hclust(dist(t(adjustedtm)), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
abline(h = 150, col = "red")#有150以上的离群样本
# 将层次聚类树中的样本分成不同的簇
clust <- cutreeStatic(sampleTree, cutHeight =150, minSize = 10)
table(clust)#有两个离群样本 剔除
# 0    1 
# 2  3046 
# 剔除离群样本，保留想要留下的样本
keepSamples <- clust == 1
wgcna<- adjustedtm[,keepSamples ] #16681 3046
#与临床信息对应
colnames(wgcna)<-gsub("-",".",colnames(wgcna))
rownames(cli)<-cli$sample
cmomsam<-intersect(colnames(wgcna),rownames(cli))
wgcnacli<-cli[cmomsam,]
save(wgcnacli,file="wgcnacli.RData")
save(wgcna,file = 'wgcna3046.RData')
##--------------------------网络构建和模块识别--------------------------------- 
###-----------------------------挑选最佳软阈值---------------------------------
# 设置 power 参数选择范围，可以自行修改设定
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# 选择最佳软阈值，获取各个阈值下的 R^2 和平均连接度
sft <- pickSoftThreshold(t(wgcna), powerVector = powers, verbose = 5)#行为样本列为基因
#看一下软阈值选择结果 
sft
#$powerEstimate
#[1] 4
#SFT.R.sq：Scale Free Topology Model Fit，表示拟合模型的拟合度,较高的值表示模型更好地拟合数据。
#slope：表示软阈值幂次与拟合模型的斜率，通常用于衡量网络的无标度拓扑结构。
#斜率绝对值越小，模型越能保持无标度拓扑结构。
# 绘制软阈值和拟合指标的关系图
cex1 = 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# 添加 R^2 水平线，使用 R^2 阈值为 0.90，官网建议最好是0.85或以上
abline(h = 0.90, col = "red")
# 绘制软阈值对平均连接度的影响图
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
##综合结果最终选择power=5

###--------------------------------构建加权共表达网络--------------------------------
# 一步法
allowWGCNAThreads()#RStudio不支持enableWGCNAThreads()
## Allowing multi-threading with up to 24 threads.
net = blockwiseModules(t(wgcna), power = 5,
                       TOMType = "unsigned", minModuleSize = 30,#每个模块最少有多少个基因
                       maxBlockSize = 16681, 
                       reassignThreshold = 0, mergeCutHeight = 0.25,#合并模块的阈值
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)
#power = sft$powerEstimate：这里指定了软阈值的取值
#maxBlockSize = nGenes：指定了最大的模块大小
#TOMType = "unsigned"：指定了共表达网络的类型，这里设置为 "unsigned"，表示构建无符号的网络
#minModuleSize = 30：指定了每个模块的最小大小
#reassignThreshold = 0：控制模块的重新分配
#mergeCutHeight = 0.25：用来合并基因模块,越大得到的模块就越小
#numericLabels = TRUE：表示模块将使用数值标签
#saveTOMs = F：用来控制是否保存共表达网络的拓扑重叠矩阵（TOM）
# 查看模块情况
table(net$colors) 
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
#6659 1318 1074  955  946  810  766  757  734  668  433  233  211  200  119  103   88   87   86   85   83 
#21   22   23   24   25 
#77   75   39   38   37 
# 根据模块中基因数目的多少，降序排列，依次编号，比如 1 为最大模块，模块中基因最多。
# 0 表示没有分入任何模块的基因。
# 将基因模块的标签转换为对应的颜色编码
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# 绘制层次聚类树
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 这里用不同的颜色来代表那些所有的模块，其中灰色默认是无法归类于任何模块的那些基因，
# 如果灰色模块里面的基因太多，那么前期对表达矩阵挑选基因的步骤可能就不太合适
#保存数据
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs,moduleLabels,moduleColors,geneTree,
     file = "wgcna-02-networkConstruction-auto.rdata")



##-----------------------------计算模块与感兴趣的性状之间的相关性和显著性----------
# 颜色标签
moduleLables <- net$colors
moduleColors <- labels2colors(net$colors)
# ME值，也就是获取eigengenes，每个ME代表一个模块
MEs <- net$MEs
head(MEs)[1:5, 1:5]
# 将基因模块与性状进行关联
# 获取eigengenes，用颜色标签计算ME值
MEList <-  moduleEigengenes(t(wgcna), colors = moduleColors)
save(MEList,file='MElist.RData')
MEs0 <- MEList$eigengenes
# 查看用颜色标签计算的ME值
head(MEs0)[1:5, 1:5] #  MEblack       MEblue      MEbrown       MEcyan   MEdarkgreen
# 排序
MEs <- orderMEs(MEs0)
head(MEs)[1:5, 1:5] #MEdarkred      MEblack        MEred MEdarkturquoise MEgreenyellow
# 计算每个模块和每个性状之间的相关性
wgcnacli<-wgcnacli[,-3]#删除样本标签和pam50分型
wgcnacli$OS<-as.numeric(wgcnacli$OS)
wgcnacli$Age<-as.numeric(wgcnacli$Age)
wgcnacli$cluster<-as.numeric(wgcnacli$cluster)
moduleTraitCor <- cor(MEs, wgcnacli, use = "p")


#这里发现相关性都不是很强，尝试合并相似模块
#计算模块的eigengenes（其实就是第一主成分）
#再通过相关性进行聚类，达到量化模块间相似性的目的，方便进行模块合并
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#选择切割高度进行合并，注意，这里的相关性 = 1-切割高度，
#比如切割高度是0.25，那就是选择相关性大于0.75的模块进行合并。
MEDissThres = 0.3 # 要根据你自己画出来的图选择
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(wgcna),moduleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
moduleTraitCor1 <- cor(mergedMEs, wgcnacli, use = "p")
#与合并之前差异不大，算了
sizeGrWindow(12, 9)
pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(moduleColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



# 计算显著性
nSamples = nrow(t(wgcna))
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)
# 可视化相关性和P值
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(1, 15, 3, 3))
# 在绘制热图之前，保存绘图参数
# 绘制热图
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(wgcnacli),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,  # 因为我们手动调整了边距
  cex.text = 0.6,
  zlim = c(-0.8, 0.8), 
  main = paste("Module-trait relationships")
)
save(moduleTraitCor,moduleTraitPvalue,file='相关性热图绘制数据.RData')
#modulered相关性最强
module = "red"
module_genes <- colnames(t(wgcna))[moduleColors==module]
module_genes<-as.data.frame(module_genes)
save(modele_genes,file='modulered.RData')


#-------------------计算GS和MM-------------------
# 获取模块名称
modNames <- substring(names(MEs), 3)
modNames
# 计算模块与基因的相关性矩阵
geneModuleMembership <- as.data.frame(cor(t(wgcna), MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "");
geneModuleMembership[1:5, 1:5]
# 计算性状与基因的相关性矩阵 
# 只有连续型性状才能进行计算，如果是离散变量，在构建样本表时就转为0-1矩阵。
geneTraitSignificance <- as.data.frame(cor(t(wgcna), wgcnacli[3], use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(wgcnacli[3]), sep = "")
names(GSPvalue) <- paste("p.GS.", names(wgcnacli[3]), sep = "")
head(geneTraitSignificance)
# 保存结果
save(geneModuleMembership, geneTraitSignificance, MMPvalue, GSPvalue, file = "GSMM.RData")
# 最后把两个相关性矩阵联合起来，指定感兴趣模块进行分析
module = "red"
pheno = "cluster"
modNames = substring(names(MEs), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno, colnames(wgcnacli))

# 获取模块内的基因
moduleGenes <- moduleColors == module
# 可视化基因与模块、表型的相关性，绘制散点图
sizeGrWindow(6, 6);
par(mfrow = c(1,1));
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, module_column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = paste("Gene significance for cluster"),
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  col ="#E64B35B2",pch = 16   # 将颜色设置为红色
)
abline(h=0.6,v=0.8,col="#8491B4B2",lty = 2,lwd=1.6)
#提取GS和MM都很高的基因
tmp1 <- rownames(geneModuleMembership)[abs(geneModuleMembership[moduleGenes, module_column])>0.8]
tmp2 <- rownames(geneTraitSignificance)[abs(geneTraitSignificance[moduleGenes, 1])>0.6]
genes <- unique(intersect(tmp1,tmp2))
length(genes)
genes<-as.data.frame(genes)
save(genes,file='285hubgene.RData')


