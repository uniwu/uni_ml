#训练6的机器学习筛选
#-------------------------mechine learning-------------------------------
work.path <- "D:/Ful数据/训练6/PrognosticML_4.0"; setwd(work.path) 
code.path <- file.path(work.path, "Codes") 
data.path <- file.path(work.path, "InputData") 
res.path <- file.path(work.path, "Results") 
fig.path <- file.path(work.path, "Figures") 
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
source(file.path(code.path, "ML.R"))
FinalModel <- c("panML", "multiCox")[2]
comsam<-intersect(rownames(coxgene),colnames(allcox))
trainos<-allcox[,c(1:2)]
#数据标准化,行为样本列为基因
trainexp<-allcox[,comsam]
trainset= scaleData(data = trainexp, centerFlags = T, scaleFlags = T) 
testset<-t(testset)
names(x = split(as.data.frame(testset), f = testos$Cohort))
testset= scaleData(data = testset, cohort = testos$Cohort, centerFlags = T, scaleFlags = T)


# Model training and validation -------------------------------------------
## method list --------------------------------------------------------
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 4 #设置最少变量数
timeVar = "OS.time"; statusVar = "OS" 

# Pre-training 
#进行预训练的主要目的是为了选择重要的变量或特征，这些特征将用于后续的模型训练。
#这个过程被称为特征选择，是机器学习和统计建模中的一个常见步骤，尤其是在处理具有大量变量的高维数据时。
#预训练可以帮助减少模型的复杂性，提高模型的解释性，同时有可能提高模型的泛化能力。
Variable = colnames(trainset)
# - 使用 strsplit() 函数按照 "+" 分割预训练方法的字符串
#得到一个列表，每个元素是一个包含预训练方法的字符向量
preTrain.method =  strsplit(methods, "\\+")
# - 使用 lapply() 函数遍历每个列表元素（即每个预训练方法字符向量）
#并将其反转（rev() 函数）并删除第一个元素（[-1]）
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
# 去除重复的方法名称
preTrain.method = unique(unlist(preTrain.method))
preTrain.method#选择部分方法进行pretrain


#os必须为数值型或逻辑值

set.seed(seed = 23) 
preTrain.var <- list()
#长时间无结果，取子集看看有没有问题
#a<-trainset[1:100,]
#b<-trainos[1:100,]  #子集可以运行，代码无误
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_expr = trainset, # 训练数据的表达式矩阵
                                 Train_surv = trainos,  # 包含生存时间和事件状态的数据框
                                 mode = "Variable",      # 模式设置为"Variable"以选择变量
                                 timeVar = "OS.time", statusVar = "OS")  
}
preTrain.var[["simple"]] <- colnames(trainset)
saveRDS(preTrain.var, file.path(res.path, "preTrain.var.rds")) 
#遍历methods列表中的每个机器学习方法
#并根据每种方法选择变量、训练模型，并最终保存模型结果
model <- list() 
set.seed(seed = 23)
for (method in methods) { 
  cat(match(method, methods), ":", method, "\n") 
  method_name = method
  method_split = strsplit(method, "\\+")[[1]] 
  
  if (length(method_split) == 1) {
    method_split <- c("simple", method_split)
  }
  
  selected.var = preTrain.var[[method_split[1]]]
  
  if (length(selected.var) <= min.selected.var) {
    cat("Method", method_name, "has less than or equal to", min.selected.var, "variables. Skipping...\n")
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method_split[2], 
                                  Train_expr = trainset[, selected.var, drop = FALSE], 
                                  Train_surv = trainos, 
                                  mode = "Model",       
                                  timeVar = "OS.time", statusVar = "OS")

    if (is.null(model[[method_name]]) || length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      cat("Model for", method_name, "is NULL or has less than or equal to", min.selected.var, "selected variables after training. Setting to NULL...\n")
      model[[method_name]] <- NULL
    }
  }
}
saveRDS(model, file.path(res.path, "mymodel.rds")) 

## Evaluate the model -----------------------------------------------------
methodsValid <- names(model)
RS_list <- list()# 创建一个空列表用于存储风险评分结果
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]],  # 使用CalRiskScore函数计算风险评分
                                    new_data = rbind.data.frame(trainset,testset), 
                                    type = "lp") 
  
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "myRS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) 

fea_list <- list() # 创建一个空列表用于存储特征变量
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) 
}


fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) 
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  
write.table(fea_df, file.path(res.path, "myfea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)


testos$OS<-as.numeric(testos$OS)
testos$OS.time<-as.numeric(testos$OS.time)
Cindexlist <- list()
for (method in methods){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = testset, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_surv = testos, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                  Train_expr = trainset, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_surv = trainos, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                  Train_name = "TCGA+METABRIC", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  cohortVar = "Cohort", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  timeVar = "OS.time", # 用于评估的生存时间，必须出现在Test_surv中
                                  statusVar = "OS") # 用于评估的生存状态，必须出现在Test_surv中
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "Cindex_mat1.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)
# Plot --------------------------------------------------------------------
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) 
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] 
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] 

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
names(CohortCol) <- colnames(Cindex_mat)


cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, 
                    avg_Cindex = avg_Cindex, 
                    CohortCol = CohortCol, 
                    barCol = "steelblue", 
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F) 


pdf(file.path(fig.path, "myheatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

#------------------------测试集提取---------------------------
library(AnnoProbe)
library(GEOquery)
library(affy)
library(Biobase)
library(ggExtra)
library(ggpubr)
library(limma)
##-----------------------------GSE20685-----------------------------------------
options( 'download.file.method.GEOquery' = 'libcurl' )
gset <- getGEO('GSE20685', destdir="D:/Ful数据/F3/测试集/GSE20685",
               AnnotGPL = T,     
               getGPL = T) 
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol` <- data.frame(sapply(gpl$`Gene symbol`, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors = FALSE)[, 1]
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
rownames(exprset_unique)[rownames(exprset_unique) == "MIR612"] <- "NEAT1"
comsam<-intersect(rownames(exprset_unique),colnames(trainset))
GSE20685<-exprset_unique[comsam,]
save(GSE20685,file = 'GSE20685.RData')
pdata<-pData(gset[[1]])
pdata1<-as.data.frame(pdata[,c(49,51,54,55,60)])
score_20685<-cbind(pdata1,score_20685)
save(score_20685,file = 'score_20685.RData')

##--------------------------------GSE48390-------------------------------------------
options( 'download.file.method.GEOquery' = 'libcurl' )
gset <- getGEO('GSE48390', destdir="D:/Ful数据/F3/测试集/GSE48390",
               AnnotGPL = TRUE,
               getGPL = TRUE)
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]), stringsAsFactors=F)[,1] 
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
rownames(exprset_unique)[rownames(exprset_unique) == "MIR612"] <- "NEAT1"
comsam<-intersect(rownames(exprset_unique),colnames(trainset))
GSE48390<-exprset_unique[comsam,]
comsam<-intersect(colnames(GSE48390),rownames(testos))
GSE48390<-GSE48390[,comsam]
save(GSE48390,file = 'GSE48390.RData')
pdata<-pData(gset[[1]])
pdata1<-as.data.frame(pdata[,c(42)])
score_88770<-cbind(pdata1,score_88770)
save(score_88770,file = 'score_88770.RData')
##-----------------------------------GSE88770-------------------------------------
gset <- getGEO('GSE88770', destdir="D:/Ful数据/F3/测试集/GSE88770",
               AnnotGPL = T,     
               getGPL = T) 
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]), stringsAsFactors=F)[,1] 
exprset <- data.frame(exprs(gset[[1]]))
exprset<-as.data.frame(exprset)
exprset$ID<-rownames(exprset)
exprset_symbol<-merge(exprset,gpl,by="ID")
exprset_symbol<-na.omit(exprset_symbol)
table(duplicated(exprset_symbol$`Gene symbol`))
exprset_unique<-avereps(exprset_symbol[,-c(1,ncol(exprset_symbol))],ID=exprset_symbol$`Gene symbol`)
table(duplicated(rownames(exprset_unique)))
exprset_unique<-as.data.frame(exprset_unique)
rownames(exprset_unique)[rownames(exprset_unique) == "MIR612"] <- "NEAT1"
comsam<-intersect(rownames(exprset_unique),colnames(trainset))
GSE88770<-exprset_unique[comsam,]
comsam<-intersect(colnames(GSE88770),rownames(testos))
GSE88770<-GSE88770[,comsam]
save(GSE88770,file = 'GSE88770.RData')
pdata<-pData(gset[[1]])
pdata1<-as.data.frame(pdata[,c(62,65,66)])
score_88770<-cbind(pdata1,score_48390)
save(score_48390,file = 'score_48390.RData')
##--------------------GSE42568-----------------------------------
rownames(GSE42568)[rownames(GSE42568) == "MIR612"] <- "NEAT1"
comsam<-intersect(rownames(GSE42568),colnames(trainset))
GSE42568<-GSE42568[comsam,]
save(GSE42568,file = 'GSE42568.RData')
gset <- getGEO('GSE42568', destdir="D:/Ful数据/F3/测试集/GSE42568",
               AnnotGPL = T,     
                    getGPL = T) 
pdata<-pData(gset[[1]])  
pdata1<-as.data.frame(pdata[,c(49,44,43,42,41)])
pdata1<-pdata1[c(18:121),]
score_42568<-cbind(pdata1,score_42568)
save(score_42568,file = 'score_42568.RData')
##--------------------GSE16446----------------------
comsam<-intersect(rownames(GSE16446),colnames(trainrsf))
GSE16446<-GSE16446[comsam,]
GSE16446<-t(GSE16446)
GSE16446<-cbind(GSE16446os,GSE16446)
GSE16446$OS.time<-as.numeric(GSE16446$OS.time)/365
test1<-cbind(test,GSE16446)
##--------------------测试集表达矩阵-----------------------------
testset<-cbind(GSE20685,GSE48390,GSE88770,GSE42568)
library(sva)
batch<-c(rep('GSE20685',327),rep('GSE48390',81),rep('GSE88770',117),rep('GSE42568',104))#批次信息
testexp_numeric <- apply(testset, 2, as.numeric)
adjustedtest <- ComBat(dat = testexp_numeric, batch = batch)
save(testset,file = '去批次未标准化的testexp.RData')
GSE42568os$Cohort<-"GSE42568"
testos<-rbind(testos,GSE42568os)
testos <- subset(testos, Cohort != "GSE37751")
save(testos,file = 'testos.RData')

#----------------------------------Coxboost-----------------------------
#10折交叉验证
# 使用optimCoxBoostPenalty寻找最佳penalty参数
trainset<-t(trainset)
set.seed(23)
sfInit(parallel = TRUE, cpus = 12)
pen <- optimCoxBoostPenalty(trainos$OS.time,
                            trainos$OS,
                            as.matrix(trainset),
                            trace=TRUE,
                            parallel = TRUE) 
sfStop()
pen$penalty #[1] 11628
#并行运算
# 初始化并行计算框架
set.seed(23)
sfInit(parallel = TRUE, cpus = 12)
cv_result <- cv.CoxBoost(time = trainos$OS.time,
                         status = trainos$OS,
                         x = as.matrix(trainset),
                         maxstepno=500,
                         penalty =pen$penalty,
                         K = 10, # 10折交叉验证
                         type="verweij",
                         parallel = TRUE)
# 关闭并行计算框架
sfStop()
cv_result$optimal.step #[1] 199
# 使用得到的最佳penalty和stepno参数构建最终的CoxBoost模型
final_model <- CoxBoost(time = trainos$OS.time,
                        status = trainos$OS,
                        x = as.matrix(trainset),
                        penalty = pen$penalty ,
                        stepno = cv_result$optimal.step)
summary(final_model)
save(final_model,file = 'coxboostmodel.RData')
#parameter estimates > 0:
#EP400NL, ZNF148, ATG12, PRAME, C12orf57, AATF, SPDL1, SIX2, TOMM34 
#parameter estimates < 0:
#IL27RA, IGSF10, CIR1, PXK, LAMA3, IQUB, CD59, SPATA7, MORN3, VPS37A 
par(mar=c(5, 5, 5, 5))  # 调整边距以确保变量名不重叠
plot(final_model, varnames=TRUE)
step.logplik<-predict(final_model,newdata=as.matrix(trainset),
                      newtime=trainos$OS.time,
                      newstatus=trainos$OS,
                      at.step=0:199,
                      type="logplik")
plot(step.logplik)
coxboost<-as.data.frame(final_model$xnames[final_model$coefficients[which.max(step.logplik),]!=0])
save(coxboost,file = 'coxboostgene.RData')
rownames(coxboost)<-coxboost[,1]
comsam<-intersect(rownames(coxboost),colnames(trainset))
trainrsf<-trainset[,comsam]

#---------------------RSF---------------------------
#调参
#合并生存和表达
trainrsf<- cbind(trainos,trainrsf)  
save(trainrsf,file = 'trainrsf.RData')
set.seed(23)
timeVar = "OS.time"; statusVar = "OS" 
fit <- rfsrc(formula = formula(paste0("Surv(", timeVar, ", ", statusVar, ")", "~.")),
             data = trainrsf,
             ntree = 1000, nodesize = 5,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T)

print(fit)
#调参
plot(fit)
which.min(fit$err.rate)# 1000
tune.nodesize(Surv(OS.time,OS) ~ ., trainrsf) #optimal nodesize: 8
set.seed(23)
fit <- rfsrc(formula = formula(paste0("Surv(", timeVar, ", ", statusVar, ")", "~.")),
             data = trainrsf,
             ntree = 1000, nodesize = 8,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T)
print(fit)
plot(fit)
save(fit,file = 'rsffit.RData')

##----------------------------------变量筛选-------------------------------------
library(tidyverse)
importance_gene <- data.frame(fit$importance) %>% 
  rownames_to_column("gene") %>% 
  arrange(- fit.importance) %>% 
  head(19)
save(importance_gene,file = 'rsfimportance_gene.RData')

library(ggRandomForests)
#首先查看VIMP法的变量重要性(其实就是importance_gene)
gg_dta <- gg_vimp(fit)
plot(gg_dta)
#最小深度法查看变量重要性
gg_dta1<- gg_minimal_depth(fit)
gg_dta1<-as.data.frame(gg_dta1[["varselect"]])
save(gg_dta1,file = 'vimp+depth.RData')
#两种结合
gg_dta <- gg_minimal_vimp(fit)

#绘制基因重要性图
color_palette <- c("#F0F0F0", "#B2182B")
ggplot(data = importance_gene, aes(x = reorder(gene, fit.importance), y = fit.importance, fill = fit.importance)) +
  geom_col() +
  labs(x = "Gene", y = "Importance", title = "Gene Importance in RSF Model") +
  scale_fill_gradient(low = color_palette[1], high = color_palette[2], guide = "none") +  # 设置颜色从略带灰色调的浅色到深色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),  # 去除网格线
        legend.position = "none") +
  coord_flip()  # 横向显示变量

#采用km法计算Brier score
bs_km <- get.brier.survival(fit,  cens.model = "km")$brier.score
#采用rfsrc法计算Brier score
bs_rsf <- get.brier.survival(fit,  cens.model = "rfsrc")$brier.score
#绘制图形并比较
pdf("brierscore.pdf")
par(pty = "s")
plot(bs_km, type = "s", col = "#E64B35B2", lwd = 3)
lines(bs_rsf, type = "s", col ="#8491B4B2", lwd = 3)
legend("bottomright",
       legend = c("cens.model" = "km",
                  "cens.moedl" = "rfs"),
       fill = c('#F39B7F', '#4DBBD5'))
dev.off()

##----------------偏依赖图---------------------------
# 计算图形布局的行数和列数，这里假设最多显示16个图形为例
# 可以根据实际变量的数量调整这些参数
nrow <- 4
ncol <- 4
nplots <- nrow * ncol
# 设置图形布局参数
par(mfrow=c(nrow, ncol))
# 获取模型中的变量名称
variable_names <- names(fit$xvar)
# 只绘制前nplots个变量的偏依赖图
for(var_name in variable_names[1:min(length(variable_names), nplots)]) {
  # 计算偏依赖
  partial_obj <- randomForestSRC::partial(fit,
                                          partial.xvar = var_name,
                                          partial.type = "mort",
                                          partial.values = fit$xvar[[var_name]],
                                          partial.time = fit$time.interest)
  
  # 提取数据
  pdta <- get.partial.plot.data(partial_obj)
  
  # 绘图
  plot(lowess(pdta$x, pdta$yhat, f = 1/3),
       type = "l", xlab = var_name, ylab = "Adjusted mortality", main = paste("PDP for", var_name))
}
# 重置图形布局参数
par(mfrow=c(1, 1))

#------------------------训练集的生存曲线和ROC----------------------------------
#计算风险评分
#day展示不清，转为year
trainrsf$OS.time<-trainrsf$OS.time / 365
library(survminer)
score_train <- data.frame(trainrsf[,c(1,2)],Score=fit$predicted)
save(score_train,file = 'scoretrain.RData')
cut <- surv_cutpoint(score_train,'OS.time','OS','Score')
summary(cut)
#      cutpoint statistic
#Score 75.05934  38.39938
save(cut,file = 'cut.RData')
#绘制cutoff值
plot(cut)
cat <- surv_categorize(cut)
save(cat,file = 'traingroup.RData')
library(survival)
survifit <- survfit(Surv(OS.time,OS)~Score,cat)#拟合生存曲线模型
#mytheme <- theme_survminer(font.legend = c(14,'plain', 'black'),
#                           font.x = c(14,'plain', 'black'),
#                          font.y = c(14,'plain', 'black')) ## 自定义主题
diff=survdiff(Surv(OS.time,OS)~Score,cat)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
##-------------------K-M曲线------------------
pdf("KM_train.pdf", width = 6, height = 6)
# KM曲线
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by = 7.5,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
# 关闭PDF设备
dev.off()
##-------------------------时间依赖性ROC---------------------
library(timeROC)
col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(score_train$OS.time,score_train$OS,score_train$Score,
                  cause = 1,weighting = 'marginal',
                  times=c(3,5,10), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_train.pdf", width = 6, height = 6)
# 设置图形参数以提高美观性
par(cex.axis = 1, cex.lab = 1.5, # 轴标签大小
    font.axis = 2, font.lab = 2, # 字体风格
    mar = c(5, 5, 2, 2) + 0.1, # 边距
    mgp = c(3, 1, 0), # 轴标题和轴线的距离
    tck = -0.02) # 轴刻度长度
# 假设ROC_rt是一个包含时间依赖性ROC数据的对象
# 绘制三条ROC曲线，分别对应不同时间点
plot(ROC_rt, time = 3, col = "#E64B35AA", main = "", lwd = 3, ylim = c(0.5, 1), xlim = c(0, 1))
plot(ROC_rt, time = 5, col = "#8491B4AA", add = TRUE, main = "", lwd = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(ROC_rt, time = 10, col = '#00A087B2', add = TRUE, main = "", lwd = 3)
# 添加图例
legend('bottomright', 
       legend = c(paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
                  paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
                  paste0('AUC at 10 years: ', sprintf("%.03f", ROC_rt$AUC[3]))), 
       col = c("#E64B35AA", "#8491B4AA", '#00A087B2'), 
       lwd = 3, 
       bty = 'n', 
       text.font = 2,
       cex = 0.9) # 图例文字大小
# 关闭图形设备
dev.off()
##----------------------------时间依赖性Cindex---------------------------------
###-----------------------------TCGA------------------------------------
setwd("D:/Ful数据/训练6/F3")
tcli<-as.data.frame(clinical.indexed[,c(2,4,19,22,23,36,40)])
tcli$submitter_id<-gsub('-','.',tcli$submitter_id)
rownames(tcga_pam50)<-gsub('-','.',rownames(tcga_pam50))
rownames(tcli)<-tcli$submitter_id
comsam<- intersect(rownames(score_train),rownames(tcli))
tcli<-tcli[comsam,]
score_train<-score_train[comsam,]
tcga_pam50<-tcga_pam50[comsam,]
tcli<-cbind(score_train,tcli)
tcli$ID<-rownames(tcli)
tcli<-merge(tcli,trainpam50,by='ID')

tcli$ajcc_pathologic_stage<-as.factor(tcli$ajcc_pathologic_stage)
sttcga<-coxph(Surv(OS.time, OS) ~ ajcc_pathologic_stage, data = tcli, x = TRUE)
tcli$ajcc_pathologic_n<-as.factor(tcli$ajcc_pathologic_n)
ntcga<-coxph(Surv(OS.time, OS) ~ ajcc_pathologic_n, data = tcli, x = TRUE)
tcli$ajcc_pathologic_t<-as.factor(tcli$ajcc_pathologic_t)
ttcga<-coxph(Surv(OS.time, OS) ~ ajcc_pathologic_t, data = tcli, x = TRUE)
tcli$ajcc_pathologic_m<-as.factor(tcli$ajcc_pathologic_m)
mtcga<-coxph(Surv(OS.time, OS) ~ ajcc_pathologic_m, data = tcli, x = TRUE)
tcli$race<-as.factor(tcli$race)
ratcga<-coxph(Surv(OS.time, OS) ~ race, data = tcli, x = TRUE)
pamtcga<-coxph(Surv(OS.time, OS) ~ pam50, data = tcli, x = TRUE)
scoretcga<-coxph(Surv(OS.time, OS) ~ Score, data = tcli, x = TRUE)
times <- seq(min(tcli$OS.time), max(tcli$OS.time), by = 2)
pk<- cindex(list(sttcga,ttcga,ntcga,mtcga,ratcga,pamtcga,scoretcga),
            Surv(OS.time, OS) ~ 1, 
            data = tcli,
            eval.times=times)
library(stats)
# 进行 Mann-Whitney U 检验
wilcox_test <- wilcox.test(pk$AppCindex$coxph, pk$AppCindex$coxph5)
print(wilcox_test)  #7.393e-15
# 保存图形到PDF文件
pdf("c-indextcga.pdf", width = 6, height = 6) 
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)
plot(pk, col = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5','#7E6148B2', '#00A087B2'),
     lwd = 2, xlab = "Time (years)", ylab = "C-index", 
     ylim = c(0.4, 0.8),
     xlim = c(1, max(tcli$OS.time)), main = "Time-dependent C-index for Models")
legend('bottomright', legend = c("Stage",'T','N','M','Race','Pam50',"Score"), 
       col =  c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5','#7E6148B2', '#00A087B2'),
       lwd = 3, bty = 'n', text.font = 2, cex = 0.9)
dev.off()


###-----------------------------METABRIC------------------------------------
mcli<-as.data.frame(metaphe[,c(1,3,5,9,13,20)])
mcli$PATIENT_ID<-gsub('-','.',mcli$PATIENT_ID)
rownames(mcli)<-mcli$PATIENT_ID
comsam<- intersect(rownames(score_train),rownames(mcli))
mcli<-mcli[comsam,]
score_train<-score_train[comsam,]
rownames(trainpam50)<-gsub('-','.',rownames(trainpam50))
trainpam50<-trainpam50[comsam,]
mcli1<-cbind(score_train,mcli)
mcli2<-cbind(mcli1,trainpam50)
save(mcli2,file="mcli2.RData")

sme<-coxph(Surv(OS.time, OS) ~ Score, data = mcli2, x = TRUE)
npime<-coxph(Surv(OS.time, OS) ~ NPI, data = mcli2, x = TRUE)
mcli2$CHEMOTHERAPY<-as.factor(mcli2$CHEMOTHERAPY)
cheme<-coxph(Surv(OS.time, OS) ~ CHEMOTHERAPY, data = mcli2, x = TRUE)
mcli2$HORMONE_THERAPY<-as.factor(mcli2$HORMONE_THERAPY)
home<-coxph(Surv(OS.time, OS) ~ HORMONE_THERAPY, data = mcli2, x = TRUE)
ageme<-coxph(Surv(OS.time, OS) ~ AGE_AT_DIAGNOSIS, data = mcli2, x = TRUE)
mcli2$RADIO_THERAPY<-as.factor(mcli2$RADIO_THERAPY)
rame<-coxph(Surv(OS.time, OS) ~ RADIO_THERAPY, data = mcli2, x = TRUE)
pamme<-coxph(Surv(OS.time, OS) ~ pam50, data = mcli2, x = TRUE)

times <- seq(min(mcli2$OS.time), max(mcli2$OS.time), by = 0.1)
pk1<- cindex(list(npime,cheme,home,ageme,rame,pamme,sme),
            Surv(OS.time, OS) ~ 1, 
            data = mcli2,
            eval.times=times)

pdf("cindexmeta.pdf", width = 6, height = 6) 
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)
plot(pk1, col = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5','#7E6148B2', '#00A087B2'),
     lwd = 2, xlab = "Time (years)", ylab = "C-index", 
     ylim = c(0.5, 0.8),
     xlim = c(min(mcli2$OS.time), max(mcli2$OS.time)), main = "Time-dependent C-index for Models")
# 添加图例
legend('bottomright', legend = c("Stage",'T','N','M','Race','Pam50',"Score"), 
       col =  c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5','#7E6148B2', '#00A087B2'),
       lwd = 3, bty = 'n', text.font = 2, cex = 0.9)
dev.off()

pdf("cindexmeta.pdf", width = 6, height = 6) 
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)


#---------------------测试集的生存和ROC曲线--------------------
comsam<- intersect(colnames(testset), colnames(trainrsf))
testset<-testset[, comsam]
test<-cbind(testos,testset)
test <- subset(test, Cohort != "GSE20711")
save(test,file = 'test.RData')
test$OS.time<-as.numeric(test$OS.time) / 365
#先删除分组数据
test<-test[,-3]
test$OS<-as.numeric(test$OS)

rsf_v<- predict(fit,newdata = test,proximity=T)
score_test <- data.frame(test[,c(1,2)],Score=rsf_v$predicted)
save(score_test,file = 'scoretest.RData')
cut <- surv_cutpoint(score_test,'OS.time','OS','Score')
plot(cut)
cat <- surv_categorize(cut)
save(cat,file = 'testgroup.RData')
survifit <- survfit(Surv(OS.time,OS)~Score,cat)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~Score,cat)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
##---------绘制测试集的K-M曲线-------
pdf("KM_test.pdf",width = 6,height =6)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by = 5,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
##---------------------测试集ROC--------------------
library(timeROC)
col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(score_test$OS.time,score_test$OS,score_test$Score,
                  cause = 1,weighting = 'marginal',
                  times=c(1,3,5), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_test.pdf", width = 6, height = 6)
# 设置图形参数以提高美观性
par(cex.axis = 1, cex.lab = 1.5, # 轴标签大小
    font.axis = 2, font.lab = 2, # 字体风格
    mar = c(5, 5, 2, 2) + 0.1, # 边距
    mgp = c(3, 1, 0), # 轴标题和轴线的距离
    tck = -0.02) # 轴刻度长度
# 假设ROC_rt是一个包含时间依赖性ROC数据的对象
# 绘制三条ROC曲线，分别对应不同时间点
plot(ROC_rt, time = 1, col = "#E64B35AA", main = "", lwd = 3, ylim = c(0.5, 1), xlim = c(0, 1))
plot(ROC_rt, time = 3, col = "#8491B4AA", add = TRUE, main = "", lwd = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(ROC_rt, time = 5, col = '#00A087B2', add = TRUE, main = "", lwd = 3)
# 添加图例
legend('bottomright', 
       legend = c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
                  paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
                  paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[3]))), 
       col = c("#E64B35AA", "#8491B4AA", '#00A087B2'), 
       lwd = 3, 
       bty = 'n', 
       text.font = 2,
       cex = 0.9) # 图例文字大小
# 关闭图形设备
dev.off()
##-----------------------测试Cindex----------------------------

###-----------------GSE20685--------------------
score_20685 <-score_test[c(1:327),]
cut <- surv_cutpoint(score_20685,'OS.time','OS','Score')
plot(cut)
cat <- surv_categorize(cut)
survifit <- survfit(Surv(OS.time,OS)~Score,cat)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~Score,cat)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

pdf("KM_20685.pdf",width = 6,height =6)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by =3.75,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()

col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(score_20685$OS.time,score_20685$OS,score_20685$Score,
                  cause = 1,weighting = 'marginal',
                  times=c(1,3,5), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_20685.pdf", width = 6, height = 6)
par(cex.axis = 1, cex.lab = 1.5, # 轴标签大小
    font.axis = 2, font.lab = 2, # 字体风格
    mar = c(5, 5, 2, 2) + 0.1, # 边距
    mgp = c(3, 1, 0), # 轴标题和轴线的距离
    tck = -0.02) # 轴刻度长度
# 假设ROC_rt是一个包含时间依赖性ROC数据的对象
# 绘制三条ROC曲线，分别对应不同时间点
plot(ROC_rt, time = 1, col = "#E64B35AA", main = "", lwd = 3, ylim = c(0.5, 1), xlim = c(0, 1))
plot(ROC_rt, time = 3, col = "#8491B4AA", add = TRUE, main = "", lwd = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(ROC_rt, time = 5, col = '#00A087B2', add = TRUE, main = "", lwd = 3)
legend('bottomright', 
       legend = c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
                  paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
                  paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[3]))), 
       col = c("#E64B35AA", "#8491B4AA", '#00A087B2'), 
       lwd = 3, 
       bty = 'n', 
       text.font = 2,
       cex = 0.9) # 图例文字大小
dev.off() 
#Cindex
library(pec)
score_20685$`age at diagnosis:ch1`<-as.numeric(score_20685$`age at diagnosis:ch1`)
score_20685$`event_metastasis:ch1`<-as.factor(score_20685$`event_metastasis:ch1`)
score_20685$`m_stage:ch1`<-as.factor(score_20685$`m_stage:ch1`)
met20685 <- coxph(Surv(OS.time, OS) ~ `event_metastasis:ch1`, data = score_20685, x = TRUE)
age20685 <- coxph(Surv(OS.time, OS) ~ `age at diagnosis:ch1`, data = score_20685, x = TRUE)
score20685<-coxph(Surv(OS.time, OS) ~ `Score`, data = score_20685, x = TRUE)
t20685<-coxph(Surv(OS.time, OS) ~ `t_stage:ch1`, data = score_20685, x = TRUE)
m20685<-coxph(Surv(OS.time, OS) ~ `m_stage:ch1`, data = score_20685, x = TRUE)
times <- seq(min(score_20685$OS.time), max(score_20685$OS.time), by = 1)
pk<- cindex(list(age20685, 
                 t20685, 
                 m20685,
                 score20685),
            Surv(OS.time, OS) ~ 1, 
            data = score_20685,
            eval.times=times)
# 保存图形到PDF文件
pdf("c-index20685.pdf", width = 6, height = 6)
# 调整绘图参数
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)
# 绘制五条 C 指数曲线
plot(pk, col = c("#E64B35AA", "#8491B4AA",  "#FFD54FAA",'#00A087B2'),
     lwd = 2, xlab = "Time (years)", ylab = "C-index", 
     ylim = c(0.4, 0.8),
     xlim = c(0.5, max(score_20685$OS.time)), main = "Time-dependent C-index for Models")
# 添加图例
legend('bottomright', legend = c("Age", "T Stage", "M Stage", "Score"), 
       col = c("#E64B35AA", "#8491B4AA", "#FFD54FAA",'#00A087B2'), 
       lwd = 3, bty = 'n', text.font = 2, cex = 0.9)
dev.off()

###-----------------GSE48390--------------------
score_48390 <-score_test[c(328:408),]
cut <- surv_cutpoint(score_48390,'OS.time','OS','Score')
plot(cut)
cat <- surv_categorize(cut)
survifit <- survfit(Surv(OS.time,OS)~Score,cat)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~Score,cat)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

pdf("KM_48390.pdf",width = 6,height =6)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by =1.5,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(score_48390$OS.time,score_48390$OS,score_48390$Score,
                  cause = 1,weighting = 'marginal',
                  times=c(1,3,5), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_48390.pdf", width = 6, height = 6)
par(cex.axis = 1, cex.lab = 1.5, # 轴标签大小
    font.axis = 2, font.lab = 2, # 字体风格
    mar = c(5, 5, 2, 2) + 0.1, # 边距
    mgp = c(3, 1, 0), # 轴标题和轴线的距离
    tck = -0.02) # 轴刻度长度
# 假设ROC_rt是一个包含时间依赖性ROC数据的对象
# 绘制三条ROC曲线，分别对应不同时间点
plot(ROC_rt, time = 1, col = "#E64B35AA", main = "", lwd = 3, ylim = c(0.5, 1), xlim = c(0, 1))
plot(ROC_rt, time = 3, col = "#8491B4AA", add = TRUE, main = "", lwd = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(ROC_rt, time = 5, col = '#00A087B2', add = TRUE, main = "", lwd = 3)
legend('bottomright', 
       legend = c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
                  paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
                  paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[3]))), 
       col = c("#E64B35AA", "#8491B4AA", '#00A087B2'), 
       lwd = 3, 
       bty = 'n', 
       text.font = 2,
       cex = 0.9) # 图例文字大小
dev.off()


###--------------------GSE88770-------------------
score_88770<-score_test[c(409:525),]
cut <- surv_cutpoint(score_88770,'OS.time','OS','Score')
plot(cut)
cat <- surv_categorize(cut)
survifit <- survfit(Surv(OS.time,OS)~Score,cat)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~Score,cat)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

pdf("KM_88770.pdf",width = 6,height =6)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by =4.5,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()

col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(score_88770$OS.time,score_88770$OS,score_88770$Score,
                  cause = 1,weighting = 'marginal',
                  times=c(3,5,10), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_88770.pdf", width = 6, height = 6)
par(cex.axis = 1, cex.lab = 1.5, # 轴标签大小
    font.axis = 2, font.lab = 2, # 字体风格
    mar = c(5, 5, 2, 2) + 0.1, # 边距
    mgp = c(3, 1, 0), # 轴标题和轴线的距离
    tck = -0.02) # 轴刻度长度
# 假设ROC_rt是一个包含时间依赖性ROC数据的对象
# 绘制三条ROC曲线，分别对应不同时间点
plot(ROC_rt, time = 3, col = "#E64B35AA", main = "", lwd = 3, ylim = c(0.5, 1), xlim = c(0, 1))
plot(ROC_rt, time = 5, col = "#8491B4AA", add = TRUE, main = "", lwd = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(ROC_rt, time = 10, col = '#00A087B2', add = TRUE, main = "", lwd = 3)
legend('bottomright', 
       legend = c(paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
                  paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
                  paste0('AUC at 10 years: ', sprintf("%.03f", ROC_rt$AUC[3]))), 
       col = c("#E64B35AA", "#8491B4AA", '#00A087B2'), 
       lwd = 3, 
       bty = 'n', 
       text.font = 2,
       cex = 0.9) # 图例文字大小
dev.off()
#Cindex
library(pec)
score_88770$`subtype:ch1`<-as.factor(score_88770$`subtype:ch1`)
type88770 <- coxph(Surv(OS.time, OS) ~ `subtype:ch1`, data = score_88770, x = TRUE)
score88770<-coxph(Surv(OS.time, OS) ~ Score, data = score_88770, x = TRUE)
times <- seq(min(score_88770$OS.time), max(score_88770$OS.time), by = 1)
pk<- cindex(list(
                 type88770,
                 score88770),
            Surv(OS.time, OS) ~ 1, 
            data = score_88770,
            eval.times=times)
library(stats)
# 进行 Mann-Whitney U 检验
wilcox_test <- wilcox.test(pk$AppCindex$coxph, pk$AppCindex$coxph1)
print(wilcox_test)  #p-value < 2.2e-16
# 保存图形到PDF文件
pdf("cindex88770.pdf", width = 6, height = 6) 
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)
plot(pk, col = c("#E64B35AA", '#00A087B2'),
     lwd = 2, xlab = "Time (years)", ylab = "C-index", 
     ylim = c(0.4, 0.8),
     xlim = c(3, max(score_88770$OS.time)), main = "Time-dependent C-index for Models")
# 添加图例
legend('bottomright', legend = c("Subtype","Score"), 
       col = c("#E64B35AA", '#00A087B2'), 
       lwd = 3, bty = 'n', text.font = 2, cex = 0.9)
dev.off()

###-------------------------GSE42568---------------------------
score_42568<-score_test[c(526:629),]
cut <- surv_cutpoint(score_42568,'OS.time','OS','Score')
plot(cut)
cat <- surv_categorize(cut)
survifit <- survfit(Surv(OS.time,OS)~Score,cat)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~Score,cat)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

pdf("KM_42568.pdf",width = 6,height =6)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by =2,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(score_42568$OS.time,score_42568$OS,score_42568$Score,
                  cause = 1,weighting = 'marginal',
                  times=c(1,3,5), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_42568.pdf", width = 6, height = 6)
par(cex.axis = 1, cex.lab = 1.5, # 轴标签大小
    font.axis = 2, font.lab = 2, # 字体风格
    mar = c(5, 5, 2, 2) + 0.1, # 边距
    mgp = c(3, 1, 0), # 轴标题和轴线的距离
    tck = -0.02) # 轴刻度长度
# 假设ROC_rt是一个包含时间依赖性ROC数据的对象
# 绘制三条ROC曲线，分别对应不同时间点
plot(ROC_rt, time = 1, col = "#E64B35AA", main = "", lwd = 3, ylim = c(0.5, 1), xlim = c(0, 1))
plot(ROC_rt, time = 3, col = "#8491B4AA", add = TRUE, main = "", lwd = 3, xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(ROC_rt, time = 5, col = '#00A087B2', add = TRUE, main = "", lwd = 3)
legend('bottomright', 
       legend = c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
                  paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
                  paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[3]))), 
       col = c("#E64B35AA", "#8491B4AA", '#00A087B2'), 
       lwd = 3, 
       bty = 'n', 
       text.font = 2,
       cex = 0.9) # 图例文字大小
dev.off()


#Cindex
score_42568$`size:ch1`<-as.numeric(score_42568$`size:ch1`)
size42568<-coxph(Surv(OS.time, OS) ~ `size:ch1`, data = score_42568, x = TRUE)
score_42568$`lymph node status:ch1`<-as.factor(score_42568$`lymph node status:ch1`)
lymph42568<-coxph(Surv(OS.time, OS) ~ `lymph node status:ch1`, data = score_42568, x = TRUE)
score_42568$`grade:ch1`<-as.factor(score_42568$`grade:ch1`)
grade42568<-coxph(Surv(OS.time, OS) ~ `grade:ch1`, data = score_42568, x = TRUE)
score_42568$`er_status:ch1`<-as.factor(score_42568$`er_status:ch1`)
er42568<-coxph(Surv(OS.time, OS) ~ `er_status:ch1`, data = score_42568, x = TRUE)
score_42568$`age:ch1`<-as.numeric(score_42568$`age:ch1`)
age42568<-coxph(Surv(OS.time, OS) ~ `age:ch1`, data = score_42568, x = TRUE)
score42568<-coxph(Surv(OS.time, OS) ~ Score, data = score_42568, x = TRUE)
times <- seq(min(score_42568$OS.time), max(score_42568$OS.time), by = 1)
pk<- cindex(list(size42568,lymph42568,grade42568,er42568,age42568,score42568),
  Surv(OS.time, OS) ~ 1, 
  data = score_42568,
  eval.times=times)
library(stats)
# 进行 Mann-Whitney U 检验
wilcox_test <- wilcox.test(pk$AppCindex$coxph, pk$AppCindex$coxph5)
print(wilcox_test)  #7.393e-15
# 保存图形到PDF文件
pdf("c-index42568.pdf", width = 6, height = 6) 
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)
plot(pk, col = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'black', '#00A087B2','#F4BBFF'),
     lwd = 2, xlab = "Time (years)", ylab = "C-index", 
     ylim = c(0.4, 0.8),
     xlim = c(1, max(score_42568$OS.time)), main = "Time-dependent C-index for Models")
# 添加图例
legend('bottomright', legend = c("Subtype","Score"), 
       col = c("#E64B35AA", "#8491B4AA", "#FFD54FAA",'black', '#00A087B2','#F4BBFF'), 
       lwd = 3, bty = 'n', text.font = 2, cex = 0.9)
dev.off()

