#临床列线图构建
#--------TCGA单因素cox回归--------
tcga<-tcga[,-(4:21)]
rownames(tcga)<-tcga$ID
rownames(tcga)<-gsub("-",".",rownames(tcga))
tcga<-tcga[,-1]
comsam<-intersect(rownames(tcga),rownames(score_train))
score_tcga<-score_train[comsam,]
tcga<-tcga[comsam,]
tcgacli<-cbind(tcga,score_tcga$Score)
colnames(tcgacli)[11]<-'Score'
setwd('D:/Ful数据/训练6/F4')
pam50_mapping <- c("Basal" = 1, "Her2" = 2, "LumB" = 3, "LumA" = 4, "Normal" = 5)
# 使用映射关系替换 PAM50 列的值
tcgacli$pam50 <- pam50_mapping[as.character(tcgacli$pam50)]
t_mapping <- c(
  "T1" = 1,
  "T1a" = 2,
  "T1b" = 3,
  "T1c" = 4,
  "T2" = 5,
  "T2a" = 6,
  "T2b" = 7,
  "T3" = 8,
  "T3a" = 9,
  "T4" = 10,
  "T4b" = 11,
  "T4d" = 12,
  "TX" = 0)
# 使用映射关系替换 ajcc_pathologic_t 列的值
tcgacli$ajcc_pathologic_t <- t_mapping[as.character(tcgacli$ajcc_pathologic_t)]
n_mapping <- c(
  "N0" = 1,
  "N0 (i-)" = 2,
  "N0 (i+)" = 3,
  "N0 (mol+)" = 4,
  "N1" = 5,
  "N1a" = 6,
  "N1b" = 7,
  "N1c" = 8,
  "N1mi" = 9,
  "N2" = 10,
  "N2a" = 11,
  "N3" = 12,
  "N3a" = 13,
  "N3b" = 14,
  "N3c" = 15,
  "NX" = 0)
tcgacli$ajcc_pathologic_n <- n_mapping[as.character(tcgacli$ajcc_pathologic_n)]
m_mapping <- c(
  "cM0 (i+)" = 1,
  "M0" = 2,
  "M1" = 3,
  "MX" = 0)
tcgacli$ajcc_pathologic_m <- m_mapping[as.character(tcgacli$ajcc_pathologic_m)]
pharmaceutical_mapping <- c(
  "no" = 1,
  "not reported" = 0,
  "yes" = 2)
tcgacli$treatments_pharmaceutical_treatment_or_therapy <- pharmaceutical_mapping[as.character(tcgacli$treatments_pharmaceutical_treatment_or_therapy)]
radiation_mapping <- c(
  "no" = 1,
  "not reported" = 0,
  "yes" = 2)
tcgacli$treatments_radiation_treatment_or_therapy<- radiation_mapping[as.character(tcgacli$treatments_radiation_treatment_or_therapy)]
tcgacli<-tcgacli[,-4]
save(tcgacli,file = 'tcga临床单因素cox.RData')
##--------TCGA单因素cox回归--------
unresults <- data.frame(covariate = character(0),
                      beta = numeric(0),
                      HR=numeric(0),
                      HR.confint.lower=numeric(0),
                      HR.confint.upper=numeric(0),
                      wald.test = numeric(0),
                      p.value = numeric(0))
covariates<-colnames(tcgacli[,(3:10)])
for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = tcgacli,iter.max = 1000)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR<-signif(summary$conf.int[,1],3)
  HR.confint.lower <- signif(summary$conf.int[,"lower .95"], 3)
  HR.confint.upper <- signif(summary$conf.int[,"upper .95"], 3)
  unresults <- rbind(unresults, data.frame(covariate = covariate,
                                       beta = beta,
                                       HR=HR,
                                       HR.confint.lower=HR.confint.lower,
                                       HR.confint.upper=HR.confint.upper,
                                       wald.test= wald.test,
                                       p.value = p.value))
}
save(unresults,file = 'tcga单因素cox结果.RData')
##-----------plot-----------------
set.seed(23)
library(forestplot)
#按照p从小到大展示
unresults <- unresults[order(unresults$p.value), ]
#标注p<0.0001
unresults$p.value[unresults$p.value < 0.0001] <- "p<0.001"
#构建tabletext
unresults[2, 1] <- 'pharmaceutical_therapy'
unresults[3, 1] <- 'radiation_therapy'
unresults[4, 1] <- 'age'
unresults[5, 1] <- 'N'
unresults[6, 1] <- 'T'
unresults[7, 1] <- 'M'
sample <- as.data.frame(unresults)
tabletext1<-as.character(sample[,1])
tabletext2<-(sample[,'p.value'])
tabletext3<-paste(round(as.numeric(sample[,'HR']),3),round(as.numeric(sample[,'HR.confint.lower']),2),sep="(")
tabletext4<-paste(tabletext3,round(as.numeric(sample[,'HR.confint.upper']),2),sep="-")
tabletext5<-paste0(tabletext4,sep=")")
tabletext<-cbind(tabletext1,tabletext2,tabletext5)

forestplot(labeltext=tabletext, #文本信息  
           mean = round(sample[,'HR'],3),##HR值
           lower = round(sample[,"HR.confint.lower"],2),##95%置信区间
           upper = round(sample[,"HR.confint.upper"],2),#95%置信区间
           boxsize = 0.6,##大小
           graph.pos=4,#图在表中的列位置
           graphwidth = unit(0.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等
           col=fpColors(box="#8491B4AA", lines="black", zero = "black"),#颜色设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           grid=T,
           lwd.xaxis=2,#X轴线宽
           title="Hazard Ratio",
           xlab="",#X轴标题
           clip=c(-Inf,3),#边界
           colgap = unit(0.5,"cm")   
)

#筛选p<0.05纳入多因素cox回归
unresults<-subset(unresults, p.value < 0.05)
tcgacli<-tcgacli[,-3]
save(tcgacli,file = 'tcga临床多因素cox.RData')
#-----------多因素cox-------------
mulcox<- coxph(Surv(OS.time, OS) ~., data =tcgacli )
mulcoxSum=summary(mulcox)
outResult=data.frame()
outResult=cbind(
  HR=mulcoxSum$conf.int[,"exp(coef)"],
  L95CI=mulcoxSum$conf.int[,"lower .95"],
  H95CIH=mulcoxSum$conf.int[,"upper .95"],
  pvalue=mulcoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
##-----------plot----------
outResult<-as.data.frame(outResult)
outResult <- outResult[order(outResult$pvalue, na.last = NA), ]
save(outResult,file = 'tcga多因素cox回归结果.RData')
#构建tabletext
outResult[1, 1] <- 'T'
outResult[2, 1] <- 'N'
outResult[3, 1] <- 'age'
outResult[4, 1] <- 'pharmaceutical_therapy'
outResult[5, 1] <- 'M'
outResult[6, 1] <- 'radiation_therapy'
outResult$HR<-as.numeric(outResult$HR)
outResult$L95CI<-as.numeric(outResult$L95CI)
outResult$H95CIH<-as.numeric(outResult$H95CIH)
outResult$pvalue<-as.numeric(outResult$pvalue)
sample <- as.data.frame(outResult)
tabletext1<-as.character(sample[,1])
tabletext2 <- sprintf("%.1e", as.numeric(sample[,'pvalue']))
tabletext3<-paste(round(as.numeric(sample[,'HR']),3),round(as.numeric(sample[,'L95CI']),2),sep="(")
tabletext4<-paste(tabletext3,round(as.numeric(sample[,'H95CIH']),2),sep="-")
tabletext5<-paste0(tabletext4,sep=")")
tabletext<-cbind(tabletext1,tabletext2,tabletext5)

forestplot(labeltext=tabletext, #文本信息  
           mean = round(sample[,'HR'],3),##HR值
           lower = round(sample[,"L95CI"],2),##95%置信区间
           upper = round(sample[,"H95CIH"],2),#95%置信区间
           boxsize = 0.6,##大小
           graph.pos=4,#图在表中的列位置
           graphwidth = unit(0.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等
           col=fpColors(box="#8491B4AA", lines="black", zero = "black"),#颜色设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           grid=T,
           lwd.xaxis=2,#X轴线宽
           title="Hazard Ratio",
           xlab="",#X轴标题
           clip=c(-Inf,3),#边界
           colgap = unit(0.5,"cm")   
)
#只有score和T分期被纳入列线图

#-------------列线图--------------
library(rms)  # 调用制作列线图的包，rms包提供了专门用于回归建模和列线图创建的函数
nomotcga<-tcgacli[,c(1,2,6,11)]
nomotcga$ajcc_pathologic_t <- gsub("T1a|T1b|T1c", "T1", nomotcga$ajcc_pathologic_t)
nomotcga$ajcc_pathologic_t <- gsub("T2a|T2b", "T2", nomotcga$ajcc_pathologic_t)
nomotcga$ajcc_pathologic_t <- gsub("T3a", "T3", nomotcga$ajcc_pathologic_t)
nomotcga$ajcc_pathologic_t <- gsub("T4b|T4d", "T4", nomotcga$ajcc_pathologic_t)
nomotcga$ajcc_pathologic_t<-as.factor(nomotcga$ajcc_pathologic_t)
colnames(nomotcga)[4]<-'Score'
save(nomotcga,file = '列线图建立数据.RData')
dd <- datadist(nomotcga)
options(datadist = "dd")
set.seed(23)
coxfit <- cph(Surv(OS.time, OS) ~ Score+ajcc_pathologic_t,
              data = nomotcga, x=T,y=T,surv = T)
# 构建生存函数，注意你的最大生存时间
surv <- Survival(coxfit) 
surv1 <- function(x) surv(3,x) # 3年OS
surv2 <- function(x) surv(5,x) # 5年OS
surv3 <- function(x) surv(10,x) # 10年OS

nom <- nomogram(coxfit,
                fun = list(surv1,surv2,surv3),
                lp = T,
                funlabel = c('3-year survival Probability',
                             '5-year survival Probability',
                             '10-year survival Probability'),
                maxscale = 100,
                fun.at = c(0.95,0.8,0.6,0.4,0.2,0.1))
# 绘制 nomogram并缩小左侧标签字体大小
plot(nom,
     xfrac=.2,#图形与变量占比
     cex.var=1,#变量字体加粗
     cex.axis=0.8,#数轴 字体的大小
     lmgp=0.2,#文字与刻度的距离
     label.every=2,#隔一个显示一个
     #col.grid=gray(c(0.6, 0.95)),#设置垂直线的颜色
     col.grid=c('steelblue',alpha('steelblue',0.2)),
     lplabel='Linear Predictorlp',#线性预测轴名字
     points.label='Points',#变量分数名字
     total.points.label='Total Points',#总分名字
     force.label=T,#强制标记的每个刻度线都绘制标签
     title='Main' )     

##--------------------------------nomo校准曲线----------------------------------
library(rms)
library(calibrate)
nomoRisk=predict(coxfit, data=nomotcga, type="lp")
tcganomo=cbind(nomotcga, Nomogram=nomoRisk)
save(tcganomo,file = 'tcganomoRisk.RData')
set.seed(23)
#3年校准曲线
f3<-  cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=tcganomo, time.inc=3)
cal3 <- rms::calibrate(f3, cmethod='KM', method="boot", u=3, m=nrow(tcganomo)/3, B=1000)
#5年校准曲线(这里的单位本来就是年，不x365)
f5 <-  cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=tcganomo, time.inc=5)
#m要根据样本量来确定，由于标准曲线通常将所有样本分成3到4组
#而m代表每个分组的样本数，因此m*3应该等于或近似等于总的样本量
#使用rms::calibrate来确保调用的是rms包中的函数，这样可以避免命名空间的冲突
cal5 <- rms::calibrate(f5,cmethod='KM', method="boot", u=5, m=nrow(tcganomo)/3, B=1000)
#10年校准曲线
f10<-  cph(Surv(OS.time, OS) ~ Nomogram, x=T, y=T, surv=T, data=tcganomo, time.inc=10)
cal10 <- rms::calibrate(f10, cmethod='KM', method="boot", u=10, m=nrow(tcganomo)/3, B=1000)
plot(cal3, lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6,
     errbar.col = c("#E64B35AA"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#E64B35AA"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1.5, col = c("#E64B35AA"), pch = 16)
mtext("")
plot(cal5,lwd = 2,lty = 1,errbar.col = c("#8491B4AA"),
     xlim = c(0,1),ylim= c(0,1),col = c("#8491B4AA"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1.5, col = c("#8491B4AA"), pch = 16)
plot(cal10,lwd = 2,lty = 1,errbar.col = c('#00A087B2'),
     xlim = c(0,1),ylim= c(0,1),col = c('#00A087B2'),add = T)
lines(cal10[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1.5, col = c('#00A087B2'), pch = 16)
abline(0, 1, lty = 3, lwd = 2, col = c(rgb(0, 118, 192, maxColorValue = 255)))
legend("topleft", legend = c("3-years","5-years","10-years"), col =c("#E64B35AA","#8491B4AA",'#00A087B2'),
       lwd = 2,cex = 1.2,bty = "n")#不显示图例边框

##------------------------------列线图的决策曲线-----------------------------------
devtools::install_github('yikeshu0611/ggDCA')
library(ggDCA)
library(survival)
library(cowplot)
#nomo模型与单独因素相比的DCA曲线
cph1 <- coxph(Surv(OS.time,OS) ~ajcc_pathologic_t, tcganomo)
cph2 <- coxph(Surv(OS.time,OS) ~Nomogram, tcganomo)
ggplot(dca(cph1,cph2),linetype =T,lwd = 0.7)+  
  theme(legend.position=c(0.8,0.8))+
  xlim(0,0.8)
# 绘制 DCA 曲线，设置实线和自定义颜色
ggplot(dca(cph1, cph2), linetype = 1, lwd = 0.9) +  
  theme(legend.position = c(0.8, 0.8)) +
  xlim(0, 0.8) +
  scale_colour_manual(values = c("#8491B4AA","#E64B35AA","#FFD54FAA",'#00A087B2')) 
##------------------列线图的时间依赖性Cindex-------------------------------------
nomo<-coxph(Surv(OS.time, OS) ~ Nomogram, data = tcganomo, x = TRUE)
tnomo<-coxph(Surv(OS.time, OS) ~ ajcc_pathologic_t,data = tcganomo, x = TRUE)
times <- seq(min(tcganomo$OS.time), max(tcganomo$OS.time), by = 2)
pk<- cindex(list(nomo,tnomo),
            Surv(OS.time, OS) ~ 1, 
            data = tcganomo,
            eval.times=times)
# 创建pdf文件并设置绘图参数
pdf("c-indextcganomo.pdf", width = 6, height = 6) 
par(cex.axis = 1.2, cex.lab = 1.5, font.axis = 2, font.lab = 2, 
    mar = c(5, 5, 2, 2) + 0.1, mgp = c(3, 1, 0), tck = -0.02)
plot(pk, col = c("#E64B35AA","#8491B4AA" ),
     lwd = 2, xlab = "Time (years)", ylab = "C-index", 
     ylim = c(0.4, 0.8),
     xlim = c(0.1, max(tcganomo$OS.time)), main = "Time-dependent C-index for Models")
# 添加图例
legend('bottomright', legend = c('Nomogram','T'), 
       col =  c("#E64B35AA","#8491B4AA"),
       lwd = 3, bty = 'n', text.font = 2, cex = 0.9)
dev.off()

##----------列线图3、5、10年ROC------------------------------
col <- c("#E64B35AA", "#8491B4AA",'#00A087B2') ## 自定义颜色
ROC_rt <- timeROC(tcganomo$OS.time,tcganomo$OS,tcganomo$Nomogram,
                  cause = 1,weighting = 'marginal',
                  times=c(3,5,10), ROC=TRUE)#,iid = T
ROC_rt
# 保存图形到PDF文件
pdf("ROC_nomotcga.pdf", width = 6, height = 6)
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


