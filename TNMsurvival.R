t1t2<-allT[allT$ajcc_pathologic_t=='T1-T2',]
cut <- surv_cutpoint(t1t2,'OS.time','OS','Score')
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
pdf("t1t2KM.pdf",width = 5,height =5)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by = 6,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()


#-------t3t4---
t3t4<-allT[allT$ajcc_pathologic_t=='T3-T4',]
cut <- surv_cutpoint(t3t4,'OS.time','OS','Score')
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
pdf("t3t4KM.pdf",width = 5,height =5)
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

#------n0n1-----------
n0n1<-allT[allN$ajcc_pathologic_n=='N0-N1',]
cut <- surv_cutpoint(n0n1,'OS.time','OS','Score')
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
pdf("n0n1KM.pdf",width = 5,height =5)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by = 6,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()


#------n2n3-------
n2n3<-allT[allN$ajcc_pathologic_n=='N2-N3',]
cut <- surv_cutpoint(n2n3,'OS.time','OS','Score')
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
pdf("n2n3KM.pdf",width = 5,height =5)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by = 6,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()

#-------m0--------
m0<-allT[allM$ajcc_pathologic_m=='M0',]
cut <- surv_cutpoint(m0,'OS.time','OS','Score')
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
pdf("m0KM.pdf",width = 5,height =5)
ggsurvplot(
  survifit, data = cat, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "Risk",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (years)",
  break.time.by = 6,
  palette = c("#E64B35AA", "#8491B4AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()


#------m1--------
m1<-allT[allM$ajcc_pathologic_m=='M1',]
cut <- surv_cutpoint(m1,'OS.time','OS','Score')
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
pdf("m1KM.pdf",width = 5,height =5)
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
