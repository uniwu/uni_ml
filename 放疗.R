#放疗
#---meta---
#----meta----
rad<-metaphe[,c(1,9)]
rownames(rad)<-rad$PATIENT_ID
rad<-rad[,-1,drop=FALSE]
rownames(rad)<-gsub('-','.',rownames(rad))
comsam<-intersect(rownames(rad),rownames(score_tcag))
rad<-rad[comsam,,drop=F]
colnames(rad)[1]<-'radio'
rad$radio<-ifelse(rad$radio== "YES", 1, ifelse(rad$radio== "NO", 0,rad$radio))
score_meta2<-cbind(score_meta1,rad)
##----km----
survifit <- survfit(Surv(OS.time,OS)~radio,score_meta2)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~radio,score_meta2)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
##------km----
pdf("KM_meta_radio.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data = score_meta2, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "radio",
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
##---boxplot----
pdf("boxplot_radio20685.pdf", width = 6, height = 5.5)
ggplot(score_meta2, aes(x =radio, y =Score, color = radio)) +
  geom_boxplot(aes(fill =radio), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2" )) +  # 设置点的颜色
  scale_fill_manual(values = c("#8491B4B2","#E64B35B2" )) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_meta$Score),max(score_meta$Score))) +
  stat_compare_means(method = "t.test")  # 添加两组之间的p值


#----tcga-----
radio<-clinical.indexed[,c(2,62)]
radio<-as.data.frame(radio)
rownames(radio)<-radio$submitter_id
radio<-radio[,-1,drop=FALSE]
rownames(radio)<-gsub('-','.',rownames(radio))
comsam<-intersect(rownames(radio),rownames(score_train))
score_tcag1<-score_train[comsam,,drop=F]
radio<-radio[comsam,,drop=F]
score_tcag1<-cbind(score_tcag1,radio)
colnames(score_tcag1)[4]<-'radio'
score_tcag1$radio<-ifelse(score_tcag1$radio== "yes", 1, ifelse(score_tcag1$radio== "no", 0,score_tcag1$radio))
score_tcag1$radio<-as.factor(score_tcag1$radio)
score_tcag1 <- score_tcag1 %>%
  filter(radio != "not reported")
survifit <- survfit(Surv(OS.time,OS)~radio,score_tcag1)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~radio,score_tcag1)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
pdf("KM_tcga_radio.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data =score_tcag1, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "radio",
  legend.labs = c("0", "1"),
  xlab = "Time (years)",
  break.time.by =6,
  palette = c( "#8491B4AA","#E64B35AA"), # 修改透明度以增加美观性
  risk.table = TRUE,
  risk.table.title = "Number at risk",
  risk.table.height = .25,
  risk.table.y.text.col = T, # 让风险表的y轴文本颜色与生存曲线匹配
  risk.table.y.text = FALSE # 不显示y轴文本，保持图表简洁
)
dev.off()
##-----boxplot----
pdf("boxplot_tcga_radioall.pdf", width = 6, height = 5.5)
ggplot(score_tcag1, aes(x = radio, y =Score, color =radio)) +
  geom_boxplot(aes(fill =radio), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2" )) +  # 设置点的颜色
  scale_fill_manual(values = c("#8491B4B2","#E64B35B2")) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_tcag1$Score),max(score_tcag1$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
dev.off()
