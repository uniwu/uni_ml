#--------化疗-------
##-----meta-------
chemo<-metaphe[,c(1,5)]
rownames(chemo)<-chemo$PATIENT_ID
chemo<-chemo[,-1,drop=FALSE]
rownames(chemo)<-gsub('-','.',rownames(chemo))
comsam<-intersect(rownames(chemo),rownames(score_train))
chemo<-chemo[comsam,,drop=F]                      
chemo$CHEMOTHERAPY<-ifelse(chemo$CHEMOTHERAPY== "YES", 1, ifelse(chemo$CHEMOTHERAPY== "NO", 0,chemo$CHEMOTHERAPY))
score_meta<-score_train[comsam,,drop=F]
score_meta<-cbind(score_meta,chemo)
colnames(score_meta)[4]<-'chemotherapy'

###----km----
survifit <- survfit(Surv(OS.time,OS)~chemotherapy,score_meta)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~chemotherapy,score_meta)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
pdf("KM_meta_chemo.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data = score_meta, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
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
###-----boxplot-----
pdf("boxplot_meta_chemoall.pdf", width = 6, height = 5.5)
ggplot(score_meta, aes(x =chemotherapy, y =Score, color = chemotherapy)) +
  geom_boxplot(aes(fill =chemotherapy), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_meta$Score),max(score_meta$Score))) +
  stat_compare_means(method = "t.test")  # 添加两组之间的p值
dev.off()


#-------tcga-----
chemo<-clinical.indexed[,c(2,50)]
chemo<-as.data.frame(chemo)
rownames(chemo)<-chemo$submitter_id
chemo<-chemo[,-1,drop=FALSE]
rownames(chemo)<-gsub('-','.',rownames(chemo))
comsam<-intersect(rownames(chemo),rownames(score_train))
score_tcag<-score_train[comsam,,drop=F]
score_tcag<-cbind(score_tcag,chemo)
chemo<-chemo[comsam,,drop=F]
score_tcag<-cbind(score_tcag,chemo)
score_tcag$treatments_pharmaceutical_treatment_or_therapy<-ifelse(score_tcag$treatments_pharmaceutical_treatment_or_therapy== "yes", 1, ifelse(score_tcag$treatments_pharmaceutical_treatment_or_therapy== "no", 0,score_tcag$treatments_pharmaceutical_treatment_or_therapy))
colnames(score_tcag)[4]<-'pharmaceutical'
score_tcag$pharmaceutical<-as.factor(score_tcag$pharmaceutical)
score_tcag <- score_tcag %>%
  filter(pharmaceutical != "not reported")
##----km----
survifit <- survfit(Surv(OS.time,OS)~pharmaceutical,score_tcag)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~pharmaceutical,score_tcag)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
pdf("KM_tcga_pharm.pdf", width = 6, height = 6)
ggsurvplot(
  survifit, data =score_tcag, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
  conf.int = TRUE,
  pval = pValue, # 使用实际的P值
  pval.size = 6,
  legend.title = "pharmaceutical",
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
pdf("boxplot_tcga_pharmaceuticalall.pdf", width = 6, height = 5.5)
ggplot(score_tcag, aes(x =pharmaceutical, y =Score, color = pharmaceutical)) +
  geom_boxplot(aes(fill =pharmaceutical), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#8491B4B2","#E64B35B2" )) +  # 设置点的颜色
  scale_fill_manual(values =  c("#8491B4B2","#E64B35B2" )) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_meta$Score),max(score_meta$Score))) +
  stat_compare_means(method = "t.test")  # 添加两组之间的p值
dev.off()
rownames(tcga_pam50)<-gsub('-','.',rownames(tcga_pam50))
comsam<-intersect(rownames(tcga_pam50),rownames(score_tcag))
tcga_pam50<-tcga_pam50[comsam,,drop=F]
score_tcag<-cbind(score_tcag,tcga_pam50)
colnames(score_tcag)[5]<-'subtype'
ggplot(score_tcag, aes(x =subtype, y =Score, color = pharmaceutical)) +
  geom_boxplot(aes(fill =pharmaceutical), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_meta$Score),max(score_meta$Score))) +
  stat_compare_means(method = "t.test")  # 添加两组之间的p值
dev.off()