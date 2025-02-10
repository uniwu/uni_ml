#内分泌治疗
#----meta----
hor<-metaphe[,c(1,9)]
rownames(hor)<-hor$PATIENT_ID
hor<-hor[,-1,drop=FALSE]
rownames(hor)<-gsub('-','.',rownames(hor))
comsam<-intersect(rownames(hor),rownames(score_train))
hor<-hor[comsam,,drop=F]
colnames(hor)[1]<-'hormone'
hor$hormone<-ifelse(hor$hormone== "YES", 1, ifelse(hor$hormone== "NO", 0,hor$hormone))
score_meta1<-cbind(score_meta,hor)
##----km----
survifit <- survfit(Surv(OS.time,OS)~hormone,score_meta1)#拟合生存曲线模型
diff=survdiff(Surv(OS.time,OS)~hormone,score_meta1)#计算差异
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
ggsurvplot(
  survifit, data = score_meta1, # 确保survifit和cat是你的survfit对象和分类变量的实际名称
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
ggplot(score_meta1, aes(x =hormone, y =Score, color = hormone)) +
  geom_boxplot(aes(fill =hormone), alpha = 0.1) +
  geom_jitter(width = 0.2) +  # 调整 width 参数以缩小点的大小
  scale_color_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +  # 设置点的颜色
  scale_fill_manual(values = c("#E64B35AA","#8491B4AA", "#FFD54FAA",'#F39B7F','#4DBBD5')) +   # 设置箱线图的填充颜色
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(min(score_meta$Score),max(score_meta$Score))) +
  stat_compare_means(method = "wilcox.test")  # 添加两组之间的p值
