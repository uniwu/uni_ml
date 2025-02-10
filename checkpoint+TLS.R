#-------------免疫检查点------------
TCGAtumor$gene<-rownames(TCGAtumor)
cpexp<-merge(TCGAtumor,免疫检查点基因集,by="gene")
cpexp<-t(cpexp)
colnames(cpexp)<-cpexp[1,]
cpexp<-cpexp[-1,]
rownames(cpexp)<-gsub('-','.',rownames(cpexp))
cpexp<-as.data.frame(cpexp)
cpexp$ID<-rownames(cpexp)
cpexp1<-merge(cpexp,group,by="ID")
library(circlize)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(grImport2)
library(gridBase)
library(dplyr)
library(data.table)
##------cor计算----------------
rownames(cpexp1)<-cpexp1$ID
cpexp2<-cpexp1[,-1]
cpexp2<-cpexp2[,-70]
cpexp2<-cpexp2[,-68]
# 将所有变量转换为数值型
cpexp2 <- mutate_all(cpexp2, as.numeric)
corResult=data.frame()
for(i in colnames(cpexp2)){
  tdc <- cor.test(cpexp2[,"Score"], cpexp2[,i], method = "spearman")
  pvalue=tdc$p.value
  corResult=rbind(corResult,
                  cbind(id=i,
                        R= tdc$estimate,
                        Low95CI=tdc$conf.int[1],
                        High95CI=tdc$conf.int[2],
                        pvalue=tdc$p.value)
  )
}
#筛选p<0.05的基因
significant_genes <- corResult[corResult$pvalue < 0.05, ]
significant_genes$R<-as.numeric(significant_genes$R)
significant_genes <- significant_genes[order(abs(significant_genes$R), decreasing = TRUE), ]
#筛选相关性top8的基因
significant_genes<-significant_genes[1:9,]
rownames(significant_genes)<-significant_genes$id
comsam<-intersect(rownames(significant_genes),colnames(cpexp2))
cpexp3<-cpexp2[,comsam]
gene<-colnames(cpexp3[,2:9])
cpexp3<-t(cpexp3)
##-------------plot-------------------
plotlst <- list()
for (i in gene) {
  x=as.numeric(cpexp3["Score",])
  y=as.numeric(cpexp3[i,])
  df2=as.data.frame(cbind(x,y))
    p=ggplot(df2, aes(x, y)) + 
      xlab('Score')+ylab(i)+
      geom_point(color = "#8491B4B2") + # 自定义散点颜色
      geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
      theme_bw() +
      stat_cor(method = 'spearman', aes(x = x, y = y))
    plotlst[[i]] <- p
}
# 拼图，4个一行
allplot <- plot_grid(plotlist = plotlst, ncol = 4, align = "hv")
allplot


#-----------三级淋巴结构相关基因------------
#参考https://pubmed.ncbi.nlm.nih.gov/31092904/
#https://www.nature.com/articles/s41586-019-1914-8#Sec1
#结合两篇文献得到50个TLS基因
tlsexp<-merge(TCGAtumor,TLS基因,by="gene")
tlsexp<-t(tlsexp)
colnames(tlsexp)<-tlsexp[1,]
tlsexp<-tlsexp[-1,]
rownames(tlsexp)<-gsub('-','.',rownames(tlsexp))
tlsexp<-as.data.frame(tlsexp)
tlsexp$ID<-rownames(tlsexp)
tlsexp1<-merge(tlsexp,group,by="ID")
##------cor计算----------------
rownames(tlsexp1)<-tlsexp1$sample
tlsexp2<-tlsexp1[,-1]
tlsexp2<-tlsexp2[,-52]
tlsexp2<-tlsexp2[,-50]
# 将所有变量转换为数值型
tlsexp2 <- mutate_all(tlsexp2, as.numeric)
corResult1=data.frame()
for(i in colnames(tlsexp2)){
  tdc <- cor.test(tlsexp2[,"Score"], tlsexp2[,i], method = "spearman")
  pvalue=tdc$p.value
  corResult1=rbind(corResult1,
                  cbind(id=i,
                        R= tdc$estimate,
                        Low95CI=tdc$conf.int[1],
                        High95CI=tdc$conf.int[2],
                        pvalue=tdc$p.value)
  )
}
#筛选p<0.05的基因
significant_genes1 <- corResult1[corResult1$pvalue < 0.05, ]
significant_genes1$R<-as.numeric(significant_genes1$R)
significant_genes1 <- significant_genes1[order(abs(significant_genes1$R), decreasing = TRUE), ]
#筛选相关性top8的基因
significant_genes1<-significant_genes1[1:9,]
rownames(significant_genes1)<-significant_genes1$id
comsam<-intersect(rownames(significant_genes1),colnames(tlsexp2))
tlsexp3<-tlsexp2[,comsam]
gene<-colnames(tlsexp3[,2:9])
tlsexp3<-t(tlsexp3)
##-------------plot-------------------
plotlst <- list()
for (i in gene) {
  x=as.numeric(tlsexp3["Score",])
  y=as.numeric(tlsexp3[i,])
  df2=as.data.frame(cbind(x,y))
  p=ggplot(df2, aes(x, y)) + 
    xlab('Score')+ylab(i)+
    geom_point(color = "#8491B4B2") + # 自定义散点颜色
    geom_smooth(method = "lm", formula = y ~ x, color = "#E64B35") + # 自定义拟合线颜色
    theme_bw() +
    stat_cor(method = 'spearman', aes(x = x, y = y))
  plotlst[[i]] <- p
}
# 拼图，4个一行
allplot <- plot_grid(plotlist = plotlst, ncol = 4, align = "hv")
allplot
