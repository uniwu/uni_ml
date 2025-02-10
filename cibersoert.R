#cibersort
res_ciber<-res_ciber[(1:22),]
rownames(group)<-group$ID
comsam<-intersect(rownames(group),colnames(res_ciber))
ciber<-res_ciber[,comsam]
ciber<-rbind(ciber,group[2,])
##筛选出low组
high <- ciber %>%
  select(which(apply(ciber == "high", 2, any)))
high$celltype<-rownames(high)
high<-high[-23,]
high1<-melt(high,id.vars='celltype')
high1$value<-as.numeric(high1$value)
colour =  c(brewer.pal(12, "Paired"),brewer.pal(8, "Dark2"),brewer.pal(12, "Set3"))
pdf("ciber_high.pdf",width=8,height=6)
ggplot(data=high1,aes(x=variable,y=value,fill=celltype))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_manual(values = colour)+
  labs(x="",y="",title="cell proportion")+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8), 
        panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title=element_text(hjust=0.8))
dev.off()
##筛选出low组
low <- ciber %>%
  select(which(apply(ciber == "low", 2, any)))
low$celltype<-rownames(low)
low<-low[-23,]
low1<-melt(low,id.vars='celltype')
low1$value<-as.numeric(low1$value)
colour =  c(brewer.pal(12, "Paired"),brewer.pal(8, "Dark2"),brewer.pal(12, "Set3"))
pdf("ciber_low.pdf",width=8,height=6)
ggplot(data=low1,aes(x=variable,y=value,fill=celltype))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_manual(values = colour)+
  labs(x="",y="",title="cell proportion")+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(ncol = 1))+
  theme_bw()+
  theme(legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8), 
        panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title=element_text(hjust=0.8))
dev.off()
