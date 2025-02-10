#甜甜圈图
library(devtools)
install_github("Blake-Madden/donutchart")
library(donutchart)
library(tidyverse)
library(gridExtra)
#读取tcga和group文件
tcga<-tcga[,-(4:21)]
tcga$ID<-gsub('-','.',tcga$ID)
duough<-merge(tcga,group[,c(1:3)],by="ID")
#--------status-----------
##-----------high-----------
high<-duough[duough$group=="high",]
#汇总OS列数据
statushigh <- high %>% 
  group_by(OS) %>% 
  summarize(count = n())
p1 <- donut_chart(
  data = statushigh,
  groupColumn = OS,
  totalsColumn = count,
  centerLabel    = "Status",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = "#8bbbca",
  endColor = "#8bbbca",
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p1
##-----------low-----------
low<-duough[duough$group=="low",]
#汇总OS列数据
statuslow <- low %>% 
  group_by(OS) %>% 
  summarize(count = n())
p2 <- donut_chart(
  data = statuslow,
  groupColumn = OS,
  totalsColumn = count,
  centerLabel    = "Status",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = "#8bbbca",
  endColor = "#8bbbca",
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p2
#---pam50-----------
##-----------high-----------
pam50high <- high %>% 
  group_by(pam50) %>% 
  summarize(count = n())
p3 <- donut_chart(
  data = pam50high,
  groupColumn = pam50,
  totalsColumn = count,
  centerLabel    = "Pam50",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor =  "#E64B35AA",
  endColor = "#E64B35AA",
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p3
##-----low-----
pam50low <- low %>% 
  group_by(pam50) %>% 
  summarize(count = n())
p4 <- donut_chart(
  data = pam50low,
  groupColumn = pam50,
  totalsColumn = count,
  centerLabel    = "Pam50",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = "#E64B35AA",
  endColor = "#E64B35AA",
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p4
#--------T分期-------
##----------high--------
table(high$ajcc_pathologic_t)
high$ajcc_pathologic_t <- gsub("T1a|T1b|T1c", "T1", high$ajcc_pathologic_t)
high$ajcc_pathologic_t <- gsub("T2a|T2b", "T2", high$ajcc_pathologic_t)
high$ajcc_pathologic_t <- gsub("T4b|T4d", "T4", high$ajcc_pathologic_t)
# 查看修改后的统计
table(high$ajcc_pathologic_t)
t_high <- high %>% 
  group_by(ajcc_pathologic_t) %>% 
  summarize(count = n())
colnames(t_high)[1] <- c('T')
p5 <- donut_chart(
  data = t_high,
  groupColumn = T,
  totalsColumn = count,
  centerLabel    = "Pam50",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = '#00A087B2',
  endColor = '#00A087B2',
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p5
##------------low-----------
table(low$ajcc_pathologic_t)
low$ajcc_pathologic_t <- gsub("T1b|T1c", "T1", low$ajcc_pathologic_t)
low$ajcc_pathologic_t <- gsub('T3a' ,"T3", low$ajcc_pathologic_t)
low$ajcc_pathologic_t <- gsub("T4b", "T4", low$ajcc_pathologic_t)
table(low$ajcc_pathologic_t)
t_low <- low %>% 
  group_by(ajcc_pathologic_t) %>% 
  summarize(count = n())
colnames(t_low)[1] <- c('T')
p6<- donut_chart(
  data = t_low,
  groupColumn = T,
  totalsColumn = count,
  centerLabel    = "Pam50",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = '#00A087B2',
  endColor = '#00A087B2',
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p6
#--------N分期-------
##----------high--------
table(high$ajcc_pathologic_n)
high$ajcc_pathologic_n <- gsub("N0 \\(i-\\)|N0 \\(i\\+\\)", "N0", high$ajcc_pathologic_n)
high$ajcc_pathologic_n <- gsub("N1a|N1b|N1c|N1mi", "N1", high$ajcc_pathologic_n)
high$ajcc_pathologic_n <- gsub("N2a", "N2", high$ajcc_pathologic_n)
high$ajcc_pathologic_n <- gsub("N3a|N3b|N3c", "N3", high$ajcc_pathologic_n)
table(high$ajcc_pathologic_n)
n_high <- high %>% 
  group_by(ajcc_pathologic_n) %>% 
  summarize(count = n())
colnames(n_high)[1] <- c('N')
p7<- donut_chart(
  data = n_high,
  groupColumn = N,
  totalsColumn = count,
  centerLabel    = "Pam50",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = "#FFD54F",
  endColor = "#FFD54F",
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p7
##------------low-----------
table(low$ajcc_pathologic_n)
low$ajcc_pathologic_n <- gsub("N0 \\(i-\\)|N0 \\(i\\+\\)|N0 \\(mol\\+\\)", "N0", low$ajcc_pathologic_n)
low$ajcc_pathologic_n <- gsub("N1a|N1b|N1c|N1mi", "N1", low$ajcc_pathologic_n)
low$ajcc_pathologic_n <- gsub("N2a", "N2", low$ajcc_pathologic_n)
low$ajcc_pathologic_n <- gsub("N3a|N3b|N3c", "N3", low$ajcc_pathologic_n)
table(low$ajcc_pathologic_n)
n_low <- low %>% 
  group_by(ajcc_pathologic_n) %>% 
  summarize(count = n())
colnames(n_low)[1] <- c('N')
p8<- donut_chart(
  data = n_low,
  groupColumn = N,
  totalsColumn = count,
  centerLabel    = "N",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor = "#FFD54F",
  endColor = "#FFD54F",
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p8
#--------M分期-------
##----------high--------
table(high$ajcc_pathologic_m)
high$ajcc_pathologic_m <- gsub("cM0 \\(i\\+\\)", "M0", high$ajcc_pathologic_m)
table(high$ajcc_pathologic_m)
m_high <- high %>% 
  group_by(ajcc_pathologic_m) %>% 
  summarize(count = n())
colnames(m_high)[1] <- c('M')
p9<- donut_chart(
  data = m_high,
  groupColumn = M,
  totalsColumn = count,
  centerLabel    = "M",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor ='#4DBBD5',
  endColor = '#4DBBD5',
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p9
##------------low-----------
table(low$ajcc_pathologic_m)
low$ajcc_pathologic_m <- gsub("cM0 \\(i\\+\\)", "M0", low$ajcc_pathologic_m)
table(low$ajcc_pathologic_m)
m_low <- low %>% 
  group_by(ajcc_pathologic_m) %>% 
  summarize(count = n())
colnames(m_low)[1] <- c('M')
p10<- donut_chart(
  data = m_low,
  groupColumn = M,
  totalsColumn = count,
  centerLabel    = "M",
  centerColor   = "white",
  centerLabelColor = "black",
  centerLabelSize = 15,
  startColor ='#4DBBD5',
  endColor = '#4DBBD5',
  outerLabelColor = "black",
  outerLabelSize = 10,
)
p10
#---拼图----
grid.arrange(p1, p3, p5, p7, p9, p2, p4, p6, p8, p10, 
             nrow = 2, ncol = 5)



#------------highlowscore-------
all<-rbind(high,low)
# 使用 gsub 函数一次性替换
table(all$ajcc_pathologic_t)
all$ajcc_pathologic_t <- gsub("T1|T2", "T1-T2", all$ajcc_pathologic_t)
all$ajcc_pathologic_t <- gsub("T3|T4", "T3-T4", all$ajcc_pathologic_t)
allT <- subset(all, ajcc_pathologic_t != "TX")
ggplot(allT, aes(x = ajcc_pathologic_t, y = Score, color = ajcc_pathologic_t)) +
  geom_boxplot(aes(fill = ajcc_pathologic_t), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  
  scale_color_manual(values = c("#83B8D7","#E64B35")) +  
  scale_fill_manual(values = c("#83B8D7","#E64B35")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("T1-T2", "T3-T4")), 
              map_signif_level = TRUE, 
              textsize = 4, vjust = 0.5,
              color = "black")

table(all$ajcc_pathologic_n)
all$ajcc_pathologic_n <- gsub("N0|N1", "N0-N1", all$ajcc_pathologic_n)
all$ajcc_pathologic_n <- gsub("N2|N3", "N2-N3", all$ajcc_pathologic_n)
allN <- subset(all, ajcc_pathologic_n != "NX")
ggplot(allN, aes(x = ajcc_pathologic_n, y = Score, color = ajcc_pathologic_n)) +
  geom_boxplot(aes(fill = ajcc_pathologic_n), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  
  scale_color_manual(values = c("#83B8D7","#E64B35")) +  
  scale_fill_manual(values = c("#83B8D7","#E64B35")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("N0-N1", "N2-N3")), 
              map_signif_level = TRUE, 
              textsize = 4, vjust = 0.5,
              color = "black")


table(all$ajcc_pathologic_m)
allM <- subset(all, ajcc_pathologic_m != "MX")
ggplot(allM, aes(x = ajcc_pathologic_m, y = Score, color = ajcc_pathologic_m)) +
  geom_boxplot(aes(fill = ajcc_pathologic_m), alpha = 0.1) + 
  geom_jitter(width = 0.2) +  
  scale_color_manual(values = c("#83B8D7","#E64B35")) +  
  scale_fill_manual(values = c("#83B8D7","#E64B35")) +   
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_signif(comparisons = list(c("M0", "M1")), 
              map_signif_level = TRUE, 
              textsize = 4, vjust = 0.5,
              color = "black")
