library(readr)
library(data.table)  
library(dplyr)    
library(limma)
setwd("C:/Users/xinyiwu/Desktop/Ful数据/训练6/F7")
#----------TCGA-----------
TCGA_rawdata <- read_tsv("TCGA-BRCA.htseq_counts.tsv")
#注意这里是log(count+1)
ID <- fread("D:/CYP27A1/gencode.v22.annotation.gene.probeMap", header = T, sep = '\t', data.table = F)
TCGA_gset <- TCGA_rawdata %>%
  inner_join(ID, by = c("Ensembl_ID" = "id")) %>%
  dplyr::select(gene, starts_with("TCGA") )
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene) )
colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)
TCGA_group_list <- ifelse(as.numeric(substring(colnames(TCGA_gset),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list) # 1104Tumor  #113normal

#去除在一半以上样本中不表达的基因
table(rowSums(TCGA_gset) == 0)
fif <- dim(TCGA_gset)[2]/2
no_exp <- data.frame(count = apply(TCGA_gset, 1, function(x) length(which(x== 0))))
no_exp$gene <- row.names(no_exp)
no_exp <- no_exp[no_exp$count > fif, ]$gene
TCGA_gset <- TCGA_gset[- which(row.names(TCGA_gset) %in% no_exp), ]
save(TCGA_gset,file = 'TCGAbrca-count.RData')

#提取TCGAnormal
normal_samples <- colnames(TCGA_gset)[which(substr(colnames(TCGA_gset), 14, 14) == '1')]
normal_matrix <- TCGA_gset[, normal_samples]
tuomr<-TCGA_gset[, -which(colnames(TCGA_gset) %in% normal_samples)]
data<-cbind(tuomr,normal_matrix)
save(data,file = 'padd+gtex.RData')
#--------------limma---------
group_list <- c(rep('Tumor',1104),rep('Normal',113))
design <- model.matrix(~0+group_list)
contrast.matrix <- makeContrasts(group_listTumor - group_listNormal ,levels = design)
fit <- lmFit(data,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
allDiff <- topTable(fit,number = Inf)
allDiff$ID <- rownames(allDiff)
allDiff$group <- case_when(allDiff$logFC > 1 & allDiff$P.Value<0.05 ~ "Up",
                           allDiff$logFC < 1 & allDiff$P.Value < 0.05 ~ "Down",
                           abs(allDiff$logFC) <= 1 ~ "None",
                           allDiff$P.Value >= 0.05 ~ "None")
filtered_allDiff <- allDiff[allDiff$group != "None", ]
sorted_deg <- filtered_allDiff[order(filtered_allDiff$logFC), ]
# 提取数值最小的50行并命名为down，最大的50行并命名为up
down <- sorted_deg[1:100, ]
up <- sorted_deg[(nrow(sorted_deg) - 99):nrow(sorted_deg), ]

down$gene<-rownames(down)
down1<-down$gene
write.table(down1, file = "down1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

up$gene<-rownames(up)
up1<-up$gene
write.table(up1, file = "up1.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#------plot------------
cmap.res <- read.delim("C:/Users/xinyiwu/Desktop/Ful数据/训练6/F7/export (1).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
drug <- c("BRD9876","importazole","sorafenib","CD−1530","LBH−589","CR−1−31B")
selres <- cmap.res[which(cmap.res$Name %in% intersect(cmap.res$Name,drug)),]
print(selres) 

## 左侧区块
# 背景颜色
dt1 <- matrix(c(0,0, # 1为深色，0为浅色
                0,0,
                1,1,
                0,0,
                1,1,
                0,1),
              ncol = 2,
              byrow = T,
              dimnames = list(c("BRD9876","importazole","sorafenib","CD−1530","LBH−589","CR−1−31B"),
                              c("Clinical status","Experimental evidence")))

# 文字标签
lb1 <- matrix(c("Preclinical","Absent",
                "Preclinical","Absent",
                "Launched","Present",
                "Preclinical","Absent",
                "Launched","Present",
                "Preclinical","Present"),
              ncol = 2,
              byrow = T)


# 文字颜色
cl1 <- matrix(c("black","black",
                "black","black",
                "white","white",
                "black","black",
                "black","black",
                "black","white"),
              ncol = 2,
              byrow = T)

# 画图
# 如果这里报错，请返回开头，看“环境设置”
hm1 <- pheatmap(mat = dt1,
                color = c("#EFFAE8","#BAE5BC"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb1,
                number_color = cl1,
                fontsize_number = 11,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)



## 右侧区块
# 背景颜色
dt3 <- matrix(c(-1, # 0为中间色
                1, # -1为最浅色
                0,
                -1, # 1为深色
                -1,
                -1),
              ncol = 1,
              byrow = T,
              dimnames = list(c("BRD9876","importazole","sorafenib","CD−1530","LBH−589","CR−1−31B"),
                              c("CMap score")))

# 文字标签
# 根据selres的结果填写
# 因为有两个药物没有CMap结果，所以是NA
lb3 <- matrix(c("NA",
                "-90.11",
                "-83.43",
                "NA",
                "NA",
                "NA"),
              ncol = 1,
              byrow = T)

# 文字颜色
cl3 <- matrix(c("black",
                "white",
                "white",
                "black",
                "black",
                "black"),
              ncol = 1,
              byrow = T)

# 画图
hm3 <- pheatmap(mat = dt3,
                color = c("#F2F2F2",'#FDCEB9',"#F26E5F"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb3,
                number_color = cl3,
                fontsize_number = 11,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

## 水平合并热图
hm <- hm1  + hm3
draw(hm) 
dev.copy2pdf(file = "heatmap1.pdf",width = 3,height = 6)


#-------plot2--------
library(ComplexHeatmap) 
drug1 <- c("temsirolimus","tanespimycin","paclitaxel","sirolimus","everolimus")
selres1 <- cmap.res[which(cmap.res$Name %in% intersect(cmap.res$Name,drug1)),]
print(selres1) 

## 左侧区块
# 背景颜色
dt1 <- matrix(c(1,1, # 1为深色，0为浅色
                0,0,
                1,1,
                1,1,
                1,1),
              ncol = 2,
              byrow = T,
              dimnames = list(c("temsirolimus","tanespimycin","paclitaxel","sirolimus","everolimus"),
                              c("Clinical status","Experimental evidence")))

# 文字标签
lb1 <- matrix(c("Launched","Present",
                "Phase 2","Terminated",
                "Launched","Present",
                "Launched","Present",
                "Launched","Present"),
              ncol = 2,
              byrow = T)


# 文字颜色
cl1 <- matrix(c("white","white",
                "black","black",
                "white","white",
                "white","white",
                "white","white"),
              ncol = 2,
              byrow = T)

# 画图
# 如果这里报错，请返回开头，看“环境设置”
hm1 <- pheatmap(mat = dt1,
                color = c("#EFFAE8","#BAE5BC"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb1,
                number_color = cl1,
                fontsize_number = 11,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)



## 右侧区块
# 背景颜色
dt3 <- matrix(c(1, # 0为中间色
                -1, # -1为最浅色
                0,
                0, # 1为深色
                1),
              ncol = 1,
              byrow = T,
              dimnames = list(c("temsirolimus","tanespimycin","paclitaxel","sirolimus","everolimus"),
                              c("CMap score")))

# 文字标签
# 根据selres的结果填写
# 因为有两个药物没有CMap结果，所以是NA
lb3 <- matrix(c("-87.02",
                "NA",
                "37.76",
                "-63.86",
                "-85.7"),
              ncol = 1,
              byrow = T)

# 文字颜色
cl3 <- matrix(c("white",
                "black",
                "black",
                "black",
                "white"),
              ncol = 1,
              byrow = T)

# 画图
hm3 <- pheatmap(mat = dt3,
                color = c("#F2F2F2",'#FDCEB9',"#F26E5F"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb3,
                number_color = cl3,
                fontsize_number = 11,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

## 水平合并热图
hm <- hm1  + hm3
draw(hm) 
dev.copy2pdf(file = "heatmap2.pdf",width = 3,height = 6)
