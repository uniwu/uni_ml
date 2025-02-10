setwd('C:/Users/xinyiwu/Desktop/Ful数据/训练6/突变/maf_data')
library(maftools)
library(dplyr)
#获取所有的压缩文件，这一步可能会需要一些时间，取决于我们研究的癌症样本规模大小
files <- list.files(pattern = '*.gz',recursive = TRUE)
all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}
#数据整理
all_mut <- read.maf(all_mut)

a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))
gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))
mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))
for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}
#所有样本突变情况汇总/排序
gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
# 大家可以修改数字20，代表TOP多少，也可选择自己感兴趣的
gene_top <- gene_count$gene[1:20] ##保存
save(mat,mat_0_1,file = "TMB.rda") ##保存为RData
write.csv(mat,"all_mut_type.csv")
write.csv(mat_0_1,"all_mut_01.csv")
#----直接绘制-------
library(TCGAbiolinks)
library(maftools)
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf <- GDCprepare(query)

#替换样本名称
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode,1,12)
maf$Tumor_Sample_Barcode <- gsub("-",".",maf$Tumor_Sample_Barcode)
#提取高低风险组
maf.High <- maf[(maf$Tumor_Sample_Barcode %in% group$ID[group$group=="high"]),]
maf.Low <- maf[(maf$Tumor_Sample_Barcode %in% group$ID[group$group=="low"]),]

library(maftools)
mut.High <- read.maf(maf=maf.High,isTCGA=T)## 读取高风险亚型的突变数据
mut.Low <- read.maf(maf = maf.Low,isTCGA = T)## 读取低风险亚型的突变数据
mut.all <- read.maf(maf = maf,isTCGA = T)## 读取总的样本突变数据
#设置颜色信息
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')
#人种
racecolors = RColorBrewer::brewer.pal(n = 4,name = 'Spectral')
names(racecolors) = c("ASIAN", "WHITE", "BLACK_OR_AFRICAN_AMERICAN",  "AMERICAN_INDIAN_OR_ALASKA_NATIVE")
##-----plot------
oncoplot(maf = mut.High,
         colors = col,#给突变配色
         top = 20)
oncoplot(maf = mut.Low,
         colors = col,#给突变配色
         top = 20)


#------ComplexHeatmap包绘制--------
#列为样本行为基因
#为不同的突变类型设置不同的颜色
library(ComplexHeatmap)
library(scales)
col <- c("Nonsense_Mutation" = "#E377C2FF", "Missense_Mutation" = "#0869a1", 
         "Splice_Site" = "#008000","Multi_Hit"="black",
         "Frame_Shift_Ins"="#8C564BFF","Frame_Shift_Del"="#aa0c0b",
         "In_Frame_Del"="#ed7620", "In_Frame_Ins"="#bbe165",
         "Nonstop_Mutation"="darkblue","Translation_Start_Site"="#f1c055")
show_col(col)
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  }, 
  In_Frame_Ins=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  Multi_Hit=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  Nonstop_Mutation=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
  },
  Translation_Start_Site=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  }
)
column_title <- "OncoPrint for High Risk"
#heatmap_legend_param 
heatmap_legend_param <- list(title = "Alternations", at = c("Nonsense_Mutation" , "Missense_Mutation", 
                                                            "Splice_Site" ,"Multi_Hit",
                                                            "Frame_Shift_Del","Frame_Shift_Ins",
                                                            "In_Frame_Del", "In_Frame_Ins",
                                                            "Nonstop_Mutation","Translation_Start_Site"), 
                             labels = c("Nonsense_Mutation" , "Missense_Mutation", 
                                        "Splice_Site" ,"Multi_Hit",
                                        "Frame_Shift_Del","Frame_Shift_Ins",
                                        "In_Frame_Del", "In_Frame_Ins",
                                        "Nonstop_Mutation","Translation_Start_Site"))
library(magick)
colnames(mat)<-gsub('-','.',colnames(mat))
rownames(group)<-group$ID
high<-group[group$group == "high", ]
comsam<-intersect(rownames(high),colnames(mat))#478
mathigh<-mat[,comsam]
oncoPrint(mathigh,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          pct_side = "right", row_names_side = "left",
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          alter_fun_is_vectorized = TRUE)

low<-group[group$group == "low", ]
comsam<-intersect(rownames(low),colnames(mat))#478
matlow<-mat[,comsam]
oncoPrint(matlow,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          pct_side = "right", row_names_side = "left",
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          alter_fun_is_vectorized = FALSE)
