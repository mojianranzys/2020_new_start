##----------class8:基因突变数据可视化-------------##
rm(list = ls())
# 安装
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
#BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library(circlize)
setwd('D:/R/R_learn/2020/')
#突变信息
mutation <- './bioinformation_start/Bioinformatics_research_analysis/class8_mutation.txt'
mutation_data <- read.table(mutation,header = T,stringsAsFactors = F,sep="\t")
rownames(mutation_data) <- mutation_data[,1]

mat = as.matrix(mutation_data[,-1])

# 热图配色
col = c('Nonsense_Mutation'='#0071C0','Missense_Mutation'='red',
        'Nonstop_Mutation'='#4196FF','In_Frame_Del'='#92D14F',
        'In_Frame_Ins'='#92960A','Frame_Shift_Del'='#00B64E',
        'Frame_Shift_Ins'='#00784E','Silent'='#000000',
        'Splice_Site'='#6A3D9A','snv'="#3A76AE",'indel'="#83A5D6")

column_title = "OncoPrint for Cancers of Gene mutation types"

#alter_fun:单个函数或函数列表，用于定义如何为不同的更改添加图形。
alter_fun = function(x,y,w,h,v) {
  n = sum(v)
  h = h*0.9
  grid.rect(x, y, w*0.9, h*1, gp = gpar(fill = "gray90", col = NA))
  if(n){
    grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
              gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
  } 
}
#简单绘图
oncoPrint(mat, alter_fun = alter_fun,
          col = col,show_column_names = T)

# 修改尺寸
oncoPrint(mat,alter_fun = alter_fun,
          col = col,show_column_names = T,
          width = unit(8,'cm'),height = unit(8,"cm"))

# 清理无突变信息
oncoPrint(mat,alter_fun = alter_fun,
          col = col,show_column_names = T,
          remove_empty_rows = T,remove_empty_columns = T)

# 修改样本的顺序
oncoPrint(mat, alter_fun = alter_fun, 
          col = col, show_column_names = TRUE,
          column_order=c(sort(colnames(mat))))

# 图形分组
oncoPrint(mat,alter_fun = alter_fun,col = col,
          show_column_names = T,
          column_split = c("odd","even","odd","even","odd","even","odd","even","odd","even","odd","even","odd","even","odd"),
          column_title=NULL)

# 修改标签
oncoPrint(mat,alter_fun = alter_fun,
          col = col,show_column_names = T,
          show_pct = F,row_names_side = 'left',
          right_annotation = NULL)
# show_pct:是否在图的左侧显示百分比值？默认显示
# row_names_side:行显示在图的哪一侧，默认为右侧
# right_annotation：注释位于图的右侧。 默认情况下，它是条形图，显示每个基因中具有一定改变的样本数量。

# 修改图例
oncoPrint(mat,alter_fun = alter_fun,col = col,
          show_column_names = T,
          heatmap_legend_param = list(title = "Alternations",
                                      at = c('Missense_Mutation','Frame_Shift_Del','In_Frame_Del','Splice_Site','Frame_Shift_Ins','Nonsense_Mutation'),
                                      labels=c('Missense_Mutation', 'Frame_Shift_Del', 'In_Frame_Del', 'Splice_Site','Frame_Shift_Ins', 'Nonsense_Mutation')))

# 增加clinic 数据
clinic <- './bioinformation_start/Bioinformatics_research_analysis/class8_clinic.txt'
clinic_data <- read.table(clinic,header = T,stringsAsFactors = F)
rownames(clinic_data) <- clinic_data[,1]
group <- clinic_data[,-1]
# 注释
annot_df <- data.frame(gender = group$gender, history = group$history, TNM = group$TNM)
col2=list(gender=c("F"="grey","M"="pink"),
         history=c("yes"="tan2","no"="olivedrab3"),
         TNM = c("I"="deepskyblue","II"="gold2","III"="brown1","IV"="plum2"))
a1=HeatmapAnnotation(df = annot_df, col = col2,annotation_name_side = "left")
oncoPrint(mat, alter_fun = alter_fun, col = col, 
          show_column_names = TRUE,
          top_annotation = a1,
          column_title=NULL)
# 分组
oncoPrint(mat, alter_fun = alter_fun, col = col, 
          show_column_names = TRUE,
          column_split = group$history,
          top_annotation = a1,
          show_pct = F,row_names_side = 'left',
          column_title=NULL)



















