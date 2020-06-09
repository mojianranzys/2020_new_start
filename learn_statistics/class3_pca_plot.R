rm(list = ls())

library(corrplot)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

# load GTEx
load("./bioinformation_start/learn_statistics/data/data/GTEx.RData")

# tissue 的注释文件
accession.df <- read.csv(file = './bioinformation_start/learn_statistics/data/data/GTEx_accession_table.csv')
annotation.df <- read.csv(file = "./bioinformation_start/learn_statistics/data/data/GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = '\t')

# make matrix 
GTEx.TPM.gene.mat.only_TPM <- log2(GTEx.TPM.gene.mat[,c(-1,-2)]+1)
rownames(GTEx.TPM.gene.mat.only_TPM) <- as.character(GTEx.TPM.gene.mat$Name)

# 计算标准差SD
SD.value <- apply(GTEx.TPM.gene.mat.only_TPM,1 , FUN = function(x){sd(x)})

# 选择变化最大的TOP500个基因
GTEx.TPM.gene.mat.only_TPM.sort <- GTEx.TPM.gene.mat.only_TPM[order(SD.value,decreasing = T),]
GTEx.TPM.gene.mat.only_TPM.select_500 <- GTEx.TPM.gene.mat.only_TPM.sort[1:500,]

# get tissue info
tissue_info <- as.character(annotation.df$SMTS[match(accession_id_vec,annotation.df$SAMPID)])

# show cor
# unique tissue
GTEx.TPM.unique_tissue <- GTEx.TPM.gene.mat.only_TPM.select_500[,!duplicated(tissue_info)]
colnames(GTEx.TPM.unique_tissue) <- tissue_info[!duplicated(tissue_info)]
corrplot(cor(GTEx.TPM.unique_tissue),type="lower",diag = F)

# all case
GTEx.TPM.all_case <- GTEx.TPM.gene.mat.only_TPM.select_500
colnames(GTEx.TPM.all_case) <- sprintf("%s.%d",tissue_info,1:length(tissue_info))
corrplot(cor(GTEx.TPM.all_case))

# run PCA
PCA_res = FactoMineR::PCA(t(GTEx.TPM.all_case),scale.unit = T,graph = F)
# 特征值：
PCA_res$eig

# var info
summary(PCA_res$var)

# plot rate of contribution
fviz_eig(X = PCA_res,addlabels = T)

# plot PCA1
fviz_contrib(X = PCA_res,choice = 'var',axes = 1)

# plot PCA2
fviz_contrib(X = PCA_res,choice = 'var',axes = 2)

# 载荷图
fviz_pca_var(X = PCA_res,col.var = 'cos2')

# PCA individual plot
# plot directly
fviz_pca_ind(X = PCA_res,pointsize = "cos2",pointshape = 21,fill=tissue_info)

fviz_pca_ind(PCA_res,pointshape=21, fill=tissue_info)

# plot with shadow
library(RColorBrewer)
RColorBrewer::display.brewer.all()
color_list <- colorRampPalette(brewer.pal(8,'Set2'))(length(unique(tissue_info)))

fviz_pca_ind(X = PCA_res,
             label="none", # hide individual labels
             habillage = factor(tissue_info), # color by groups
             palette = color_list,
             addEllipses = T) # Concentration ellipses













