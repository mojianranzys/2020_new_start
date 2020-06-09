rm(list = ls())
#install.packages('fpc')
#-------------1.hcluster-------------#

load('./bioinformation_start/learn_statistics/data/data/GTEx.RData')

# tissue的注释文件
accession_df <- read.csv(file = './bioinformation_start/learn_statistics/data/data/GTEx_accession_table.csv')
annotaion.df <- read.csv(file = "./bioinformation_start/learn_statistics/data/data/GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

# make matrix
GTEx.TPM.gene.mat.only_TPM <- log2(GTEx.TPM.gene.mat[,c(-1,-2)]+1)
rownames(GTEx.TPM.gene.mat.only_TPM) <- as.character(GTEx.TPM.gene.mat$Name)

# get index
col_index <- limma::strsplit2(colnames(GTEx.TPM.gene.mat.only_TPM),split = "\\.")[,5]

# tissue info
tissue_info <- as.character(annotaion.df$SMTS[match(accession_id_vec,annotaion.df$SAMPID)])
new_colname <- sprintf("%s_%s",tissue_info,col_index)
colnames(GTEx.TPM.gene.mat.only_TPM) = new_colname 

# run hclust
library('fpc')
dis.mat = dist(t(GTEx.TPM.gene.mat.only_TPM),method = "euclidean")
png("./bioinformation_start/learn_statistics/class4_hclust_plot.png",width = 1200,height = 600)
plot(hclust(dis.mat))
dev.off()

# plot with cor
library(pheatmap)
png("./bioinformation_start/learn_statistics/class4_hclust_pheatmap.png",width = 1000,height = 1000)
pheatmap(cor(GTEx.TPM.gene.mat.only_TPM))
dev.off()

#--------------2.MDS-------------#
# run MDS
dist.mat = dist(t(GTEx.TPM.gene.mat.only_TPM), method = "minkowski", p=5)
#dist.mat = sqrt(2*(1 - abs(cor(GTEx.TPM.gene.mat.only_TPM,method = "kendall"))))

MDS_mat <- cmdscale(d = dis.mat,k = 2,eig = T)
MDS_xy.df <- data.frame(MDS_1 = MDS_mat$points[,1],
                        MDS_2 = MDS_mat$points[,2],
                        tissue = tissue_info)
head(MDS_xy.df)

library(ggplot2)
png("./bioinformation_start/learn_statistics/class4_mds_plot.png",width = 800,height = 800)
ggplot(data = MDS_xy.df,aes(x=MDS_1,y=MDS_2,color=factor(tissue))) +
  geom_point()+
  geom_text(aes(label = tissue))+
  theme_bw()
dev.off()

#----------------3.t-SNE-----------------#
library(Rtsne)

SD.vector <- apply(GTEx.TPM.gene.mat[,c(-1,-2)], 1, FUN = function(x){sd(x)})
GTEx.TPM.gene.mat.log2 <- log2(GTEx.TPM.gene.mat[,c(-1,-2)]+1)
GTEx.TPM.gene.mat.log2.top3000 <- GTEx.TPM.gene.mat.log2[order(SD.vector,decreasing = T),][1:3000,]

set.seed(20200608)
tSNE_res <- Rtsne(X = t(GTEx.TPM.gene.mat.log2.top3000),dims = 3,perplexity = 10,pca = T)

# plot region
tSNE_res.df <- as.data.frame(tSNE_res$Y)
colnames(tSNE_res.df) <- c('tSNE1','tSNE2','tSNE3')
tSNE_res.df$tissue <- tissue_info
library(ggplot2)
#  tSNE1 & tSNE2
ggplot(data = tSNE_res.df,aes(tSNE1,tSNE2,color = tissue))+
  geom_point()+
  geom_text(aes(label=tissue))+
  theme_bw()

#  tSNE1 & tSNE3
ggplot(data = tSNE_res.df,aes(tSNE1,tSNE3,color = tissue))+
  geom_point()+
  geom_text(aes(label=tissue))+
  theme_bw()

#  tSNE2 & tSNE3
ggplot(data = tSNE_res.df,aes(tSNE2,tSNE3,color = tissue))+
  geom_point()+
  geom_text(aes(label=tissue))+
  theme_bw()

# use plotly
library(plotly)
length(unique(tissue_info))
color_list <- c(RColorBrewer::brewer.pal(12, "Set3"),RColorBrewer::brewer.pal(6, "Set2"))
p <- plot_ly(data = tSNE_res.df,x=~tSNE1,y = ~tSNE2,z=~tSNE3,color = ~tissue,colors = color_list) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'tSNE1'),
                      yaxis = list(title = 'tSNE2'),
                      zaxis = list(title = 'tSNE3')))

p
setwd("./bioinformation_start/learn_statistics/")
htmlwidgets::saveWidget(as_widget(p), "class4_tSNE_plotly.html")
setwd('../../../2020/')

#---------------4.SVD--------------------#
# run SVD
GTEx.mat <- as.matrix(GTEx.TPM.gene.mat.only_TPM)
dim(GTEx.mat)
svd.res <- svd(GTEx.mat)

U <- as.matrix(svd.res$u)
D <- diag(svd.res$d)
V <- as.matrix(svd.res$v)
GTEx.svd <- U %*% D %*% t(V)
dim(GTEx.svd)
plot(x=as.vector(GTEx.mat)[1:1000], y=as.vector(GTEx.svd)[1:1000])

# 压缩
i= 20
svd.mat = svd.res$u[,1:i] %*% diag(svd.res$d[1:i]) %*% t(svd.res$v[,1:i])

plot(x=as.vector(GTEx.mat)[1:1000], y=as.vector(svd.mat)[1:1000])




