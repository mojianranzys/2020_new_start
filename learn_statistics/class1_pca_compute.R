# make pca
rm(list = ls())

# load data
data(iris)
head(iris)
# 1.R包计算
# load package
library(FactoMineR)
iris.mat <- iris[,1:4]
iris.pca <- FactoMineR::PCA(iris.mat,graph = F)

# 特征值
iris.pca$eig

# PCA
iris.pca$ind$coord[,1]
iris.pca$ind$coord[,2]
# 2.手动运算
# 归一化
iris.mat.scale <- scale(iris.mat,center = T)
# 协方差矩阵
iris.mat.cov <- cov(iris.mat.scale)
# 求解特征值
iris.mat.cov.eigen <- eigen(iris.mat.cov)
iris.mat.cov.eigen$values
# 特征向量
iris.mat.cov.eigen$vectors

# 求pc1,pc2等
pc1 <- iris.mat.scale %*% iris.mat.cov.eigen$vectors[,1]
pc2 <- iris.mat.scale %*% iris.mat.cov.eigen$vectors[,2]

# 3.判断包计算出来PCA结果和手动计算出来的结果是否一致
plot(x = pc1,y = iris.pca$ind$coord[,1])
plot(x = pc2,y = iris.pca$ind$coord[,2])
# 看图是否为一条直线，如果是直线，则结果一样
