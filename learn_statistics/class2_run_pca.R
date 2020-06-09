rm(list = ls())
# 1.pca plot
library(corrplot)
library(FactoMineR)
library(factoextra)

data("iris")
iris_df <- iris[,-5]

# plot correlation
# corrplot:圆圈的大小和颜色深度表示相关性
corrplot(cor(iris_df))

# diag = F：去掉中间自身之间的correlation
corrplot(cor(iris_df),diag = F)  

corrplot(cor(iris_df),diag = F,type = 'upper')
corrplot(cor(iris_df),diag = F,type = 'lower')

# 2.Run PCA
iris.pca <- FactoMineR::PCA(X = iris_df,graph = F)

# 特征值：
iris.pca$eig

# 贡献率
fviz_eig(X = iris.pca,addlabels = T)
# 因为是正交的，所以所有贡献率之和为1

fviz_pca_ind(X = iris.pca,
             pointsize = 'coord',
             pointshape = 21,
             fill= iris$Species)
# 颜色
library(RColorBrewer)
RColorBrewer::display.brewer.all()
color_list <- RColorBrewer::brewer.pal(3,"Set1")

fviz_pca_ind(X = iris.pca,
             habillage = iris$Species,
             palette = color_list,
             addEllipses = T,
             label = 'none')

# 载荷图
fviz_pca_var(X = iris.pca,col.var = 'coord')
# 高权重的数值比较重要,即越接近+/-1，越重要


