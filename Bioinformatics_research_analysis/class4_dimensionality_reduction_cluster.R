##------------class4:降维分析&聚类分析-----------##
## 1. PCA
rm(list =ls())
# 1.1 模拟数据
data.matrix <- matrix(nrow =100,ncol =10)
#样本命名
colnames(data.matrix) <- c(
  paste("wt",1:5,sep=""),
  paste("ko",1:5,sep =""))
#基因命名
rownames(data.matrix) <- paste("gene",1:100,sep ="")
#模拟copy number
# rpois:从已知的平均发生率(rate)的泊松分布中生成泊松随机变量
# x=5:（非负整数）分位数的向量;lambda:（非负）均值的向量
# 分位数：quantile(x)
for (i in 1:100){
  wt.values <- rpois(5,lambda = sample(x=10:1000,size = 1))
  ko.values <- rpois(5,lambda = sample(x=10:1000,size = 1))
  data.matrix[i,] <- c(wt.values,ko.values)
}
#查看前5行,
head(data.matrix)
#1.2 run pca。prcomp函数要求样本为行，基因为列
pca<-prcomp(t(data.matrix),scale = T)
#prcomp函数返回3个值：x(cov(x)就是矩阵对角元素(sdev^2)),sdev(标准差),rotation(特征向量矩阵)

#1.3 作图
#使用前2列绘图
#x 包含绘制图表的主成分PCs,
plot(pca$x[,1],pca$x[,2])
#sdev返回标准差，查看每个PC解释的方差百分比
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per,main="Scree Plot",
        xlab = "Principa Component",
        ylab = "Percent Variation") 

#也可以用ggplot美化图片
library(ggplot2)
#按照ggplot2的格式要求转换数据
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data
ggplot(data=pca.data,aes(x=X,y=Y,label=Sample))+
  geom_text()+
  xlab(paste("PC1 - ", pca.var.per[1],"%",sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2],"%",sep=""))+
  theme_bw()+
  ggtitle("PCA Graph")

#1.4 查看loading score
#rotation返回Loading score
loading_scores <- pca$rotation[,1]
loading_scores
#由于对两组数据都感兴趣，所有用绝对值排序
gene_scores <- abs(loading_scores)
gene_scores_ranked <- sort(gene_scores,decreasing=T)
top_10_genes <- names(gene_scores_ranked[1:10])
top_10_genes

#看这些基因中哪一个是/负正的loading score
pca$rotation[top_10_genes,1]

## 2. K-means
#1.1 选择数据集
#使用著名的鸢尾花数据集，通过测定花萼等信息分为3类
newiris <- iris
#k-means 不支持因子类变量
newiris$Species <- NULL
kc <- kmeans(newiris, 3)
#2.2 比较聚类效果与实际数据
table(iris$Species, kc$cluster)

#2.3 作图
plot(newiris[c("Sepal.Length", "Sepal.Width")], col=kc$cluster)
points(kc$centers[,c("Sepal.Length", "Sepal.Width")], col=1:3, pch=8, cex=2)


###################################################
## 第三讲.作业,使用以下数据练习PCA分析，并绘碎石图。
###################################################
test<-data.frame(
  X1=c(148, 139, 160, 149, 159, 142, 153, 150, 151, 139,
       140, 161, 158, 140, 137, 152, 149, 145, 160, 156,
       151, 147, 157, 147, 157, 151, 144, 141, 139, 148),
  X2=c(41, 34, 49, 36, 45, 31, 43, 43, 42, 31,
       29, 47, 49, 33, 31, 35, 47, 35, 47, 44,
       42, 38, 39, 30, 48, 36, 36, 30, 32, 38),
  X3=c(72, 71, 77, 67, 80, 66, 76, 77, 77, 68,
       64, 78, 78, 67, 66, 73, 82, 70, 74, 78,
       73, 73, 68, 65, 80, 74, 68, 67, 68, 70),
  X4=c(78, 76, 86, 79, 86, 76, 83, 79, 80, 74,
       74, 84, 83, 77, 73, 79, 79, 77, 87, 85,
       82, 78, 80, 75, 88, 80, 76, 76, 73, 78)
)
pca1 <- prcomp(test,scale = T)
#x 包含绘制图表的主成分PCs,
plot(pca1$x[,1],pca1$x[,2])
#sdev返回标准差，查看每个PC解释的方差百分比
pca.var1 <- pca1$sdev^2
pca.var.per1 <- round(pca.var1/sum(pca.var1)*100,1)
# 碎石图
barplot(pca.var.per1,main="Scree Plot",
        xlab = "Principa Component",
        ylab = "Percent Variation") 


