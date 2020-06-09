## --------------------class6:R 语言绘图基础----------------#
#散点图 + 低级绘图函数
plot(c(1:26), c(1:14,13:2), main="plot(x,y)", sub="Type of Point and Line",
     xlab="point symbols", ylab="line ypes", ylim=c(0,15), 
     type="p", pch=0:25, col="darkblue", cex=1.5, cex.lab=1.2, font=2, bg="yellow")		
axis(1, c(1,3), c("A","C"))			# 自定义坐标轴 
axis(4, 5, "five")              # side：1下 2左 3上 4右
abline(h=1:6,lty=1:6,lwd=1:6,col="grey")	# 添加线条

# 柱状图
head(VADeaths)  # 查看数据集
barplot(VADeaths, beside = TRUE, 
        col = c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk"), 
        legend = rownames(VADeaths), ylim = c(0, 100))
barplot(VADeaths, col = c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk"),
        legend = rownames(VADeaths))

# red pocket
gender <- c("M","F","F","F","M","M","M","M","F","M","M","F")
age <- c(20,21,18,25,27,28,20,18,22,27,20,20)
money <- c(5,6,1,2,0.5,0.9,10,22,2,9,10,11)
result = cbind(gender, age, money)

# 饼图
table(result[,1])   # 统计性别
matrix(table(result[,1]))   # matrix转换
pie(c(matrix(table(result[1,1])),matrix(table(result[2,1]))),
    labels=c("girls","boys"),col=c("lightblue","lightgrey"))

# 直方图
hist(age,breaks=5)

# 箱线图
boxplot(money, col="grey")  # 单变量
boxplot(money~gender, col=topo.colors(2))   # 样本分组

# ggplur
library(ggpubr)
data <- cbind.data.frame(gender, age, money)  # 合并为数据框
p <- ggboxplot(data,x="gender",y="money",color="gender",add="jitter")
print(p)
p + stat_compare_means()  # method 设置其他检验方式

data("ToothGrowth")
p <- ggboxplot(ToothGrowth, x="dose", y="len", color = "dose", 
               palette = c("#00AFBB", "#E7B800", "#FC4E07"), add = "jitter", shape="dose")
print(p)
my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2")) # 指定比较组
p+stat_compare_means(comparisons = my_comparisons)

# 聚类图
hc <- hclust(dist(USArrests), "ave")  # dist()函数计算变量间距离, hclust()进行聚类
plot(hc, hang = -1) 
# hang 等于数值，表示标签与末端树杈之间的距离，若是负数，则表示末端树杈长度是0，即标签对齐。
plot(hc, hang = -1, cex = 0.5) 

# 热图
library(pheatmap)
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
head(test)
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# 作业
library(pheatmap)
DEG <- read.table("./class6_DEG.txt",sep="\t",header=T,row.names=1)
head(DEG)
pheatmap(DEG[2:7],show_rownames=T,scale="row",cex=0.6)



