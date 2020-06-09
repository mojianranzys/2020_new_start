## -------------class5:传统机器学习算法--------##
## 1. SVM
rm(list=ls())
library(e1071) #实现SVM
# iris；著名的（费舍尔或安德森氏）虹膜数据集以厘米为单位,
#分别测量了3种虹膜中每种花的50朵花的萼片长度和宽度以及花瓣的长度和宽度
data(iris)
plot(iris)
plot(iris$Sepal.Length,iris$Sepal.Width,col=iris$Species)
plot(iris$Petal.Length,iris$Petal.Width,col=iris$Species)

# 1.1 准备数据
#随机抽取100样本
s<-sample(150,100)
col <- c("Petal.Length","Petal.Width","Species")
#生成训练集和测试集
iris_train <- iris[s,col]
iris_test <- iris[-s,col]

#1.2模型训练
svmfit <- svm(Species ~.,data = iris_train,kernel = "linear",cost = .1,scale = F)
print(svmfit)

plot(svmfit,iris_train[,col])

#1.3调试代价函数，
# cost能实现SVM对分类误差及分离边界的控制，
# 如果cost比较小，分离间隔会比较大（软间隔），产生比较多的被错分样本；
# 相反当加大cost时，会缩小分类间隔，从而减小错分样本。

# kernel <- c("linear","polynomial","radial","sigmoid")
# 核函数的选择：一般用线性核和高斯核，也就是Linear核与RBF核
# 时间上RBF会耗费更多，但是效果不会比线性差
tuned <- tune(svm,Species~.,data = iris_train,
              kernel = "linear",
              ranges = list(cost=c(0.001,0.01,0.1,1,10,100)))
summary(tuned)
# tune:计算出一个最优的cost值

svmfit <- svm(Species ~.,data = iris_train,kernel = "linear",cost = 10,scale = F)
print(svmfit)

#1.4预测
p <- predict(svmfit,iris_test[,col],type="class")
plot(p)

#1.5评价
table(p,iris_test[,3])
mean(p == iris_test[,3])

## 2. 条件推断树(conditional inference tree)
#条件决策树是经典传统决策树的一种重要变体。
#条件推断树与传统决策树类似，但变量和分割的选取是基于显著性检验的，
#而不是纯净度或同质性一类的度量。显著性检验是置换检验。
rm(list=ls())
#install.packages("party")
library(party)
#readingSkills的R内置数据集来创建决策树。 
#它描述了某人的readingSkills的分数，如果我们知道变量“年龄”，“shoesize”，“分数”，以及该人是否为母语者。
print(head(readingSkills))
# Create the input data frame.
input.dat <- readingSkills
# Create the tree.
output.tree <- ctree(
  nativeSpeaker ~ age + shoeSize + score, 
  data = input.dat)
# Plot the tree.
plot(output.tree)
####################################
## 3.随机森林 
####################################
rm(list=ls())
library(ggplot2)
library(cowplot)
#install.packages("randomForest")
library(randomForest)
#3.1 载入数据
url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data"
data <- read.csv(url,header=F)
colnames(data) <- c("age","sex","cp","trestbps","chol","fbs","restecg","thalach","exang","oldpeak","slope","ca","thal","hd")
str(data)

# 3.2清洗数据
data[data == "?"] <- NA
data[data$sex == 0,]$sex <- "F"
data[data$sex == 1,]$sex <- "M"

data$sex <- as.factor(data$sex)
data$cp <- as.factor(data$cp)
data$fbs <- as.factor(data$fbs)
data$restecg <- as.factor(data$restecg)
data$exang <- as.factor(data$exang)
data$slope <- as.factor(data$slope)

data$ca <- as.integer(data$ca)
data$ca <- as.factor(data$ca)

data$thal <- as.integer(data$thal)
data$thal <- as.factor(data$thal)

data$hd <- ifelse(test=data$hd == 0 ,yes="Healthy",no="Unhealthy")
data$hd <- as.factor(data$hd)
str(data)

#3.3建模
set.seed(123)
#set.seed():设定生成随机数的种子，种子是为了让结果具有重复性。
#如果不设定种子，生成的随机数无法重现。
#rfImpute()函数可为存在缺失值的数据集进行插补（随机森林法），得到最优的样本拟合值
data.imputed <- rfImpute(hd~.,data = data,iter=6)

model <- randomForest(hd ~.,data=data.imputed,proximity=T)
model

#3.4 调试
#3.4.1确定500棵树是否足够
# oob表示没有用于构建数的数据
# err.rate:预测错误率
oob.err.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate),times=3),
  Type=rep(c("OOB","Healthy","Unhealthy"),each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
          model$err.rate[,"Healthy"],
          model$err.rate[,"Unhealthy"]))
png("./class5_randomforest_ooberr_ggplot.png")
ggplot(data=oob.err.data,aes(x=Trees,y=Error))+
  geom_line(aes(color=Type))
dev.off()

#3.4.2调试构树时的随机变量
oob.values<-vector(length=10)
for(i in 1:10){
  temp.model <- randomForest(hd~.,data=data.imputed,mtry=i,ntree=1000)
  oob.values[i]<-temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values



