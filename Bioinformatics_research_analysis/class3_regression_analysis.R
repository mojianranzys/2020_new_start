#------------------class:回归分析-------------------#
rm(list=ls())

library(car)
data(women)

# 1. 线性回归
# 1.1 简单的线性回归
fit <- lm(weight ~ height,data = women)
fit$coefficients
summary(fit)
#结果解读：
# Residuals: Min(最小值),1Q(四分位数1/4),Median(中位数),3Q(四分位数3/4),Max(最大值)
# coefficients: Estimate(估值),Std.Error(标准误差),t value(T值),Pr(>|t|)(p-value)，P<0.05则 拒绝原假设
# Multiple R-squared:拟合优度,	Adjusted R-squared:修正的拟合优度，值越高，拟合度越好
# F-statistic:F检验,判断整体的显著性检验，p-value < 0.05 表示显著

#模型诊断 
shapiro.test(residuals(fit))
par(mfrow=c(2,2))
plot(fit)

#1.2 用二次项等式优化
fit2 <- lm(weight ~ height + I(height^2), data=women)
par(mfrow=c(2,2))
plot(fit2)

#这第二组图表明多项式回归拟合效果比较理想，基本符合了线性假设、残差正态性(除了 观测点13)和同方差性(残差方差不变)。
#观测点15看起来像是强影响点(根据是它有较大的 Cook距离值)，删除它将会影响参数的估计。
#事实上，删除观测点13和15，模型会拟合得会更好
newfit <- lm(weight~ height + I(height^2), data=women[-c(13,15),])
par(mfrow=c(2,2))
plot(newfit)

#1.3 多元线性回归
#以基础包中的state.x77数据集为例，我们想探究一个州的犯罪率和其他因素的关系，
#包括 人口、文盲率、平均收入和结霜天数(温度在冰点以下的平均天数)。

#创建了一个名为states的数据框，包含了我们感兴趣的变量
states <- as.data.frame(state.x77[,c("Murder", "Population",
                                     "Illiteracy", "Income", "Frost")])
#多元回归分析中，第一步最好检查一下变量间的相关性。
#cor()函数提供了二变量之间的相 关系数.
cor(states)
#car包中scatterplotMatrix()函数则会生成散点图矩阵。
library(car)
scatterplotMatrix(states,main="Scatter Plot Matrix")
fit3 <- lm(Murder ~ Population + Illiteracy + Income + Frost, data=states)
summary(fit3)
#结果解释
#文盲率上升1%，谋杀率将会上升4.14%，它的系数在p<0.001的水平下显著不为0。
#相反，Frost的系数没有显著不为0(p=0.954)，表明当控制其他变量不变时，Frost与Murder 不呈线性相关。
#总体来看，所有的预测变量解释了各州谋杀率57%的方差。

#1.4 模型比较
states <- as.data.frame(state.x77[,c("Murder", "Population", "Illiteracy", "Income", "Frost")])
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost,
           data=states)
fit2 <- lm(Murder ~ Population + Illiteracy, data=states)
anova(fit2, fit1)
#预测
predict(fit1, states[1:5,])
predict(fit2, states[1:5,])

#1.5有显著交互项的多元回归
fit4 <- lm(mpg ~ hp + wt + hp:wt, data=mtcars)
fit4$coefficients
#展示
library(effects)
plot(effect("hp:wt", fit4,xlevels=list(wt=c(2.2,3.2,4.2))), multiline=TRUE)


# 2. 逻辑回归
# 用于预测分类型变量
rm(list=ls())
library(rcompanion)
library(rms)
library(ResourceSelection)
library(VGAM)

#结果为2分类是使用glm
#glm() -- Generalized linear model
setwd("D:/R/R_learn/2020/bioinformation_start/Bioinformatics_research_analysis")
#手动导入数据框
#缺失值处理
Framingham <- read.csv('./class3_Framingham.csv',stringsAsFactors = F,header = T)
dim(Framingham)
any(is.na(Framingham))
Framingham <- na.exclude(Framingham)
any(is.na(Framingham))

Framingham$sex <- factor(Framingham$sex)
# attach:是对what添加路径索引，避免重复输入what名称
attach(Framingham)

#2.1 建模
#逻辑回归算法使用glm函数
mod1 <- glm(chdfate~sex,family = binomial())
summary(mod1)

mod2 <- glm(chdfate~.,data=Framingham,family = binomial())
summary(mod2)

mod3 <- glm(chdfate~.,data = Framingham[,-3],family = binomial())
summary(mod3)

#2.2 模型比较
anova(mod2,mod3,test = 'Chisq')

#2.3 作图
#数据打包
ddist <- datadist(Framingham)
options(datadist = 'ddist')
# lrm:构建logisitc回归模型
mod4 <- lrm(chdfate~sex + sbp + scl + age +bmi)
# 绘制logisitc回归的风险预测值的nomogram图(列线图)
nom <- nomogram(mod4,fun=plogis,lp=F)
plot(nom)
# detach:是撤销attach()建立的路径索引，往往二者配套使用。
detach(Framingham)

## ----------------作业----------------------
rm(list=ls())
# 题目：使用以下数据建立 身高预测体重模型 ，使用模型预测身高=170cm时，体重=？
# 数目如下：
height <- c(151, 174, 138, 186, 128, 136, 179, 163, 152, 131)
weight <- c(63, 81, 56, 91, 47, 57, 76, 72, 62, 48)
data1 <- cbind(height,weight)
data1 <- as.data.frame(data1)
end <- c(170,0)
data1 <- rbind(data1,end)
attach(data1)
mod5 <- glm(weight~height)
summary(mod5)
predict(mod5, data1[11,])
# 结果：当身高=170时，体重等于76

