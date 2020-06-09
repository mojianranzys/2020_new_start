## --------- class2:方差分析 -----------
rm(list = ls())
library(ggpubr)
library(ggsci)
library(ggsignif)
library(car)
library(userfriendlyscience)
data("ToothGrowth")

# 将计量转为因子
ToothGrowth$dosef <- factor(ToothGrowth$dose,ordered = T)
# 1.1:正态性检验
with(ToothGrowth,tapply(len, dosef, shapiro.test))

# 1.2 方差齐性检验，三组及以上可以用levene检验,value > 0.5 可认为等方差。
leveneTest(len~dosef,ToothGrowth)

# 1.3 单因素方差分析
# anova:Summary()函数会输出残差和模型
AOV1 <- aov(len~dosef,ToothGrowth)
summary(AOV1)

# 1.4 诊断模型
# 查看残差正态性

res1 <- residuals(AOV1)
res1 <- AOV1$residuals
shapiro.test(res1)
ggqqplot(res1)

#各个计量组之前有差异，是可信的，
#但是知道具体是那些组之间的差异。
#需要使用额外的检验。
#1.5 两两比较
Tu <- TukeyHSD(AOV1)
P <- Tu$dosef[,4]
format(P,scientific = F)

# 1.6 作图
# 箱线图
ggboxplot(ToothGrowth,x = 'dosef',y = 'len',fill = 'dosef',palette = 'npg') + stat_compare_means(method='anova')

comp <- list(c('0.5','1'),c('0.5','2'),c('1','2'))

ggplot(ToothGrowth,aes(x=dosef,y=len,fill=dosef))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot()+
  geom_point(size=2,alpha = 0.6)+ #alpha 设置透明度,范围是0到1，全透明到不透明
  geom_signif(comparisons = comp,y_position = c(36,40,38),
              annotations = format(unname(P),scientific = T,digits = 3))+
  theme_classic2(base_size = 16) # 创建具有轴线的经典主题

# 1.7 符合正态分布，不符合方差齐性
# 使用 welch ANOVA
oneway.test(len~dosef,ToothGrowth,var.equal = F)
# 两两比较
GH <- with(ToothGrowth,posthocTGH(len,dosef,method='games-howell'),digits=3)
GH$output

# 2.不符合正态分布使用Kruskal wallis检验
kruskal.test(len~dosef,data = ToothGrowth)
with(ToothGrowth,pairwise.wilcox.test(len,dosef,exac = F))
# 3.数据集中计量和给药途径都应该考虑
# 3.1 双因素方差分析
AOV2 <- aov(len~dosef*supp,ToothGrowth,)
AOV2<- aov(len~dosef + supp:dosef,ToothGrowth) # 同上
summary(AOV2)

# 3.2 诊断模型
# 查看残差正态性
res2 <- AOV2$residuals
shapiro.test(res2)
ggqqplot(res2)
leveneTest(res2~dosef*supp,ToothGrowth)

#3.3两两比较
TK <- TukeyHSD(AOV2,'dosef:supp')
TKint <- as.data.frame(TK$`dosef:supp`)
TKint$pair <- rownames(TKint)
TKint$sig <-cut(TKint$`p adj`,c(0,0.05,1),labels =c('p<0.05','NS'),right=F)

#3.4 作图
ggplot(ToothGrowth,aes(x=dosef,y=len,fill=supp))+
  stat_boxplot(geom='errorbar')+
  geom_boxplot()

# 画两比较
plot(TK)
mytheme <- theme_classic2(base_size = 16)
ggplot(TKint,aes(x=`p adj`,y=pair,color=sig))+
  geom_errorbarh(aes(xmin = lwr,xmax = upr),height=0.3)+
  geom_vline(xintercept=0,lty=2)+
  mytheme

