# class 1:假设检验与回归分析
rm(list=ls())
library('ggpubr')

##-----------------假设检验-------------------##
# 1. 卡方检验，样本独立性检验

library(help='vcd') #查看该包的基础内容
library('vcd')
mytable<- table(Arthritis$Treatment,Arthritis$Improved)
chisq.test(mytable) # 卡方检验

# 2. 差异检验: t检验，wilcoxon检验
library(help='datasets')
data("ToothGrowth")

# 如果数据分布符合正态分布使用t-test,不符合用wilcoxon
# 密度图
ggdensity(ToothGrowth$len[ToothGrowth$supp == 'OJ'])

ggdensity(ToothGrowth,x='len',color = 'supp',
          fill = 'supp',palette = c('red','darkblue'),
          rug=T,add='median')

# 画直方图
gghistogram(ToothGrowth,x = 'len',color = 'supp',
            fill = 'supp',palette = c('grey','darkblue'),
            rug = T,add= 'median',bins = 100)

# 画QQ图
ggqqplot(ToothGrowth,x = 'len',color = 'supp')

# 用函数查看数据分布
# 统计检验函数查看数据是否符合正态分布，常用的是shapiro.test
# P-value > 0.05表示符合正态分布，小于则不符合，需使用Wilcoxon检验
shapiro.test(ToothGrowth$len)
shapiro.test(ToothGrowth$len[ToothGrowth$supp == 'OJ'])
# tapply 查看supp不同分组的结果
with(ToothGrowth,tapply(len, supp, shapiro.test))
# 结果都不符合正态分布，用wilcoxon检验,
# 若 p-value <0.05,原假设不成立。意味两者分布不同。警告“无法精確計算带连结的p值“这是因为数据中存在重复的值，一旦去掉重复值，警告就不会出现。
wilcox.test(len~supp,data = ToothGrowth)
# alternative表示被则假设，two.sided(缺省)，双边检验(H1:μ≠H0),less表示单边检验(H1:μ<μ0),greater表示单边检验(H1:μ>μ0)
wilcox.test(len~supp,data = ToothGrowth,exact = F,alternative = 't')
# 配对样本
wilcox.test(len~supp,data = ToothGrowth,exact = F, alternative = 'two.sided',paired = T)
# 单样本
wilcox.test(len~supp,data = ToothGrowth,exact = F, alternative = 't')

# 如果两者符合正态分布，则使用T-test,使用之前可以进行方差补齐检验，var.test只适用于2组样本
var.test(len~supp,data=ToothGrowth)
# 独立样本检验
t.test(len~supp,data=ToothGrowth,alternative = 't', var.equal = T)

# 配对的t-test,要求两组差值符合正态分布
diff <- with(ToothGrowth,len[supp == 'OJ'] - len[supp == 'VC'])
shapiro.test(diff)
t.test(ToothGrowth$len,mu = 15,alternative = 't')

# 作图
#箱线图
ggboxplot(ToothGrowth,x = 'supp',y = 'len',
          color = 'supp',xlab = 'Supplements',
          ylab = 'Tooth length',title = 'Tooth Growth')


ggboxplot(ToothGrowth,x="supp",y = 'len',color = 'supp',
          xlab = 'Supplements',ylab = "Tooth Length",title = "Tooth Growth") +
  stat_compare_means(method = "t.test",method.args = list(alternative = "g", var.equal =F),
                     label.x = 2,label.y = 38)

P <- ggpaired(ToothGrowth,x="supp",y = 'len',color = 'supp',
              xlab = 'Supplements',ylab = "Tooth Length",title = "Tooth Growth") +
  stat_compare_means(method = "t.test",method.args = list(alternative = "g", var.equal =F),
                     label.x = 2,label.y = 38)

# ggpar 把legend 挪到右侧
ggpar(P,legend = 'right'
      ,xlab = "Supplements",ylab = "Tooth Length"
      ,title = "Tooth Growth")








