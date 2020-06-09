rm(list = ls())
#----------生存分析和cox回归分析---------#
#install.packages("survival")
#install.packages("survminer")
library(survival)
library(survminer)

infile <- './bioinformation_start/Bioinformatics_research_analysis/class11_clinic.csv'
dat <- read.csv(infile,header = T,stringsAsFactors = F,row.names = 1)

# 1.K-M曲线：研究生存时间的分布特点，估计生存率及标准误差，绘制生存曲线
#survfit：创建KM生存曲线
fit <- survfit(Surv(Time,Status) ~ Smoke,data= dat)
#查看fit的结果
res.sum <- surv_summary(fit,data = dat)
res.sum
# survdiff:用于不同组的统计检验
survdiff(Surv(Time,Status) ~ Smoke,data = dat)

# 结果可视化
ggsurvplot(fit, data = dat,
           size = 1,
           conf.int = F, # 有无置信区间
           censor.shape = 124,censor.size = 2, # 删失数据的表现形状
           xlim = c(0,72), #横轴范围
           xlab = "Time (Months)",
           break.time.by = 10, #横坐标
           pval = TRUE,pval.size = 4,  #是否显示pvalue
           palette = c('black','red'),  #生存曲线颜色
           # 添加风险表
           risk.table = T, risk.table.col = "strata",
           risk.table.title="risk table",risk.table.height = 0.16,
           tables.theme = theme_cleantable(),risk.table.fontsize=4,
           ggtheme = theme_bw(),risk.table.y.text.col = T)
png('./bioinformation_start/Bioinformatics_research_analysis/class11_survial_ggsurvplot.png',width = 800,height = 800)
ggsurvplot(fit, data = dat,
           size = 1,legend.title = "Smoke",legend.labs = c("No", "Yes"),
           conf.int = T,censor.shape = 124,censor.size = 2,
           xlim = c(0,72),xlab = "Time (Months)",
           break.time.by = 10,pval = TRUE,pval.size = 4,
           palette = c('black','red'),
           risk.table = T,risk.table.col = "strata",risk.table.title="risk table",
           risk.table.height = 0.16,tables.theme = theme_cleantable(),risk.table.fontsize=4,
           ggtheme = theme_bw(),risk.table.y.text.col = T)
dev.off()

# 2. Cox回归分析
# 单因素回归
fit2 <- coxph(Surv(Time,Status) ~ Smoke,data=dat)

# 多因素回归
fit3 <- coxph(Surv(Time,Status) ~ Gender+Age+Smoke+EGFR+KRAS+TP53,data=dat)
coxph_sum <- summary(fit3)
coxph_sum
# coef:回归系数; exp(coef):风险比(HR); se(coef):标准误差；
# z:coef与se(coef)的比; Pr(>|z|):p值; Lower .95,upper .95:风险比的95％置信区间的上限和下限

# 森林图：
ggforest(fit3, data = dat)
# 修改标签
dat2 <- within(dat, {
  Gender <- factor(Gender, labels = c("female", "male"))
  Smoke <- factor(Smoke, labels = c("no", "yes"))
  EGFR <- factor(EGFR, labels = c("None","Mutation"))
  KRAS <- factor(KRAS, labels = c("None","Mutation"))
  TP53 <- factor(TP53, labels = c("None","Mutation"))
})
fit4 <- coxph(Surv(Time, Status) ~ Gender+Age+Smoke+EGFR+KRAS+TP53, data = dat2)
png('./bioinformation_start/Bioinformatics_research_analysis/class11_cox_ggforest.png',width = 800,height = 800)
ggforest(fit4,data=dat2)
dev.off()
# Hazard ratio(HR)：HR=1:无影响;HR<1:保护因素;HR>1:危险因素;

# Nomogram图
library(rms)
#构建Cox比例风险回归模型
fit5 <- cph(Surv(Time,Status) ~ Gender+Age+Smoke+EGFR+KRAS+TP53,data=dat,x = T,y = T,surv = T)
## 打包数据-根据rms包要求
dd=datadist(dat) 
options(datadist="dd") 

## 绘制列线图 
surv<- Survival(fit5) 
# 建立生存函数 
surv1 <- function(x)surv(1*36,lp=x)
surv2 <- function(x)surv(1*60,lp=x)
nom <- nomogram(fit5,fun=list(surv1,surv2),lp = F,
                funlabel = c('3-Year Survival','5-Year Survival'),
                maxscale = 100,
                fun.at = c('0.9','0.85','0.8','0.7','0.6','0.5','0. 4','0.3','0.2','0.1'))
plot(nom, xfrac=0.2,cex.axis=0.7,cex.var=0.8)


## 构建logisitc回归模型
fit6 <- lrm(Surv(Time,Status) ~ Gender+Age+Smoke+EGFR+KRAS+TP53,data=dat) 

## 绘制logisitc回归的风险预测值的nomogram图
nom <- nomogram(fit6, fun= function(x)1/(1+exp(-x)),
                lp=F, funlabel="Risk")
plot(nom)


