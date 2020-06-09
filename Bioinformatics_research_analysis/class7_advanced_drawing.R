##-------------class7:高级绘图-----------##
rm(list = ls())
library(ggplot2)
data(mpg)
head(mpg) # 耗油量数据集
#hwy：高速公路行驶记录每加仑行驶的英里数（miles per gallon，mpg）
#displ：发动机排量（L）
#drv：动力传动系统（前轮f，后轮r，四轮4）
#class：描述汽车种类的变量（双座，SUV，紧凑型等）

# 散点图
p <- ggplot(data=mpg, aes(x=displ, y=hwy)) + geom_point() 
# 创建p图层，数据集mpg，aes进行数据映射，+ geom_point() 添加散点图图层
print(p) # 绘制图层
p + geom_point(colour="darkblue")  # colour控制颜色属性
p + geom_point(aes(colour=class))     # 不同class进行颜色映射
# 标尺scale 颜色、形状
p + geom_point(aes(colour=class)) + scale_color_brewer(palette = "Set1")  # 修改颜色标尺

# 整合
ggplot(data=mpg, aes(x=displ, y=hwy)) + 
  geom_point(aes(colour=class, shape=class)) + 
  scale_color_brewer(palette = "Set1") + # 注意添加图层“+”号的置于末尾
  scale_shape_manual(values=rep(1:7))  # 修改形状

# 标尺scale 坐标刻度
print(p)
p + scale_x_log10()  # x轴取对数
p + scale_x_continuous(trans="log10") # 同上
p + scale_x_reverse() # X轴颠倒
p + scale_x_sqrt()  # x轴开根号
p + coord_trans(x="sqrt") # 同上

# 主题
ggplot(data=mpg, aes(x=displ, y=hwy)) + geom_point() +
  theme_bw() + # 白色背景
  theme(panel.grid.minor.y =element_blank(), # 删掉网格线
        axis.text.x = element_text(size = 15, color = "darkblue",
                                   face = "bold", vjust = 0.5,
                                   hjust = 0.5, angle = 45))
# face ("plain", "italic", "bold", "bold.italic"):“普通”，“斜体”，“粗体”，“ bold.italic”
#vjust:水平对齐[0，1]；hjust:垂直对齐[0,1];angle:角度[0,360]
ggplot(data=mpg, aes(x=displ, y=hwy)) + geom_point() + theme_classic()
#theme_classic():具有x和y轴线且无网格线的经典外观主题
ggplot(data=mpg, aes(x=displ, y=hwy)) + geom_point() + theme_minimal()
#theme_minimal():没有背景注释的简约主题。

# 图例
ggplot(data=mpg, aes(x=displ, y=hwy, colour=class)) + 
  geom_point() + scale_color_brewer(palette = "Set1") +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        legend.key.size=unit(1,'cm'),
        legend.background=element_rect(colour="black",fill="white"))

# 标签
ggplot(mpg, aes(x=displ, y=hwy)) + xlim(2,6) + geom_point() + ylab("") + 
  geom_hline(yintercept=20) +
  labs(x="displ(L)", title="mpg") + 
  geom_text(x=3, y=40, label="Text 2")
# geom_abline添加斜线， geom_hline添加水平线，geom_vline添加垂直线

# 分面facet
p + geom_point(aes(colour=class)) + facet_wrap(~class)

# 点线图
df <- data.frame(supp=rep(c("VC", "OJ"), each=3),
                 dose=rep(c("D0.5", "D1", "D2"),2),
                 len=c(6.8, 15, 33, 4.2, 10, 29.5))
head(df)
ggplot(df, aes(x=dose, y=len, group=supp)) +
  geom_line(aes(linetype=supp, color = supp))+
  geom_point(aes(shape=supp, color = supp))

# 柱状图/条形图
p <- ggplot(mpg,aes(x=manufacturer,fill=manufacturer)) + geom_bar()
p 
p + theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45)) # 修改坐标轴刻度

ggplot(mpg,aes(x=manufacturer, fill=class)) + 
  geom_bar() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45))

ggplot(mpg,aes(x=manufacturer, fill=class)) + 
  geom_bar(position=position_dodge()) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45))

# 饼图
ggplot(mpg,aes(x=factor(1),fill=manufacturer)) + geom_bar() + 
  coord_polar(theta="y")  # 极坐标转换

# 直方图
ggplot(mpg) + geom_histogram(aes(hwy))
ggplot(mpg) + geom_histogram(aes(x=hwy,fill=class)) + 
  scale_fill_brewer() #根据变量填充颜色
ggplot(mpg) + geom_histogram(aes(x=hwy,fill=class),binwidth=2) + 
  scale_fill_brewer() #binwidth设置组距

ggplot(mpg) + geom_histogram(aes(x=hwy,fill=class),position="fill") +
  scale_fill_brewer() # 按百分比填充
#position="fill":堆叠图形元素并将高度标准为1

ggplot(mpg) + geom_histogram(aes(x=hwy,fill=class),position="dodge") + 
  scale_fill_brewer() + xlim(25,28)
# position="dodge":避免重叠，并排放置

# 密度图
ggplot(mpg) + geom_density(aes(x=hwy,colour=class,fill=class),alpha=0.5)

# 箱线图
ggplot(mpg) + geom_boxplot(aes(x=class,y=hwy,fill=class))

# 提琴图
ggplot(mpg) + geom_violin(aes(x=class,y=hwy,fill=class))

ggplot(mpg,aes(x=class,y=hwy,fill=class)) + 
  geom_violin()+
  stat_summary(fun = mean,geom="point") # stat 统计均值

# homework
library(ggplot2)
ggplot(mpg,aes(x=class,y=hwy,fill=class)) + 
  geom_violin() +
  stat_summary(fun = mean,geom="point") +
  stat_summary(fun = mean,geom="line",group = 1) +
  labs(title="mpg") +
  theme(plot.title = element_text(hjust = 0.5))


