## --------------------class10:突变特征及分子分型----------------##
rm(list = ls())
#1. single signature
#install.packages('deconstructSigs')
library("deconstructSigs")
library("BSgenome.Hsapiens.UCSC.hg19")

dat <- read.csv('./bioinformation_start/Bioinformatics_research_analysis/class10_test.input.csv',stringsAsFactors = F,head=T)

#mut.to.sigs.input ：输出96种点突变的统计文件
sigs.input <- mut.to.sigs.input(mut.ref = dat,sample.id="Sample",
                                chr="chr",pos="pos",ref="ref",alt="alt",
                                bsg=BSgenome.Hsapiens.UCSC.hg19)

#whichSignatures : 推断signature的组成
test = whichSignatures(tumor.ref = sigs.input,
                       signatures.ref = signatures.cosmic,
                       sample.id = 1,contexts.needed = TRUE,
                       tri.counts.method = 'exome2genome',
                       signature.cutoff=0)
#contexts.needed：输入数据只包含突变上下三碱基的突变数目，则必须进行数据标准化，选TURE,
#tri.counts.method：WES: exome2genome

write.csv(test$weights,file='./bioinformation_start/Bioinformatics_research_analysis/class10_test.signature.weights.csv',quote=F)

#pdf("./bioinformation_start/Bioinformatics_research_analysis/class10_plotsignatures.pdf")
plotSignatures(test,sub="test")
#dev.off()

makePie(test,sub="test")

#2. hclust
library(dendextend)
df <-  data.frame(num=c(16.9,38.5,39.5,80.8,82,34.6,116.1))
row.names(df) <- c('A','B','C','D','E','F','G')
## 计算欧式距离
out.dist <- dist(df, method = "euclidean")
## 根据距离聚类
out.hclust=hclust(out.dist,method="average")
plot(out.hclust)
## 分成两组
ct = cutree(out.hclust,2)
## 可视化分组
rect.hclust(out.hclust, k=2)
## dendextend R 包 它提供了剪切树状图的功能，以及着色标签和分支
dct = as.dendrogram(out.hclust)
dct = color_branches(dct, k = 2)
plot(dct)

# 3. ConsensusClusterPlus
#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

dat2 <- read.table(file = './bioinformation_start/Bioinformatics_research_analysis/class10_all_signature.weights.xls',sep = '\t',header = T,row.names = 1)
dat2 <- as.matrix(dat2)
setwd('D:/R/R_learn/2020/bioinformation_start/Bioinformatics_research_analysis/')
Signature.res = ConsensusClusterPlus(dat2,maxK=9,reps=1000,pItem=0.8,
                                     pFeature=1,title="SignatureConsensus",
                                     clusterAlg='hc',distance='euclidean',
                                     seed=666666,plot="pdf",innerLinkage="ward.D",
                                     finalLinkage="ward.D",writeTable=TRUE)
setwd('D:/R/R_learn/2020/')
# maxK:分型采用的最大聚类数目
# reps:抽样重复次数，建议1000
# pltem:抽样包含的样本比例，建议0.8
# pFeature:抽样包含的特征比例，建议1
# clusterAlg:聚类算法：c('hc', 'pam', 'km', 'kmdist')
# distance:聚类采用的距离矩阵, c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski" or 'custom')
# innerLinkage :抽样样本聚类的联接方法
# finaLinkage:一致性矩阵聚类的联接方法
# seed:设置随机种子，以获得可重复的结果

# 文献参考参数：
# RNA -> hc, pearson(皮尔森), ward.D innerlinkage algorithm;
# CNV -> hc, euclidean(欧氏距离), ward.D innerlinkage algorithm; 
# Signature -> hc, euclidean, ward.D innerlinkage algorithm

# ComplexHeatmap:对结果可视化
library(ComplexHeatmap)
library(dendextend)
library(circlize)
cluster_color = c("#FF3232","#00D9F9","#AAD870","#F774FF","#0055FF","#FFA358","#FF61B4","#C15555","#FFEB00","#43B6D5")
cluster_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Cmatrix = read.csv('./bioinformation_start/Bioinformatics_research_analysis/SignatureConsensus/SignatureConsensus.k=6.consensusMatrix.csv',header = T,row.names = 1)
group = read.csv('./bioinformation_start/Bioinformatics_research_analysis/SignatureConsensus/SignatureConsensus.k=6.consensusClass.csv',header = F)
Rawdata=read.table("./bioinformation_start/Bioinformatics_research_analysis/class10_all_signature.weights.xls", header=T, row.names=1,sep="\t")
## 权重数据中心化及标准化
## 中心化：将数据减去均值
## 标准化：在中心化的基础上,再除以数据的标准差
Zdata<-scale(Rawdata)
dat3 <-as.matrix(Zdata)

## 对Cmatirx行名列名重命名，行名列名均为样本名
colnames(Cmatrix)=colnames(Rawdata)
rownames(Cmatrix)=colnames(Rawdata)
## 重命名group列名
colnames(group) = c("Sample","Cluster")
## 将group格式转为data.frame 并把Sample列作为行名
df = data.frame(Subtype=group$Cluster)
rownames(df)=group$Sample

## hclust构建聚类树
## 查看ConsensusClusterPlus聚类 
## method 和 finalLinkage="complete"要一致
column_hc = hclust(as.dist(1 - Cmatrix), 
                   method = "ward.D")
## cutree 分成六组
column_ct = as.matrix(cutree(column_hc,6))
# 提取聚类树样本顺序
orderS=column_hc$labels[column_hc$order]
# 对cutree分组按照聚类树样本顺序排序
orderM = as.matrix(column_ct[orderS,])
colnames(orderM)="orderG"
orderG = as.vector(orderM[,1])
## 聚类树着色
column_dend = as.dendrogram(column_hc)
column_dend = color_branches(column_dend, k = 6,
                             col = cluster_color[1:6])

ug=unique(orderG)
uc=cluster_color[1:6]
#设置column bar分组颜色 和聚类数着色需要一致
anno_color=colorRamp2(ug,uc)

column_Anno = HeatmapAnnotation(df = df,
                                col = list(Subtype =anno_color),
                                show_legend =FALSE, annotation_height = 0.2,
                                show_annotation_name = TRUE)
#设置图例，从图左到右顺序命名
col_fun = cluster_fun
lgd1 = Legend(at = c(-2,0,2),col_fun = col_fun,title='Signature Weight')
lgd2 = Legend(labels = c("C1","C2","C3","C4","C5","C6"),
              title ="Mutation Subtype", 
              legend_gp= gpar(fill = cluster_color[1:6]))
lgd_list = list(lgd1,lgd2)

ht <- Heatmap(dat3, name = "heatmap", col = col_fun,
              cluster_columns = column_dend,
              column_dend_reorder=FALSE,
              column_dend_height = unit(3, "cm"),
              row_dend_width = unit(2, "cm"),
              top_annotation =column_Anno,
              show_heatmap_legend = FALSE,
              show_column_names = FALSE,
              show_row_names = T)
png('./bioinformation_start/Bioinformatics_research_analysis/class10_mutation_Subtype.heatmap.png',
    width=4600,height=3600,type="cairo-png",res=72*5)
draw(ht, annotation_legend_list = lgd_list)
dev.off()
