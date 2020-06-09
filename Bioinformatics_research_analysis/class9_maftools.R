## --------------------class9:maftools----------------##
rm(list = ls())
library(grid)
if (!require("BiocManager")){
  install.packages("BiocManager")
}
#BiocManager::install('maftools')
library(maftools)
maf <- './bioinformation_start/Bioinformatics_research_analysis/class9_Total.Filter.maf'
laml <- read.maf(maf)
# 1.plotmafsummary:概览maf文件
plotmafSummary(maf = laml)
plotmafSummary(maf = laml,rmOutlier = T,
               addStat = 'median',dashboard = T,
               titvRaw = F)
# titvRaw: 默认为T,T显示总数，F显示百分比
plotmafSummary(laml, dashboard = F)


#2.oncoplot:瀑布图
oncoplot(maf = laml ,top = 20 )

#3.Oncostrip:draw any number of genes using top or genes arguments.
oncostrip(maf = laml, genes=c("TTN","TP53","PTEN","PHLPP1"))

#4.plotTiTv:箱线图
#one.function classifies SNPs into Transitions and Transversions 
#two.returns a list of summarized tables in various ways.
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

#5.lollipopPlot：showing mutation spots on protein structure,eg:TP53
lollipopPlot(maf = laml, gene = 'TP53', AACol = 'AAchange', 
             showMutationRate = TRUE)

#6.rainfallPlot:characterized by genomic loci with localized hyper-mutations
#if detectChangePoints = TRUE：rainfall plot also highlights regions where potential changes in inter-event distances are located.
rainfallPlot(maf = laml, pointSize = 0.6,detectChangePoints = TRUE)

#7.tcgaCompare:draws distribution of variants compiled from over 10,000 WXS samples across 33 TCGA landmark cohorts.
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML')

#8.plotVaf:plots Variant Allele Frequencies as a boxplot
pdf('./bioinformation_start/Bioinformatics_research_analysis/class9_maftools_plotvaf.pdf')
plotVaf(maf = laml, vafCol = 't_AF')
dev.off()
#9.geneCloud:plot word cloud plot for mutated genes
geneCloud(input = laml, minMut = 3)

somaticInteractions(laml,top= 10)

output <- somaticInteractions(laml)
write.csv(output$pair, file="./bioinformation_start/Bioinformatics_research_analysis/class9_somaticInteractions.pairwise.csv", 
            quote=FALSE, row.names=FALSE)

#加入临床数据
clinic <- './bioinformation_start/Bioinformatics_research_analysis/class9_clinical.xls'
laml_clinic <- read.table(clinic,header = T)

oncoplot(maf=laml,annotationDat=laml_clinic,clinicalFeatures='Class')

# 在读入maf文件的同时读入临床信息文件
lamlclin = read.maf(maf= maf, clinicalData = clinic) 
oncoplot(maf=lamlclin,clinicalFeatures=c('Class','Gender'))
# 修改注释的颜色
oncoplot(maf=lamlclin,clinicalFeatures=c('Gender','Class'),
         annotationColor=list(Class=c('C1'='#66C2A5','C2'='#FC8D62','C3'='#8DA0CB'),
                              Gender=c('Male'='#FBB4AE','Female'='#B3CDE3')))
#按照注释临床信息排序（默认按照第一个临床信息）
oncoplot(maf=lamlclin,clinicalFeatures=c('Gender','Class'),
         annotationColor=list(Class=c('C1'='#66C2A5','C2'='#FC8D62','C3'='#8DA0CB'),
                              Gender=c('Male'='#FBB4AE','Female'='#B3CDE3')),
         sortByAnnotation = T)

# 下载安装TCGAmutations的包  ERROR 下载不下来

#----------------- maftools Signature analysis——————————————————————————##
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
# 构建突变矩阵
laml.tnm = trinucleotideMatrix(maf = laml, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

#extractSignatures-用于从96种替换类型中提取突变特征。
#这里nTry=8设定尝试n的最大值，程序将调用NMF计算共生性相关系数并选取最合适的值
#BiocManager::install('NMF')
library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 8)
#结果显示：预测的最佳值是y轴上的相关值显著下降的值。在这种情况下，它在n = 3处

# 绘制突变特征图
laml.sig = extractSignatures(mat = laml.tnm, n = 3)
pdf('./bioinformation_start/Bioinformatics_research_analysis/class9_plotsignatures.pdf')
maftools::plotSignatures(nmfRes = laml.sig, title_size = 0.8)
dev.off()
#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")


