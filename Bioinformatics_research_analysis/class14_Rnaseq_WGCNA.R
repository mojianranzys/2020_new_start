##---------------class14:RNA高级分析:WGCNA--------------##
rm(list = ls())
#BiocManager::install("WGCNA")
#BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
library(WGCNA)
#允许R程序最大线程运行
enableWGCNAThreads()
setwd("D:/R/R_learn/2020/bioinformation_start/Bioinformatics_research_analysis/")

# 1.数据的载入和处理
# 导入RNAseq的FPKM表达矩阵；
RNAseq_fpkm <- read.table('class14_wgcna_exp.txt',header = T)
dim(RNAseq_fpkm)

#根据MAD(绝对离差中位数)筛选top5000个基因的表达矩阵，并转置
datExpr = t(RNAseq_fpkm[order(apply(RNAseq_fpkm,1,mad),decreasing = T)[1:5000],])
#apply(X, MARGIN, FUN, ...):MARGIN:c(1,2),1:rows,2:columns
datExpr[1:4,1:4]

# 导入临床表型数据
datTraits <- read.table('class14_match_group.txt',header = T)

# 针对样本做一个系统聚类
datExpr_tree <- hclust(dist(datExpr),method = 'average')
par(mar = c(0,5,2,0))
plot(x = datExpr_tree,
     main = 'sample clustering',
     sub = '',xlab = '',cex.lab = 1.5)

# 添加对应的颜色
sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)),
                                colors = c("thistle3","wheat2","plum","skyblue","azure"),
                                signed = F)
png('./class14_datExpr_tree_plotDendroAndColors.png',width = 800,height = 800)
plotDendroAndColors(datExpr_tree,
                    colors = sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1,4,3,1),
                    cex.rowText = 0.01,
                    main = 'sample dendrogram and traint heatmap')
dev.off()

##-------------如果需要剔除离群样本-----------##
clust = cutreeStatic(datExpr_tree,cutHeight = 300000,minSize = 10)
#cutHeight参数即聚类树被砍掉的高度，minSize设定最小分支数量
table(clust)
# 0：表示需要剔除的样本
keepsample <- (clust == 1)
datExpr2 <- datExpr[keepsample,]
dim(datExpr2)
# 聚类
datExpr_tree2 <- hclust(dist(datExpr2),method = 'average')
plot(x = datExpr_tree2,
     main = 'sample clustering',
     sub = '',xlab = '',cex.lab = 1.5)
##------------------------——-----------------##

# 确定表型性状与样本名称，关联表达矩阵和表型信息
sampleNames <- rownames(datExpr)
traitRows = match(sampleNames,datTraits$gsm)
rownames(datTraits) = datTraits[traitRows, 1] 

# 2.选择最佳的beta值，筛选软阈值,使样本变成无尺度分布，即少部分基因处于绝对优势的位置(表达量高)
powers <- c(c(1:10),seq(from = 12,to = 20,by = 2))

#pickSoftThreshold:针对多种软阈值能力的无标度拓扑分析
sft <- pickSoftThreshold(datExpr,powerVector = powers,verbose = 5)

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = 'soft threshold(power)',
     ylab = 'Mean Connectivity',
     type = 'n',
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels = powers,cex = 0.9,col = 'red')

# 查看推荐的最佳软阈值
sft$powerEstimate

# 3.构建加权共表达网络
# 网络构建

net = blockwiseModules(datExpr = datExpr,
                       power = sft$powerEstimate,
                       maxBlockSize = 6000,
                       minModuleSize = 30,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = T,
                       pamRespectsDendro = F,
                       saveTOMs = F,
                       verbose = 3)

#查看模块聚类情况
table(net$colors)

# 4.模块可视化

mergedColors = labels2colors(net$colors)
table(mergedColors)
png('./class14_module_plotDendroAndColors.png',width = 800,height = 800)
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    'Module colors',
                    dendroLabels = F,
                    hang = 0.03,
                    addGuide = T, 
                    guideHang = 0.05)
dev.off()
#注：灰色默认是无法归类于任何模块的那些基因，如果灰色模块里面的基因太多，
#那么前期对表达矩阵挑选基因的步骤可能就不太合适。

# 5.模块和形状的关系
# 是针对于连续变量（数值型特征），如果是分类变量，需要转换成0-1矩阵的形式方可使用是针对于连续变量（数值型特征），
#（1表示属于此组或有此属性，0表示不属于此组或无此属性）

# 查看分组情况
table(datTraits$subtype)

# 明确基因和样本的数量
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# 把性状转为0-1矩阵
design <- model.matrix(~0 + datTraits$subtype)
colnames(design) <- levels(datTraits$subtype)
moduleColors <- labels2colors(net$colors)

# moduleEigengenes:计算给定单个数据集中模块的模块本征基因（第一主成分）。
Mes0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
head(Mes0)
MEs <- orderMEs(MEs = Mes0)
head(MEs)
moduleTraitCor <- cor(MEs, design , use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix <- paste(signif(moduleTraitCor,2),"\n(",
                    signif(moduleTraitPvalue,1),")",
                    sep = '');
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6,8.5,3,3))
png('./class14_module_trait_relationships.png',width = 1000,height = 800)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module trait relationships")
dev.off()

# 6.感兴趣性状的模块的具体基因分析
# 看step5图发现：Luminal亚型最相关的是brown模块

# 计算模块和基因的相关性矩阵
modNames <- substring(names(MEs),3)
geneModuleMembership <- as.data.frame(cor(datExpr,MEs,use = 'p'))

# 计算每个模块和基因的皮尔森相关系数矩阵(MM)
# MEs是每个模块在每个样本里面的;datExpr是每个基因在每个样本的表达量
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples = nSamples))
names(geneModuleMembership) <- paste("MM",modNames,sep = '')
names(MMPvalue) <- paste('p.MM', modNames,sep='')

# 计算性状与基因的相关性矩阵(GS)
Luminal <- as.data.frame(design[,3])
names(Luminal) <- 'Luminal'

geneTraitSignificance <- as.data.frame(cor(datExpr,Luminal,use = 'p'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) = paste("GS.", names(Luminal), sep="")
names(GSPvalue) = paste("p.GS.", names(Luminal), sep="")

# 把两个相关性矩阵联合起来，指定感兴趣的模块进行分析
module <- 'brown'
column <- match(module,modNames)
moduleGenes = moduleColors==module

par(mfrow = c(1,1));
verboseScatterplot(x = abs(geneModuleMembership[moduleGenes,column]),
                   y = abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Luminal",
                   main =  paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab=1.2,cex.axis=1.2,col = module)

# 7.可视化

# 明确基因和样本数量
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
geneTree <- net$dendrograms[[1]]

# 计算距离矩阵
dissTom <- 1- TOMsimilarityFromExpr(datExpr = datExpr,power = 6)
plotTOM <- dissTom^7
#diag:提取或替换矩阵的对角线，或构造对角线矩阵。
diag(plotTOM) = NA

# 随机选择400个基因画拓扑重叠热图
nSelect = 400
set.seed(10)     # 设计随机种子
select = sample(nGenes,size = nSelect)
selectTOM <- dissTom[select,select]

# 重新聚类
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];

# 打开绘画窗口
sizeGrWindow(9,9)
plotDiss <- selectTOM^(sft$powerEstimate)
diag(plotDiss) <- NA

# TOM plot:行列同时做层次聚类，比较耗时
png('./class14_TOMplot.png',width = 800,height = 800)
TOMplot(dissim = plotDiss,
        dendro = selectTree,
        Colors = selectColors,
        main = "Network heatmap plot, selected genes")
dev.off()

# Recalculate module eigengenes(基因和样本构成矩阵)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 只有连续型性状才能分析
MET = orderMEs(cbind(MEs, Luminal))

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
sizeGrWindow(5,7.5)
par(cex = 0.9)
png('./class14_module_relationships_plotEigengeneNetworks.png',width = 800,height = 800)
plotEigengeneNetworks(multiME = MET,
                      setLabels = '',
                      marDendro = c(0,4,1,2),
                      marHeatmap = c(3,4,1,2),
                      cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()

# 单独绘制模块系统进化树(层次聚类),不带热图
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(multiME = MET,
                      setLabels = "Eigengene dendrogram",
                      marDendro = c(0,4,2,0),
                      plotHeatmaps = F)

# 单独绘制模块之间的相关性热图
par(cex = 1.0)
plotEigengeneNetworks(multiME = MET,
                      setLabels = "Eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = F,
                      xLabelsAngle = 90)

# 8.提取指定模块的基因名
module <- 'brown'
probes <- colnames(datExpr)
inModule <- moduleColors == module
modprobes <- probes[inModule]

# 9.模块的导出
# 主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape,VisANT等网络可视化软件
TOM <- TOMsimilarityFromExpr(datExpr = datExpr,power = 6)
modTOM <- TOM[inModule,inModule]
dimnames(modTOM) <- list(modprobes,modprobes)

## 导出到VisANT的结果提取
vis <- exportNetworkToVisANT(adjMat = modTOM,
                             file = paste("class14_VisANTInput-", module, ".txt", sep=""),
                             weighted = T,
                             threshold = 0)

# 导入到cytoscape的结果提取
cyt <- exportNetworkToCytoscape(adjMat = modTOM,
                                edgeFile = paste("class14_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("class14_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                weighted = T,
                                threshold = 0.02,
                                nodeNames = modprobes,
                                nodeAttr = moduleColors[inModule])













