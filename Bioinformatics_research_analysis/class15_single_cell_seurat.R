##-----------class15:sincle cell analysis-----------##
rm(list = ls())
#BiocManager::install('multtest')
#install.packages('Seurat')

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

# 1.数据
setwd("./bioinformation_start/Bioinformatics_research_analysis/")
pbmc.data <- Read10X(data.dir = 'hg19')

# 初始化
# 每个细胞表达的基因数要大于等于200
# 每个基因在不少于三个细胞中表达
pbmc <- CreateSeuratObject(pbmc.data,
                          min.cells = 3,
                          min.features = 200,
                          project= "10X PBMC")

# 2.质控
# 线粒体:低质量/垂死细胞通常表现出广泛的线粒体污染
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc,pattern = "^MT-")
VlnPlot(object = pbmc,features = c('nFeature_RNA','nCount_RNA','percent.mt'))

# FeatureScatter:在一组单个单元格上创建两个要素（通常是要素表达）的散点图。 单元格按其标识类别着色。 
# 这两个特征之间的皮尔逊相关性显示在图的上方,相关系数的绝对值越大，相关性越强,相关系数越接近于1或-1，相关度越强，相关系数越接近于0，相关度越弱。
png('./class15_ncount_mt_featureScatter.png',width = 800,height = 800)
FeatureScatter(object = pbmc,feature1 = 'nCount_RNA',feature2 = 'percent.mt')+
  geom_hline(aes(yintercept = 5),colour = '#BB0000',linetype = 'dashed')+
  ggtitle("QCmt")
dev.off()

png('./class15_ncount_nfeature_featureScatter.png',width = 800,height = 800)
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ 
  geom_hline(aes(yintercept=200),colour="#BB0000", linetype="dashed") +
  geom_hline(aes(yintercept=2500),colour="#BB0000", linetype="dashed") +
  ggtitle("QCng")
dev.off()

# Features数量:过滤具有超过2,500或小于200的独特特征计数的细胞(低质量细胞或空液滴通常只有很少的基因,细胞双峰或多重峰可能表现出异常高的基因计数)，
# 线粒体计数> 5％的细胞
pbmc <- subset(x = pbmc,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt > -Inf & percent.mt < 5)

# 3.标准化：LogNormalize，归一化值存储在pbmc[["RNA"]]@data
pbmc <- NormalizeData(object = pbmc,
                      normalization.method = 'LogNormalize')

# 4.识别高度可变的Features
pbmc <- FindVariableFeatures(object = pbmc,mean.function = ExpMean,dispersion.function = LogVMR,nfeatures = 2000)
head(x = HVFInfo(object = pbmc))
top10 <- head(VariableFeatures(pbmc),10)
# 5.scaledata:线性回归分析,中心化和比例化
pbmc <- ScaleData(object = pbmc,vars.to.regress = c("nCounts_RNA", "percent.mt"))
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10,repel =F)
png('./class15_variablefeatures.png',width = 800,height = 800)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# 5.PCA降维
pbmc <- RunPCA(object = pbmc,npcs = 20,verbose = F)
DimPlot(object = pbmc,reduction = 'pca')
FeaturePlot(object = pbmc,features = 'LYZ')
# 热图显示不同PCA的差异基因
png('./class15_pca_20_dimheatmap.png',width = 800,height = 800)
DimHeatmap(object = pbmc,
           dims = 1:20,
           nfeatures = 30,
           cells = 400,
           reduction = 'pca',
           balanced = T)
dev.off()
png('./class15_pca_2_dimheatmap.png',width = 800,height = 800)
DimHeatmap(object = pbmc,
           dims = 2,
           nfeatures = 30,
           cells = 400,
           reduction = 'pca',
           balanced = T)
dev.off()
# 选择后续使用的PCA数量
png('./class15_pca_ElbowPlot.png',width = 800,height = 800)
ElbowPlot(object = pbmc)
dev.off()
# 由图可以观察PC9-10周围的“肘部”，这表明大多数真实信号都是在前10PC中捕获的。

# 6.聚类细胞
pbmc <- FindNeighbors(object = pbmc,
                      reduction = 'pca',
                      dims = 1:10)
pbmc <- FindClusters(object = pbmc,
                     resolution = 0.7,
                     algorithm = 1)
# 生成降维图：TNSE和UMAP
pbmc <- RunTSNE(object = pbmc,dim.use = 1:10,do.fast = T)

png('./class15_tsne_dimPlot.png',width = 800,height = 800)
DimPlot(object = pbmc,reduction = 'tsne')
dev.off()

pbmc <- RunUMAP(object = pbmc,reduction = 'pca',dims = 1:10)
png('./class15_umap.png',width = 800,height = 800)
DimPlot(object = pbmc,reduction = 'umap')
dev.off()

# 7.寻找分群的Marker基因
cluster1.markers <- FindMarkers(object = pbmc,ident.1 = 3,min.pct = 0.25)
# ident.1：选择的cluster

head(cluster1.markers, n =5)
pbmc.makers <- FindAllMarkers(object = pbmc,min.pct = 0.25,only.pos = T,thresh.use = 0.25)

# 展示不同的cluster下的前N个差己基因
pbmc.makers %>% group_by(cluster) %>% top_n(2, avg_logFC)
# %>% :dplyr包的管道函数，将前一步的结果直接传参给下一步的函数，从而省略了中间的赋值步骤，可以大量减少内存中的对象，节省内存。

VlnPlot(object = pbmc,  c("TCL1A", "CD79A"))
FeaturePlot(object = pbmc,
            features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),
            cols = c("grey", "blue"),
            reduction = "tsne")

# 为不同的cluster标记marker
current.cluster.ids <- c(0,1,2,3,4,5,6,7)
new.cluster.ids <- c("Naive CD4 T","Memory CD4 T", "CD14+ Mono", 
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC")

pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident,from = current.cluster.ids,to = new.cluster.ids)
png("./class15_new_clusterids_dimplot.png",width = 800,height = 800)
DimPlot(object = pbmc,reduction = 'tsne',label = T,pt.size = 0.5)
dev.off()

