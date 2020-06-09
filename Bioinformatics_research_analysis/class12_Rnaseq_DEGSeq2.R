##---------------class12:RNA基因差异表达分析----------##
rm(list=ls())
library(DESeq2)

dat <- read.table(file = './bioinformation_start/Bioinformatics_research_analysis/class12_Readscount.txt',header = T,row.names = 1)
# 删除基因名称列
database <- round(as.matrix(dat))
# 设置分组信息
condition <- factor(c(rep("C",3),rep('T',3)))
coldat <- data.frame(row.names = colnames(database),condition)

# 构建dds(DESeqDataSet Object)对象
dds <- DESeqDataSetFromMatrix(countData = database,colData = coldat,design = ~condition)
dim(dds)

# 过滤低质量数据:将所有样本基因表达量之和小于1的基因过滤掉
dds_filter <- dds[rowSums(counts(dds)) > 1,]

# DESeq2:差异表达分析
dds2 <- DESeq(dds)

# 提取差异表达分析结果
res1 <- results(dds2)  #默认
res2 <- results(dds2,pAdjustMethod = 'fdr', alpha = 0.05)  #根据padj < 0.01定义显著性
res <- res2[order(res2$padj),]

summary(res)
table(res$padj < 0.05)
write.csv(res,file = './bioinformation_start/Bioinformatics_research_analysis/class12_all_result.csv')
#合并
resdat <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=F)),by = 'row.names',sort=F)

# 标准化
dds_norm <- estimateSizeFactors(dds)
resdat_norm <- merge(as.data.frame(res), as.data.frame(counts(dds_norm, normalized=T)),by="row.names",sort=FALSE)

# 设定阈值，筛选差异基因
diff_gene_deseq2 <- subset(res,padj<0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file= "./bioinformation_start/Bioinformatics_research_analysis/class12_DEG_TVSC.csv")  
# 将筛选的差异基因和表达量结果合并
res_diff_dat <- merge(as.data.frame(diff_gene_deseq2),as.data.frame(counts(dds2,normalize=F)),by='row.names',sort=F)
write.table(res_diff_dat,file= "./bioinformation_start/Bioinformatics_research_analysis/class12_DEG_recount_result.xls",sep="\t", row.names=FALSE,col.names=TRUE,quote = F)

# 1.绘制MA图
library(ggplot2)
png('./bioinformation_start/Bioinformatics_research_analysis/class12_Rnaseq_MA.png',width = 800,height = 800)
plotMA(res,ylim=c(-5,5))
dev.off()
# 2.gene聚类热图
library(genefilter)
library(pheatmap)
rownames(res_diff_dat)=res_diff_dat[,1]
pheatmap(res_diff_dat[8:13],show_rownames=T,scale="row",cex=0.6)
png('./bioinformation_start/Bioinformatics_research_analysis/class12_Rnaseq_heatmap.png',width = 800,height = 800)


pheatmap(res_diff_dat[8:13],show_rownames=T,scale="row",
         cex=0.7,cellheight = 8,cellwidth = 80)
dev.off()
# 3.火山图
df <- data.frame(resdat$padj,resdat$log2FoldChange)
colnames(df) <- c('FDR',"FC")

# 下调基因
df.G <- subset(df,FC < -1 & FDR < 0.05)
df.G <- cbind(df.G,rep("Down",nrow(df.G)))
colnames(df.G)[3] <- "Style"

# 无差异基因
df.B <- subset(df,abs(FC) <=1 | FDR >= 0.05)
df.B <- cbind(df.B,rep("Normal",nrow(df.B)))
colnames(df.B)[3] <- "Style"

# 上调基因
df.R <- subset(df,FC > 1  & FDR < 0.05)
df.R <- cbind(df.R,rep("Up",nrow(df.R)))
colnames(df.R)[3] <- "Style"

# 合并
df.t <- rbind(df.G,df.B,df.R)

# 火山图
png('./bioinformation_start/Bioinformatics_research_analysis/class12_Rnaseq_Volcano.png',width = 800,height = 800)
ggplot(data=df.t, aes(x = FC, y = -log10(FDR), color = Style)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'Normal')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', linetype="longdash",size = 0.3) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', linetype="longdash",size = 0.3) +
  xlim(-5, 5) + ylim(0, 10) +
  labs(x = 'Log2 FC', y = '-Log10 (FDR)', color = '', title = 'CaseVSControl')
dev.off()
