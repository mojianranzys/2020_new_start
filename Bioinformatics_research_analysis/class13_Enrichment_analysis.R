##-------------基因功能富集分析-----------##
rm(list = ls())
# 1.安装富集分析需要的包
# 富集分析
#BiocManager::install('clusterProfiler')
library(clusterProfiler)

# 人类基因注释包
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)

# GO富集分析
#BiocManager::install('topGO')
library(topGO)

# 反应组学
#BiocManager::install('ReactomePA')
require(ReactomePA)

# 通路
#BiocManager::install('pathview')
require(pathview)

require(DOSE)
library(DO.db)
library(enrichplot)
library(ggplot2)
# 2.分析
# 查看数据库的ID类型
keytypes(org.Hs.eg.db)
gene <- read.table('./bioinformation_start/Bioinformatics_research_analysis/class13_input.txt',header = T,row.names = 1)
symbol <- rownames(gene)

# 将Symbol 转换为EntrezID
id <- bitr(symbol,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Hs.eg.db',drop = T)

# 3.enrichGO：GO富集分析
go <- enrichGO(gene = id$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = 'ALL',
               pvalueCutoff = 0.05,
               pAdjustMethod = 'BH')
#ont:c("BP", "MF", "CC"), or "ALL" for all three

write.table(go,"./bioinformation_start/Bioinformatics_research_analysis/class13_go_Enrich.xls",sep="\t",row.names=FALSE,col.names=TRUE,quote=F)
table(go$ONTOLOGY)

# 可视化：有向无环图，从上至下所定义的功能范围变小，选取 GO 富集分析的前10个结果作为图的主节点。
# 从黄色过滤到红色，对应p值从大到小，颜色越深代表富集程度越高

# 柱状图
barplot(go,showCategory = 10)

# 散点图
dotplot(go,showCategory = 10)
dotplot(go,showCategory=10,title="GO enrichment Top10")

# 分别对BP.CC.MF做富集分析
egoBP<-enrichGO(gene=id$ENTREZID, 
                OrgDb = org.Hs.eg.db, 
                ont='BP',
                pAdjustMethod = 'BH',
                pvalueCutoff = 0.05)
egoCC<-enrichGO(gene=id$ENTREZID, 
                OrgDb = org.Hs.eg.db, 
                ont='CC',
                pAdjustMethod = 'BH',
                pvalueCutoff = 0.05)
egoMF<-enrichGO(gene=id$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont='MF',
                pAdjustMethod = 'BH',
                pvalueCutoff = 0.05)
# 绘制DAG图
p1 <- plotGOgraph(egoCC)
p2 <- plotGOgraph(egoMF)
p3 <- plotGOgraph(egoBP)

# 绘制富集图
# 提取CC中top10 gene
t1 <- egoCC[order(egoCC[,5]),]
picCC <- t1[1:10,c(2,3)]
picCC$pic <- parse_ratio(picCC[,2])
picCC$type<-"CC"

# MF
t2 <- egoMF[order(egoMF[,5]),]
picMF <- t2[1:10,c(2,3)]
picMF$pic <- parse_ratio(picMF[,2])
picMF$type<-"MF"

# BP
t3 <- egoBP[order(egoBP[,5]),]
picBP <- t3[1:10,c(2,3)]
picBP$pic <- parse_ratio(picBP[,2])
picBP$type<-"BP"

# 合并CC,MF,BP的Top10的富集结果
picture_dat <- rbind(picCC,picMF,picBP)
picture_dat$Description <- as.factor(picture_dat$Description)
picture_dat$type <- as.factor(picture_dat$type)

library(RColorBrewer)
r <- brewer.pal(3,'Set2')

# ggplot 可视化
g1 <- ggplot(picture_dat,aes(x = picture_dat$Description,y = picture_dat$pic))+
  geom_bar(stat="identity",aes(fill=picture_dat$type))+
  labs(fill="Type")+
  xlab("GO term")+
  ylab("GeneRatio")+
  scale_fill_manual(values = r)+
  theme(panel.grid.major=element_line(colour=NA))+
  facet_wrap(~picture_dat$type,scales="free")+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(panel.border=element_blank())+
  theme(strip.background = element_rect(colour="white", fill="white"))
print(g1)

# 4.KEGG富集分析
kegg <- enrichKEGG(gene = id$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH')

# 可视化
barplot(kegg,showCategory = 10)
dotplot(kegg,showCategory = 10)

require(pathview)
a <- as.data.frame(cbind(row.names(gene),gene$log2FC))
colnames(a) = c('SYMBOL',"LOGFOLD")

# 根据SYMBOL合并转换后的ID合并
r2 = merge(x=id,y=a,by='SYMBOL',all.x=T)
genelist <- as.numeric(r2$LOGFOLD)
names(genelist) <- r2$ENTREZID

#可视化
setwd('./bioinformation_start/Bioinformatics_research_analysis/')
hsa04062 <- pathview(gene.data = genelist,
                      pathway.id = "hsa04062", 
                      species="hsa",
                      plot.col.key=FALSE, 
                      limit=list(genelist = max(abs(genelist))))

setwd('../../../2020/')
# reactome富集分析
require(ReactomePA)
reac <- enrichPathway(gene = id$ENTREZID,
                      organism = 'human',
                      pvalueCutoff = 0.05,
                      pAdjustMethod = 'BH',
                      readable = T)
head(reac)
barplot(reac,showCategory = 10)
dotplot(reac,showCategory = 10)

# 5.对差异基因进行GSEA分析
library(enrichplot)
if(F){
  gene2 <- read.table('./bioinformation_start/Bioinformatics_research_analysis/class13_input2.txt',header = T)
  #提取log2FC的值，并排序
  gene2 <- gene2[order(gene2$log2FC,decreasing = T),] 
  row.names(genelist) <- genelist[,1]
  genelist <- as.data.frame(gene2$log2FC)
  row.names(genelist) <- gene2$GeneName
  genelist <- as.numeric(genelist)
}
gene2 <- read.table('./bioinformation_start/Bioinformatics_research_analysis/class13_input2.txt',header = T)
##提取log2FC的值
gene2<- gene2[order(gene2$log2FC,decreasing = T),] 
##匹配基因名称
genelist <- gene2$log2FC
names(genelist) <- gene2$GeneName
##根据log2FC进行排序
genelist

# gesgo做富集分析
gsemf <- gseGO(genelist, 
               OrgDb = org.Hs.eg.db, 
               keyType = "SYMBOL", ont="MF", 
               pAdjustMethod = "BH" , 
               nPerm = 1000, 
               minGSSize = 10, 
               maxGSSize = 1000, 
               pvalueCutoff = 0.05)
write.table(gsemf,"./bioinformation_start/Bioinformatics_research_analysis/class13_go_gsea.xls",sep="\t",row.names=FALSE,col.names=TRUE,quote=F)

# 可视化

# gseaplot
gseaplot(gsemf,geneSetID = 'GO:0003674')

# gseaplot2 -------ERROR
gseaplot2(gsemf, geneSetID="GO:0003674",pvalue_table = T)
gseaplot2(gsemf, geneSetID=c("GO:0003674","GO:0005515","GO:0005488"))


