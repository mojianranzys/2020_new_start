##-------------class16:单细胞拟时分析：monocle3

# 更新R到最新版(用RGui不是RStudio)
if(F){
  install.packages("installr")
  library("installr")
  updateR()
}

#安装相关包
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))

#安装monocle3
#install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install('monocle3')
