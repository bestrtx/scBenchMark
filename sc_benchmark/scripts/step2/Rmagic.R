source("R_conf.R")
Sys.setenv(RETICULATE_PYTHON = "/mnt/data1/LJR/software/miniconda3/envs/bestgtx/bin/python3")
system("bash ./sc_benchmark/scripts/kill_user_usage.sh")

# BiocManager::install("Rmagic")
library(Rmagic)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
load(datapath)
dat=t(dat)

bmmsc <- dat
bmmsc[1:5,1:10]

system("bash ./sc_benchmark/scripts/start_monitor.sh -n ./sc_benchmark/output/step2/Rmagic/usage -t 5")
start_time <- proc.time()

# keep genes expressed in at least 10 cells
keep_cols <- colSums(bmmsc > 0) > 10
bmmsc <- bmmsc[,keep_cols]
# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(bmmsc)), bins=50) +
  geom_vline(xintercept = 1000, color='red')
# keep cells with at least 1000 UMIs
keep_rows <- rowSums(bmmsc) > 1000
bmmsc <- bmmsc[keep_rows,]
# Normalizing data
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)
# Using MAGIC for downstream analysis
bmmsc_MAGIC <- magic(bmmsc, genes="all_genes",t=4)

end_time <- proc.time()
execution_time <- end_time - start_time
system("bash ./sc_benchmark/scripts/kill_user_usage.sh")

save(bmmsc_MAGIC, file = "./sc_benchmark/output/step2/Rmagic/MAGIC_data.RData")
load("./sc_benchmark/output/step2/Rmagic/MAGIC_data.RData")
as.data.frame(bmmsc_MAGIC)[1:5, 1:10]
dat = as.data.frame(bmmsc_MAGIC)
dat=t(dat)


load(lablepath)

scRNA <- CreateSeuratObject(counts = dat, project = "Grubman", min.cells = 1, min.features = 1); scRNA

#使用从MT-作为线粒体基因集
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "MT_")
#qc指标可视化
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & nCount_RNA > 500)

scRNA <- ScaleData(scRNA)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- RunPCA(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:10)
scRNA <- RunUMAP(scRNA, dims = 1:10)

for (res in c(0.05, 0.085, 0.1, 0.15, 0.25, 0.5)) {
  scRNA <- FindClusters(scRNA, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.

plot_grid(ncol = 3, DimPlot(scRNA, reduction = "umap", group.by = "RNA_snn_res.0.15") + ggtitle("louvain_0.15"), 
          DimPlot(scRNA, reduction = "umap", group.by = "RNA_snn_res.0.25") + ggtitle("louvain_0.25"), 
          DimPlot(scRNA, reduction = "umap", group.by = "RNA_snn_res.0.5") + ggtitle("louvain_0.5"))

png("./sc_benchmark/output/step2/Rmagic/MAGIC_umap.png", height=768, width=768)
DimPlot(scRNA, reduction = "umap", group.by = "RNA_snn_res.0.25") + ggtitle("louvain_0.25")
dev.off()

c = paste(scRNA@assays[["RNA"]]@counts@Dimnames[[2]], scRNA@meta.data[["RNA_snn_res.0.25"]], sep="  ")
c = as.matrix(c)

write.table(c,"./sc_benchmark/temp.txt",row.names = FALSE, quote = FALSE)
lable <- read.table("./sc_benchmark/temp.txt", quote="\"", comment.char="")


#根据样本名合并标签
m0 <- merge(lable_ture, lable, by.x = "X", by.y = "row.names")
output = m0[,c(2,3)]
write.table(output,"./sc_benchmark/temp.txt", row.names=F, quote=FALSE, col.names = FALSE)

#轮廓系数模块
test=scRNA@reductions[["umap"]]@cell.embeddings
dist_matrix <- dist(test)

library(fpc)
stats <- cluster.stats(dist_matrix, lable$V1)
silhouette_coefficient <- stats$avg.silwidth

output_txt="./sc_benchmark/output/step2/imputation_output.txt"
write.table("-----Rmagic_imputation-----",output_txt, row.names=F, quote=FALSE, col.names = FALSE, append = TRUE)
write.table("运行时间为(秒):",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = " ")
write.table(execution_time[1],output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("silhouette_coefficient:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = " ")
write.table(silhouette_coefficient,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)

source("./sc_benchmark/scripts/DE_evaluation.R")