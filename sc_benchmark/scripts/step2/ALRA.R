source("R_conf.R")
system("bash ./sc_benchmark/scripts/kill_user_usage.sh")

library(ALRA)
load(datapath)

system("bash ./sc_benchmark/scripts/start_monitor.sh -n ./sc_benchmark/output/step2/ALRA/usage -t 5")
start_time <- proc.time()

b_nk_example=as.matrix(dat)
b_nk_example=t(b_nk_example)

A_norm <- normalize_data(b_nk_example)

k_choice <- choose_k(A_norm)

A_norm_completed <- alra(A_norm,k=k_choice$k)[[3]]
A_norm_completed=t(A_norm_completed)
colnames(A_norm_completed) <- colnames(dat)
dat=A_norm_completed

end_time <- proc.time()
execution_time <- end_time - start_time
system("bash ./sc_benchmark/scripts/kill_user_usage.sh")

save(dat, file = "./sc_benchmark/output/step2/ALRA/ALRA_data.RData")
load("./sc_benchmark/output/step2/ALRA/ALRA_data.RData")

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

load(lablepath)

scRNA <- CreateSeuratObject(counts = dat, project = "Grubman", min.cells = 1, min.features = 1); scRNA

#使用从MT-作为线粒体基因集
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "MT_")
#qc指标可视化
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & nCount_RNA > 500)

scRNA <- ScaleData(scRNA)
#scRNA <- NormalizeData(scRNA)
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

png("./sc_benchmark/output/step2/ALRA/ALRA_umap.png", height=768, width=768)
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
write.table("-----ALRA_imputation-----",output_txt, row.names=F, quote=FALSE, col.names = FALSE, append = TRUE)
write.table("运行时间为(秒):",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = " ")
write.table(execution_time[1],output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)
write.table("silhouette_coefficient:",output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE,eol = " ")
write.table(silhouette_coefficient,output_txt, row.names=F, quote=FALSE, col.names = FALSE,append = TRUE)

source("./sc_benchmark/scripts/DE_evaluation.R")