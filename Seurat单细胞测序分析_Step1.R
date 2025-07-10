library(Seurat)
library(dplyr)
library(DESeq2)

################################################################################
################################################################################

# 读取数据
path <- "/Users/ekeulseuji/Downloads/PBMC_Original/"
mtx <- Matrix::readMM(file = paste0(path, "matrix.mtx.gz"))
features <- read.table(file = paste0(path, "features.tsv.gz"))
barcodes <- read.table(file = paste0(path, "barcodes.tsv.gz"))
metadata <- read.table(file = paste0(path, "metadata.tsv.gz"), header = TRUE)

# 整理到Seurat object里面
rownames(mtx) <- make.unique(features$V1)
colnames(mtx) <- barcodes$V1
rownames(metadata) <- metadata$X
rna <- CreateSeuratObject(counts = mtx, assay = "RNA", meta.data = metadata)

# 检查是否关联上了metadata
head(rna@meta.data)  

# 查看数据 独特基因量 UMI量 线粒体量
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
# rna <- subset(rna, subset = nFeature_RNA > 200
#               & nFeature_RNA < 2500
#               & percent.mt < 5)

VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

################################################################################
################################################################################

# rna_sct <- rna  
# rna_sct <- SCTransform(rna_sct, verbose = FALSE)
# ------------------------------------------------------------- SCT
# 安装依赖包（若未安装）
# if (!require("sctransform")) install.packages("sctransform")

# 一步完成标准化 高变基因选择和缩放
# rna <- SCTransform(
#   rna,
#   assay = "RNA",          # 原始数据所在的assay
#   new.assay.name = "SCT", # 结果存储到SCT assay
#   vars.to.regress = NULL, # 可选（校正批次/线粒体含量等）
#   verbose = FALSE
# )

# ------------------------------------------------------------- 传统
# # Normalize 消除细胞间测序深度差异
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# # 挑出表达有差异的2000个基因
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
# 
# # Scale线性变换 使基因在所有细胞中的平均表达为0+方差为1 增加可比性
# # all.genes <- rownames(rna)
rna <- ScaleData(rna) # features = all.genes
# 
# # top 10 高变基因 减少数据维度 提高计算效率
# # top10 <- head(VariableFeatures(rna), 10)
# # # plot variable features with (p1) and without labels (p2)
# # plot1 <- VariableFeaturePlot(rna)
# # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



################################################################################
################################################################################

# PCA识别数据的主要变异来源 PC1和PC2通常代表最重要的生物学信号
rna <- RunPCA(rna, features = VariableFeatures(object = rna))  # 传统
print(rna[["pca"]], dims = 1:5, nfeatures = 5)
# rna <- RunPCA(rna, assay = "SCT", npcs = 50)  # SCT 指定assay

# 查看主成分 选择变异来源数量
VizDimLoadings(rna, dims = 1:2, reduction = "pca")
DimPlot(rna, reduction = "pca") + NoLegend()

DimHeatmap(rna, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(rna, dims = 1:15, cells = 500, balanced = TRUE)

# 找拐点确定数量
ElbowPlot(rna)


################################################################################
################################################################################

# Louvain聚类 选了10个主成分
rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna, resolution = 0.5)
# head(Idents(rna), 4)  # cluster ID分在Idents里面

# UMAP可视化聚类
rna <- RunUMAP(rna, dims = 1:10)
# DimPlot(rna, reduction = "umap")
# 原数据有两个treatment所以分开展示聚类结果
DimPlot(rna, split.by = "treatment", label = TRUE, repel = TRUE)


################################################################################
################################################################################

# 找出各个Cluster相对来说的高表达基因 标注细胞簇的类型
# find markers for every cluster compared to all remaining cells
# report only the positive ones
rna.markers <- FindAllMarkers(rna, only.pos = TRUE)
rna.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)  

# 精细过滤
rna.markers.filtered <- rna.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  mutate(Cell_type = Idents(rna)[cluster] %>% as.character()) # 添加细胞类型注释

# 用Azimuth标注细胞簇的类型 一般作为参考 注意要是SCT Assay
library(Azimuth)
# rna <- RunAzimuth(rna, reference = "pbmcref")
# 可能跑不动 那就下载原count matrix之后去网站标注
# https://app.azimuth.hubmapconsortium.org/app/human-pbmc
# tmp <- LayerData(rna, assay = "RNA", layer = "counts")
# saveRDS(tmp, file="rna_count.rds")

# 查看各cluster的top标记基因
top_markers <- rna.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  pull(gene)

print(top_markers)

rna_1G <- subset(rna, subset = treatment == "1G")  # 分类出正常重力
rna_uG <- subset(rna, subset = treatment == "uG")  # 和微重力

# 看经典PBMC基因对应高表达的Cluster
FeaturePlot(rna_1G, features = c("IL7R", "CCR7", "CD14", "LYZ", "IL7R", 
                              "S100A4", "MS4A1", "CD8A", "FCGR3A", "CX3CR1", 
                              "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", 
                              "PPBP", "IL3RA", "CD3E", "HBA1", "HBB", 
                              "CLEC4C", "GZMA", "CD1C", "SIRPA", "CD44"), 
            cols = c("lightgrey", "blue"), pt.size = 1)

FeaturePlot(rna_uG, features = c("IL7R", "CCR7", "CD14", "LYZ", "IL7R", 
                                 "S100A4", "MS4A1", "CD8A", "FCGR3A", "CX3CR1", 
                                 "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", 
                                 "PPBP", "IL3RA", "CD3E", "HBA1", "HBB", 
                                 "CLEC4C", "GZMA", "CD1C", "SIRPA", "CD44"),
            cols = c("lightgrey", "red"), pt.size = 1)

# 获取当前所有聚类ID
current_clusters <- levels(Idents(rna))

# 定义新细胞类别标签（数量必须与目前的clusters数量一致）
new.cluster.ids <- c("CD4 Naive", "CD4 TEM", "CD8 TCM", "CD8 TEM", "NK", "CD14 Mono", 
                     "cDC2", "B Naive", "CD16 Mono", "cDC1", "B Memory", "CD4 CTL",
                     "CD4 TCM", "Eryth", "Treg", "pDC")[1:length(current_clusters)]

# 命名必须严格匹配数量
names(new.cluster.ids) <- current_clusters
# 执行重命名
rna <- RenameIdents(rna, new.cluster.ids)

# 查看注释结果
DimPlot(rna, split.by = "treatment", label = TRUE, repel = TRUE, reduction = "umap")
# 单独查看基因
# VlnPlot(rna, features = c("MS4A1", "CD79A"), pt.size = 0)

# 手动类型注释存入 predicted.celltype.l2 
rna$predicted.celltype.l2 <- Idents(rna) 


################################################################################
################################################################################

# 调整p值 除掉零值p_val_adj
rna.markers.filtered$p_val_adj[rna.markers.filtered$p_val_adj == 0] <- sort(
  unique(rna.markers.filtered$p_val_adj))[2]

# 根据uG和1G之间的差别看基因表达差异
# library(MAST)
# library(DESeq2)
rna.marker.treat <- FindMarkers(rna, group.by = "treatment",
                                ident.1 = "uG", ident.2 = "1G", 
                                min.pct=0, logfc.threshold = 0,
                                test.use = "wilcox")

# 根据log2FC大小排列
rna.marker.treat <- rna.marker.treat[order(rna.marker.treat$avg_log2FC, decreasing = T),] 


# 在差异表达基因（DEGs）结果中过滤掉那些可能不感兴趣的基因
rna.marker.treat.notaxid<-rna.marker.treat[grep("taxid|MT-", row.names(rna.marker.treat),
                                                invert = T),]
library(EnhancedVolcano)

EnhancedVolcano(rna.marker.treat.notaxid,
                lab = rownames(rna.marker.treat.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim =c(-1.5,2),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "uG vs 1G",
                subtitle = bquote(italic("unstimulated PBMC")))

ggsave("uGvs1G_unstimulated_Volcano.pdf", width=10, height=10)


library(ggplot2)
library(scales)

# 选择要展示的基因（例如top 20）
top_genes <- rna.markers.filtered %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique() %>%
  head(20)

# 绘制点图
rna.markers.filtered %>%
  filter(gene %in% top_genes) %>%
  ggplot(aes(
    x = gene, 
    y = factor(Cell_type, levels = unique(Cell_type)),  # 保持细胞类型顺序
    color = avg_log2FC, 
    size = -log10(p_val_adj))
  ) +
  geom_point() +
  scale_size_area(
    max_size = 5,
    breaks = c(0, 25, 50, 100),
    labels = c("0", "25", "50", "100+"),
    guide = "legend"
  ) +
  colorspace::scale_color_continuous_divergingx(
    palette = 'RdBu',
    rev = TRUE,
    limits = c(-1, 1),
    name = 'Log2FC'
  ) +
  labs(
    x = "Genes",
    y = "Cell Types",
    title = "Top Differential Genes by Cell Type"
  ) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10)
  )

# ggsave("DotPlot_Markers_by_CellType.pdf", width = 11, height = 6)


# 调整p值 除掉零值p_val_adj
rna.markers.filtered$p_val_adj[rna.markers.filtered$p_val_adj == 0] <- sort(
  unique(rna.markers.filtered$p_val_adj))[2]

rna.celltype0 <- subset(rna, seurat_clusters == 0)
rna.marker.treat0 <- FindMarkers(rna.celltype0, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")


# 根据uG和1G之间的差别看基因表达差异
# library(MAST)
# library(DESeq2)

# 根据log2FC大小排列
rna.marker.treat0.sorted <- rna.marker.treat0[order(rna.marker.treat0$avg_log2FC, decreasing = T),] 


# 在差异表达基因（DEGs）结果中过滤掉那些可能不感兴趣的基因
rna.marker.treat0.notaxid<-rna.marker.treat0.sorted[grep("taxid|MT-", 
                                                         row.names(rna.marker.treat0.sorted),
                                                         invert = T),]
library(EnhancedVolcano)

# 在各细胞簇里绘制火山图，直观展示同类细胞不同处理组间的多基因表达差异（上调/下调）
EnhancedVolcano(rna.marker.treat0.notaxid,
                lab = rownames(rna.marker.treat0.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim =c(-1.5,2),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "uG vs 1G",
                subtitle = bquote(italic("unstimulated PBMC")))

ggsave("uGvs1G_unstimulated_Volcano.pdf", width=10, height=10)



# 各个Cluster里面做细分看基因表达量对比
rna.celltype1 <- subset(rna, seurat_clusters == 1)
rna.marker.treat1 <- FindMarkers(rna.celltype1, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype2 <- subset(rna, seurat_clusters == 2)
rna.marker.treat2 <- FindMarkers(rna.celltype2, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype3 <- subset(rna, seurat_clusters == 3)
rna.marker.treat3 <- FindMarkers(rna.celltype3, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype4 <- subset(rna, seurat_clusters == 4)
rna.marker.treat4 <- FindMarkers(rna.celltype4, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype5 <- subset(rna, seurat_clusters == 5)
rna.marker.treat5 <- FindMarkers(rna.celltype5, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype6 <- subset(rna, seurat_clusters == 6)
rna.marker.treat6 <- FindMarkers(rna.celltype6, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype7 <- subset(rna, seurat_clusters == 7)
rna.marker.treat7 <- FindMarkers(rna.celltype7, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype8 <- subset(rna, seurat_clusters == 8)
rna.marker.treat8 <- FindMarkers(rna.celltype8, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype9 <- subset(rna, seurat_clusters == 9)
rna.marker.treat9 <- FindMarkers(rna.celltype9, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype10 <- subset(rna, seurat_clusters == 10)
rna.marker.treat10 <- FindMarkers(rna.celltype10, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype11 <- subset(rna, seurat_clusters == 11)
rna.marker.treat11 <- FindMarkers(rna.celltype11, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype12 <- subset(rna, seurat_clusters == 12)
rna.marker.treat12 <- FindMarkers(rna.celltype12, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype13 <- subset(rna, seurat_clusters == 13)
rna.marker.treat13 <- FindMarkers(rna.celltype13, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype14 <- subset(rna, seurat_clusters == 14)
rna.marker.treat14 <- FindMarkers(rna.celltype14, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")
rna.celltype15 <- subset(rna, seurat_clusters == 15)
rna.marker.treat15 <- FindMarkers(rna.celltype15, group.by = "treatment", 
                                 ident.1 = "uG", ident.2 = "1G")

# 在各细胞簇里查看差异表达的基因
# 用小提琴图直观展示同类细胞基因在不同处理组间差异表达的水平
View(rna.marker.treat5)
VlnPlot(rna.celltype5, features = "CCL2", split.by = "treatment", pt.size = 0)




