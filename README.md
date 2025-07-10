# Single-Cell-Sequencing-Analysis-Tutorial

### 起始文件：barcodes.tsv, features.tsv, matrix.mtx

### 1. 数据质量控制：
- 在R Studio (v4.4.3) 里进行三个起始文件以及metadata注释文件的读入
- 为count matrix添加行名与列名之后 把数据整理到Seurat object里
  rna <- CreateSeuratObject(counts = mtx, assay = "RNA", meta.data = metadata)

- 查看数据的独特基因量nFeatures UMI量nCounts 属线粒体的基因量percent.mt
  VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

- 将不合格的细胞剔除
  rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10)
  基因数过少可能是空液滴或低质量细胞 基因数过多可能是一个油滴不小心捕获了多个细胞
  而且线粒体基因高表达通常表示细胞死亡或损伤
  所以应针对基因数在200到3500之间 线粒体基因在10%以内的细胞进行后续分析

### 2. 细胞亚群识别与聚类：
- 对过滤后的数据进行标准化（可以很好地消除细胞之间测序深度的差异）
  rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

- 再对数据进行线性变换 保证各基因在所有细胞中的表达均数=0 方差=1（引入不同处理组间的可比性）
  all.genes <- rownames(rna)
  rna <- ScaleData(rna, features = all.genes)

- 挑选出2000个高变基因（在不同细胞之间表达差异较大的基因）用于减少数据维度并提高计算效率
  rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

- 用这些高变基因进行主成分分析 识别数据的主要变异来源
  rna <- RunPCA(rna, features = VariableFeatures(object = rna))

- 根据细胞间高维距离进行细胞亚群聚类（PCA结果可视化后根据拐点我们选了10个PC进行后续聚类）
  rna <- FindNeighbors(rna, dims = 1:10)
  rna <- FindClusters(rna, resolution = 0.5)

### 3. 差异表达分析（不同处理组之间细胞类型分布的差异）：
- 根据Seurat里面Azimuth工具生成的分类 结合FeaturePlot高表达基因 对簇的细胞种类进行标注
  rna <- RunAzimuth(rna, reference = "pbmcref")

  rna_1G <- subset(rna, subset = treatment == "1G")  # 正常重力
  FeaturePlot(rna_1G, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A7", "GNLY"))
  rna_uG <- subset(rna, subset = treatment == "uG")  # 微重力
  FeaturePlot(rna_uG, features = c("S100A4", "MS4A1", "CD8A", "FCGR3A", "CX3CR1"))
  
  new.cluster.ids <- c("CD4 Naive", "CD4 TEM", "CD8 TCM", "CD8 TEM", "NK", "CD14 Mono",
                     "cDC2", "B Naive", "CD16 Mono", "cDC1", "B Memory", "CD4 CTL",
                     "CD4 TCM", "Eryth", "Treg", "pDC")[1:length(current_clusters)]  
  names(new.cluster.ids) <- current_clusters
  rna <- RenameIdents(rna, new.cluster.ids)  # 对各簇进行细胞种类标注

- 用UMAP投影并排显示不同处理组（微重力-正常重力）的细胞类型分布
  DimPlot(rna, split.by = "treatment", label = TRUE, repel = TRUE, reduction = "umap")

### 4. 差异表达分析（不同处理组之间细胞类型分布的差异）：
- 用Wilcox模型对各类别细胞不同处理组（微重力-正常重力）之间的差异表达基因进行筛选
  rna.celltype5 <- subset(rna, seurat_clusters == 5)
  rna.marker.treat5 <- FindMarkers(rna.celltype5, group.by = "treatment", 
                                   ident.1 = "uG", ident.2 = "1G"))

- 在各细胞簇里绘制火山图 直观展示同类细胞不同处理组间的多基因表达差异（上调/下调）
  EnhancedVolcano(rna.marker.treat5.notaxid,
                  lab = rownames(rna.marker.treat5.notaxid),
                  x = 'avg_log2FC', y = 'p_val_adj', xlim =c(-1.5,2),
                  FCcutoff = 0.25, labSize = 5, 
                  drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                  title = "uG vs 1G", subtitle = bquote(italic("unstimulated PBMC")))

- 在各细胞簇里查看差异表达的基因 用小提琴图直观展示同类细胞单个基因在不同处理组间差异表达的水平
  View(rna.marker.treat5)
  VlnPlot(rna.celltype5, features = "CCL2", split.by = "treatment", pt.size = 0)

### 5. 后续分析？
