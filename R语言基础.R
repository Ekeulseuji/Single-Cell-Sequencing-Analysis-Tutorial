library(DESeq2)

download.file(url="https://comics.med.sustech.edu.cn/summer2025/human_read_counts.csv", destfile="hrc.csv")

# count table
hrc_dat <- read.csv("hrc.csv", row.names=1)
head(hrc_dat)

# coldata
col_data <- data.frame(
  row.names = colnames(hrc_dat),
  tissue = c(rep('A',10),rep('B',10)))

# DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = hrc_dat,
                              colData = col_data,
                              design = ~ tissue)
dds <- DESeq(dds)

# PCA 
vsd <- vst(dds)
plotPCA(vsd, intgroup='tissue')

# results
res <- results(dds, lfcThreshold = 1)
summary(res)
