## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  install.packages("GeoTcgaData")

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  if(!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  devtools::install_github("huerqiang/GeoTcgaData")

## -----------------------------------------------------------------------------
library(GeoTcgaData)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  library(TCGAbiolinks)
#  query <- GDCquery(project = "TCGA-ACC",
#                    data.category = "DNA Methylation",
#                    data.type = "Methylation Beta Value",
#                    platform = "Illumina Human Methylation 450")
#  GDCdownload(query, method = "api", files.per.chunk = 5, directory = Your_Path)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  merge_result <- Merge_methy_tcga(Your_Path_to_DNA_Methylation_data)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  library(ChAMP)
#  diff_gene <- methyDiff(cpgData = merge_result, sampleGroup = sample(c("C","T"),
#      ncol(merge_result[[1]]), replace = TRUE))

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  diff_gene$p.adj <- p.adjust(diff_gene$pvalue)
#  genes <- diff_gene[diff_gene$p.adj < 0.05, "gene"]
#  library(clusterProfiler)
#  library(enrichplot)
#  library(org.Hs.eg.db)
#  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
#  dotplot(ego)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  library(TCGAbiolinks)
#  query <- GDCquery(project = "TCGA-LGG",
#                    data.category = "Copy Number Variation",
#                    data.type = "Gene Level Copy Number Scores")
#  
#  GDCdownload(query, method = "api", files.per.chunk = 5, directory = Your_Path)
#  
#  data <- GDCprepare(query = query,
#                     directory =  "Your_Path")

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  class(data) <- "data.frame"
#  cnvData <- data[, -c(1,2,3)]
#  rownames(cnvData) <- data[, 1]
#  sampleGroup  = sample(c("A","B"), ncol(cnvData), replace = TRUE)
#  diffCnv <- diff_CNV(cnvData, sampleGroup)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  pvalues <- diffCnv$pvalue * sign(diffCnv$odds)
#  genes <- rownames(diffCnv)[diffCnv$pvalue < 0.05]
#  library(clusterProfiler)
#  library(enrichplot)
#  library(org.Hs.eg.db)
#  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
#  dotplot(ego)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  library(TCGAbiolinks)
#  query <- GDCquery(project = "TCGA-ACC",
#                    data.category = "Simple Nucleotide Variation",
#                    data.type = "Masked Somatic Mutation",
#                    workflow.type = "MuSE Variant Aggregation and Masking")
#  
#  GDCdownload(query, method = "api", files.per.chunk = 5, directory = Your_Path)
#  
#  data_snp <- GDCprepare(query = query,
#                     directory =  "Your_Path")

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  samples <- unique(data_snp$Tumor_Sample_Barcode)
#  sampleType <- sample(c("A","B"), length(samples), replace = TRUE)
#  names(sampleType) <- samples
#  pvalue <- diff_SNP_tcga(snpData = data_snp, sampleType = sampleType)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  pvalue2 <- sort(pvalue, decreasing = TRUE)
#  library(clusterProfiler)
#  library(enrichplot)
#  library(org.Hs.eg.db)
#  gsego <- gseGO(pvalue2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
#  dotplot(gsego)

## ---- message=FALSE, warning=FALSE--------------------------------------------
aa <- c("MARCH1","MARC1","MARCH1","MARCH1","MARCH1")
bb <- c(2.969058399,4.722410064,8.165514853,8.24243893,8.60815086)
cc <- c(3.969058399,5.722410064,7.165514853,6.24243893,7.60815086)
file_gene_ave <- data.frame(aa=aa,bb=bb,cc=cc)
colnames(file_gene_ave) <- c("Gene", "GSM1629982", "GSM1629983")
result <- gene_ave(file_gene_ave, 1)

## -----------------------------------------------------------------------------
aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3",
        "MARCH3 /// MARCH4","MARCH1")
bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
input_file <- data.frame(aa=aa,bb=bb,cc=cc)
rep1_result <- rep1(input_file," /// ")
rep2_result <- rep2(input_file," /// ")

## -----------------------------------------------------------------------------
id_conversion_vector("symbol", "ensembl_gene_id", c("A2ML1", "A2ML1-AS1", "A4GALT", "A12M1", "AAAS")) 

## -----------------------------------------------------------------------------
result <- id_conversion(profile)

## -----------------------------------------------------------------------------
library(clusterProfiler)
bitr(c("A2ML1", "A2ML1-AS1", "A4GALT", "A12M1", "AAAS"), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db, drop = FALSE)

## -----------------------------------------------------------------------------
lung_squ_count2 <- matrix(c(1,2,3,4,5,6,7,8,9),ncol=3)
rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
jieguo <- countToFpkm_matrix(lung_squ_count2)

## ---- message=FALSE, warning=FALSE--------------------------------------------
lung_squ_count2 <- matrix(c(0.11,0.22,0.43,0.14,0.875,0.66,0.77,0.18,0.29),ncol=3)
rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
jieguo <- countToTpm_matrix(lung_squ_count2)

## ---- message=FALSE, warning=FALSE--------------------------------------------
tcga_cli <- tcga_cli_deal(system.file(file.path("extdata","tcga_cli"),package="GeoTcgaData"))

