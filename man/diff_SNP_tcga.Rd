% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNP.R
\name{diff_SNP_tcga}
\alias{diff_SNP_tcga}
\title{Do difference analysis of SNP data downloaded from TCGAbiolinks}
\usage{
diff_SNP_tcga(snpData, sampleType)
}
\arguments{
\item{snpData}{data.frame of SNP data downloaded from TCGAbiolinks}

\item{sampleType}{vector of sample group}
}
\description{
Do difference analysis of SNP data downloaded from TCGAbiolinks
}
\examples{
\dontrun{
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-ACC",
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation",
                  workflow.type = "MuSE Variant Aggregation and Masking")

GDCdownload(query, method = "api", files.per.chunk = 5, directory = Your_Path)

data_snp <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  "Your_Path") 
samples <- unique(data_snp$Tumor_Sample_Barcode)
sampleType <- sample(c("A","B"), length(samples), replace = TRUE)
names(sampleType) <- samples
pvalue <- diff_SNP_tcga(snpData = data_snp, sampleType = sampleType)
}
}
