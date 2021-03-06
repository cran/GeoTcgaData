#' Get the differentially expressioned genes using DESeq2 package
#'
#' @param profile_input a data.frame
#'
#' @return a data.frame, a intermediate results of DESeq2
#' @export
#'
#' @examples
#' profile2 <- classify_sample(kegg_liver)
classify_sample <- function(profile_input){
    biaoqian<-profile_input[1,]
    dada<-biaoqian
    for(i in 2:length(biaoqian)){
      dada[i]<-unlist(strsplit(biaoqian[i],"-"))[4]
      dada[i]<-unlist(strsplit(dada[i],"_"))[1]
    }
  
    file2<-profile_input
    file2[1,]<-dada
    return(file2)
}


#' Get the differentially expressioned genes using DESeq2 package
#'
#' @param profile2_input a result of classify_sample
#'
#' @return a matrix, information of differential expression genes
#' @export
#'
#' @examples
#' profile2 <- classify_sample(kegg_liver)
#' jieguo <- diff_gene(profile2)
diff_gene<-function(profile2_input){
  if(requireNamespace("DESeq2", quielty = TRUE)) {
  database <- profile2_input[-1,-1]

  database <- matrix(as.numeric(database),nrow=nrow(database))
  rownames(database) <- profile2_input[-1, 1]
  condition <- profile2_input[1, -1]
  database <- round(as.matrix(database))
  rownames(database) <- as.character(rownames(database))
  condition <- as.numeric(condition)
  condition <- round(condition)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData=database, 
      colData=S4Vectors::DataFrame(condition), design=~condition)

  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  table(res$padj <0.05)
  res <- res[order(res$padj),]
  return(res)
} else { message("please install the package DESeq2!")}
}