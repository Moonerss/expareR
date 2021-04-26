#' @title log2 transformation of expression matrix
#' @description Do the expression matrix needs log2 transformation
#' @usage log2expr(exprMat)
#' @param exprMat gene expression matrix with row as genes and sample in the column
#' @param gene_col column name or number of gene id, default is NULL, choose the first column
#' @return a gene expression matrix after log2 transformation
#' @details The function will first check whether the expression matrix has undergone
#'   log2 transformation; if the expression matrix have done log2 transformation, return
#'   the raw expression matrix, else do the log2 transformation
#' @importFrom cli cli_alert_info
#' @name log2expr
#' @export
#'
#' @examples
#'   mat <- matrix(sample(1:10, 30, replace = T), ncol = 5)
#'   log2_mat <- log2expr(exprMat = mat)
#'   res <- log2expr(exprMat = log2_mat)
log2expr <- function(exprMat = NULL, gene_col = NULL){

  ## get gene expression
  if (is.data.frame(exprMat)) {
    if (is.null(gene_col)) {
      message("`gene_col = NULL`, choose first column as gene id.")
      gene_col <- 1
      colnames(exprMat)[1] <- "genes_id"
    } else if (is.character(gene_col)) {
      colnames(exprMat)[which(colnames(exprMat) == gene_col)] <- "genes_id"
    } else if (is.numeric(gene_col)) {
      colnames(exprMat)[gene_col] <- "genes_id"
    }
    exprMat <- exprMat %>% column_to_rownames(var = "genes_id") %>% as.matrix()
  }

  ## check whether done log2 transformation
  qx <- as.numeric(quantile(exprMat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  loged <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

  if (loged) {
    exprMat[exprMat < 0] <- 0
    exprMat <- log2(exprMat+1)
    cli::cli_alert_info("log2 transformation finished!")
  } else {
    cli::cli_alert_info("log2 transformation is not necessary.")
  }

  return(exprMat)
}
