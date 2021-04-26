#' @title Transform the expression matrix to normal quality
#' @description Transform the expression matrix  with bad quality to normal quality
#' @name transfer_data
#' @param data a gene or probe data frame like \code{arrayMat} in \code{\link{anno_array}} or \code{exprMat} in \code{\link{remove_duplicate}}
#' @param gene_col The column name or number of geneid in \code{exprMat}. Dafult is NULL, choose the first column as geneid.
#' @param transfer_type the data type to change, NA, Inf or all
#' @param data_to transform data into 0 or 1, or others, default is 0
#' @return return minipulated data
#' @export
#'
transfer_data <- function(exprMat, gene_col = NULL, data_type = c("all", "NA", "Inf"), data_to = 0) {
  data_type <- match.arg(data_type)

  ## get genes
  if (is.null(gene_col)) {
    message("`gene_col = NULL`, choose first column as gene id.")
    gene_col <- 1
    colnames(exprMat)[1] <- "gene_id"
  } else if (is.character(gene_col)) {
    colnames(exprMat)[which(colnames(exprMat) == gene_col)] <- "gene_id"
  } else if (is.numeric(gene_col)) {
    colnames(exprMat)[gene_col] <- "gene_id"
  }
  genes <- exprMat[, "gene_id"]
  exprMat <- exprMat %>% select(!gene_id) %>% as.matrix()

  if (data_type == "all") {
    res <- apply(exprMat, 2, function(x) {ifelse(is.na(x) | is.infinite(x), data_to, x)})
  } else if (data_type == "NA") {
    res <- apply(exprMat, 2, function(x) {ifelse(is.na(x), data_to, x)})
  } else if (data_type == "Inf") {
    res <- apply(exprMat, 2, function(x) {ifelse(is.infinite(x), data_to, x)})
  }

  res <- cbind.data.frame(genes, res)
  return(res)
}
