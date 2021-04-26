#' @title Remove duplicate gene symbol
#' @description Remove duplicate gene symbol on a gene expression data
#' @name remove_duplicate
#' @param exprMat A data frame with geneid as as first column and other columns are expression of gene in each samples.
#' @param gene_col The column name or number of geneid in \code{exprMat}. Dafult is NULL, choose the first column as geneid.
#' @param method method used to filter duplicate genes; filter duplicated gene by the mean of median value, or combine the value by mean or median value instead
#' @param value use mean of median value as reference
#' @details We use the mean or median value of the duplicated genes as the expression value of this gene.
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr group_by ungroup mutate pull
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom stats median
#' @export
#' @return return a expression matrix
#'
remove_duplicate <- function(exprMat, gene_col = NULL, method = c("order", "combine"), value = c("mean", "median")) {

  ## get gene id
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

  ## check whether have duplicate
  dup <- nrow(exprMat) - length(unique(genes))
  if (dup == 0) {
    exprMat <- exprMat %>% column_to_rownames(var = "gene_id") %>% as.matrix()
    return(exprMat)
  } else {
    meta_expr <- exprMat %>%
      group_by(gene_id) %>% nest() %>% ungroup()

    method <-  match.arg(method)
    value <- match.arg(value)
    # filter genes
    if (method == "order") {
      if (value == "mean") {
        meta_expr <- meta_expr %>% mutate(expr = purrr::map(.data$data, function(x) {
          order_x <- x[order(apply(x, 1, mean, na.rm = T), decreasing = T),]
          order_x <- order_x[1,]
        }))
      } else if (value == "median") {
        meta_expr <- meta_expr %>% mutate(expr = purrr::map(.data$data, function(x) {
          order_x <- x[order(apply(x, 1, median, na.rm = T), decreasing = T),]
          order_x <- order_x[1,]
        }))
      }
    } else if (method == "combine") {
      if (value == "mean") {
        meta_expr <- meta_expr %>% mutate(expr = purrr::map(data, ~ apply(.x, 2, mean, na.rm = T)))
      } else if (value == "median") {
        meta_expr <- meta_expr %>% mutate(expr = purrr::map(data, ~ apply(.x, 2, median, na.rm = T)))
      }
    }

    res <- cbind.data.frame(genes = meta_expr$gene_id, do.call(rbind, meta_expr$expr))
    return(res)
  }
}
