#' @title Remove duplicate gene symbol
#' @description Remove duplicate gene symbol on a gene expression data
#' @name remove_duplicate
#' @param exprMat A expression matrix, with geneid as rownames and sample as columns.
#' @param method method used to filter duplicate genes; default is mean value
#' @details We use the mean or median value of the duplicated genes as the expression value of this gene.
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr group_by ungroup mutate pull
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @export
#' @return return a expression matrix
#' @examples
remove_duplicate <- function(exprMat, method = c("mean", "median")) {

  ## check whether have duplicate
  dup <- nrow(exprMat) - length(unique(rownames(exprMat)))
  if (dup == 0) {
    return(exprMat)
  } else {
    meta_expr <- exprMat %>% as.data.frame() %>%
      rownames_to_column(var = "genes") %>%
      group_by() %>% nest() %>% ungroup()

    method = match.arg(method)
    # filter genes
    if (method == "mean") {
      meta_expr <- meta_expr %>% mutate(median_expr = purrr::map(data, ~ apply(.x, 2, mean, na.rm = T)))
    } else if (method == "median") {
      meta_expr <- meta_expr %>% mutate(median_expr = purrr::map(data, ~ apply(.x, 2, median, na.rm = T)))
    }
    res <- do.call(rbind, meta_expr$median_expr) %>% as.matrix()
    rownames(res) <- meta_expr %>% pull(genes)
    return(res)
  }
}
