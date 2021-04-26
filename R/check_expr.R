#' @title Check the data quality of expression matrix
#' @description Check the expression matrix whether need other process
#' @name check_expr
#' @param data a gene or probe data frame like \code{arrayMat} in \code{\link{anno_array}} or \code{exprMat} in \code{\link{remove_duplicate}}
#' @param id_col The column name or number of geneid or probe id. Dafult is NULL, choose the first column as id.
#' @param verbose output other useful information, default TRUE.
#' @details The expression matrix will be checked for three task:
#'    \itemize{
#'      \item check whether have NA
#'      \item check whether have Inf or -Inf
#'      \item check the standard deviation(sd) whether is 0
#'      \item check whether have duplicated genes
#'   }
#'   if have problem, return FALSE, else no problem return TRUE.
#' @return return a logical value.
#' @export

check_expr <- function(data, id_col = NULL, verbose = TRUE) {

  ## get probe in expression data
  if (is.null(id_col)) {
    message("`id_col = NULL`, choose first column as id.")
    id_col <- 1
    colnames(data)[1] <- "ids"
  } else if (is.character(id_col)) {
    colnames(data)[which(colnames(data) == id_col)] <- "ids"
  } else if (is.numeric(id_col)) {
    colnames(data)[id_col] <- "ids"
  }
  exprMat <- data %>% select(!ids) %>% as.matrix()
  all_id <- data[, "ids"]

  ## check whether have NA
  if (verbose) cli::cli_alert_info("Checking whether have `NA` ...")
  na_res <- ifelse(sum(is.na(exprMat)) > 0, FALSE, TRUE)
  if (verbose) cli::cli_alert_success("Done!")

  ## check whether have Inf ot -Inf
  if (verbose) cli::cli_alert_info("Checking whether have `Inf` or `-Inf` ...")
  inf_res <- ifelse(min(exprMat, na.rm = T) == -Inf | max(exprMat, na.rm = T)== Inf, FALSE, TRUE)
  if (verbose) cli::cli_alert_success("Done!")

  ## check whether have duplicated genes
  if (verbose) cli::cli_alert_info("Checking whether have duplicated genes ...")
  id_res <- ifelse(length(all_id) != length(unique(all_id)), FALSE, TRUE)
  if (verbose) cli::cli_alert_success("Done!")

  ## check standard deviation(sd)
  if (verbose) cli::cli_alert_info("Checking the standard deviation ...")
  sd <- apply(exprMat, 1, function(x) sd(x, na.rm = T) == 0)
  sd_res <- ifelse(nlevels(as.factor(sd))>1, FALSE, TRUE)
  if (verbose) cli::cli_alert_success("Done!")

  ## get the result
  if (all(c(na_res, inf_res, sd_res))) {
    if (verbose) cli::cli_alert_info("No problem~~~~")
    return(TRUE)
  } else {
    cli::cli_alert_warning(col_yellow("Have some problem:"))
    lid <- cli_ol()
    if (!na_res) cli_li(("Have `NA` in matrix"))
    if (!inf_res) cli_li("Have `Inf` or `-Inf` in matrix")
    if (!id_res) cli_li("Have duplicated genes")
    if (!sd_res) cli_li("Some vairables in matrix have no variance between samples")
    cli_end(lid)
    return(FALSE)
  }
}
