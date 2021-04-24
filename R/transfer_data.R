#' @title Transform the expression matrix to normal quality
#' @description Transform the expression matrix  with bad quality to normal quality
#' @name transfer_data
#' @param exprMat the expression matrix
#' @param data_type the data type to change, NA, Inf or both
#' @param data_to transform data into 0 or 1, or others, default is 0
#' @return return minipulated data
#' @export
#' @example
#'
transfer_data <- function(exprMat, data_type = c("all", "NA", "Inf"), data_to = 0) {
  data_type <- match.arg(data_type)
  if (data_type == "all") {
    res <- apply(exprMat, 2, function(x) {ifelse(is.na(x) | is.infinite(x), data_to, x)})
  } else if (data_type == "NA") {
    res <- apply(exprMat, 2, function(x) {ifelse(is.na(x), data_to, x)})
  } else if (data_type == "Inf") {
    res <- apply(exprMat, 2, function(x) {ifelse(is.infinite(x), data_to, x)})
  }
  return(res)
}
