#' @title Annotating array probe expression matrix
#' @description Annotating array probe expression matrix and remove duplicated genes
#' @name anno_expr
#' @param arrayMat a data frame with probe id as first column and other columns are expression of gene in each samples
#' @param anno annotation file with gene symbol and other ids matched to probe id of \code{arrayMat},
#'    at least contain two column of probe id and symbol id.
#' @param probe_col The column name or number of probe id in \code{arrayMat}. Dafult is NULL, choose the first column as probe id.
#' @param symbol column name of gene symbol in annotation file
#' @param probe column name of relevant id in annotation file
#' @importFrom cli cli_alert_info
#' @importFrom dplyr filter left_join select pull
#' @importFrom tibble rownames_to_column column_to_rownames
#' @return return a gene expression matrix
#' @export
#'
anno_expr <- function(arrayMat, anno, probe_col = NULL, symbol = "symbol", probe = "probe_id") {

  ## get annotation file
  colnames(anno)[which(colnames(anno)==symbol)] <- "symbol"
  colnames(anno)[which(colnames(anno)==probe)] <- "probe_id"
  anno <- anno[,c("probe_id","symbol")]
  anno <- anno[!anno$symbol=="NA_NA",]
  anno <- anno[!is.na(anno$symbol),]

  ## get probe in expression data
  if (is.null(probe_col)) {
    message("`probe_col = NULL`, choose first column as probe id.")
    probe_col <- 1
    colnames(arrayMat)[1] <- "probe_id"
  } else if (is.character(probe_col)) {
    colnames(arrayMat)[which(colnames(arrayMat) == probe_col)] <- "probe_id"
  } else if (is.numeric(probe_col)) {
    colnames(arrayMat)[probe_col] <- "probe_id"
  }
  probes <- arrayMat %>% pull(probe_id)

  ## annotate probe
  anno_count <- length(intersect(probes, anno$probe_id))/nrow(arrayMat)
  cli_alert_info(paste0(paste0(sprintf(">>> %1.2f%%", 100*anno_count)," of probe in expression set was annotated")))

  inter_probe <- intersect(anno$probe_id, probes)
  anno <- dplyr::filter(anno, is.element(probe_id, inter_probe))
  arrayMat <- arrayMat[match(pull(anno, probe_id), pull(arrayMat, probe_id)),]

  ## trans data
  exprMat <- arrayMat %>%
    dplyr::left_join(anno, by = "probe_id") %>%
    dplyr::select(!probe_id) %>%
    dplyr::select(symbol, !symbol)
  colnames(exprMat)[1] <- "genes"

  return(exprMat)
}
