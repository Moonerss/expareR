#' @title Annotating gene expression matrix
#' @description Annotating gene expression matrix and remove duplicated genes
#' @name anno_expr
#' @param exprMat gene expression matrix with probe id as rownames and sample as columns.
#' @param anno annotation file with gene symbol and other ids matched to rowname of expression matrix,
#'    at least contain two column of probe id and symbol id.
#' @param symbol column name of gene symbol in annotation file
#' @param probe column name of relevant id in annotation file
#' @importFrom cli cli_alert_info
#' @importFrom dplyr filter left_jion select
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
#' @example
#'
anno_expr <- function(exprMat, anno, symbol = "symbol", probe = "probe_id", method = c("mean", "median")) {

  ## get annotation file
  colnames(anno)[which(colnames(anno)==symbol)] <- "symbol"
  colnames(anno)[which(colnames(anno)==probe)] <- "probe_id"
  anno <- anno[,c("probe_id","symbol")]
  anno <- anno[!anno$symbol=="NA_NA",]
  anno <- anno[!is.na(anno$symbol),]

  ## annotate probe
  anno_count <- length(intersect(rownames(exprMat), anno$probe_id))/nrow(exprMat)
  cli::cli_alert_info(paste0(paste0(sprintf(">>> %1.2f%%", 100*anno_count)," of probe in expression set was annotated")))

  inter_probe <- intersect(anno$probe_id, rownames(exprMat))
  anno <- dplyr::filter(anno, is.element(inter_probe, probe_id))
  exprMat <- exprMat[inter_probe,]
  exprMat <- exprMat %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "probe_id") %>%
    dplyr::left_join(anno, by = "probe_id") %>%
    dplyr::select(!probe_id) %>%
    tibble::column_to_rownames(var = "symbol") %>%
    as.matrix()

  ## remove duplicate
  res <- remove_duplicate(exprMat = exprMat, method = match.arg(method))
  return(res)
}
