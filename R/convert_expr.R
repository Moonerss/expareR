#' @title Convert read counts to transcripts per million (TPM)
#' @description Convert the count matrix to TPM matrix
#' @name convert_expr
#' @param exprMat A data frame with geneid as as first column and other columns are expression of gene in each samples.
#' @param idType Type of gene id. "ensembl" "entrez","symbol"
#' @param gene_col The column name or number of geneid in \code{exprMat}. Dafult is NULL, choose the first column as geneid.
#' @param from The type of matrix, support: "count", "tpm", "fpkm", "rpkm"
#' @param to Which type will be converted to, support: "tpm", "fpkm", "rpkm"
#' @param org Organism type, support human and mouse. We used the gen length information get from ensembl website (www.ensembl.org),
#'    if you use the custom gene length information, you can ignore it.
#' @param effLength The gene length information which contains gene id `id` and gene length `length`. Dafult is NULL, we use the
#'    information stored in this package.
#' @export
#'
convert_expr <- function(exprMat, idType = c("ensembl", "entrez", "symbol"), gene_col = NULL, from = "count", to = "tpm",
                        org = c("human", "mouse"), effLength = NULL) {

  ## get the Organism
  org <- match.arg(org)
  idType <- match.arg(idType)

  ## get gene expression
  if (is.null(gene_col)) {
    message("`gene_col = NULL`, choose first column as gene id.")
    gene_col <- 1
    colnames(exprMat)[1] <- "genes_id"
  } else if (is.character(gene_col)) {
    colnames(exprMat)[which(colnames(exprMat) == gene_col)] <- "genes_id"
  } else if (is.numeric(gene_col)) {
    colnames(exprMat)[gene_col] <- "genes_id"
  }
  gene_id <- exprMat %>% pull(genes_id)

  ## check whther have gene length information and format the ref data
  if (is.null(effLength)) {
    if (org == "human") {
      ref <- expareR::human
    } else if (org == "mouse") {
      ref <- expareR::mouse
    } else {
      stop("We only support `human` and `mouse`, if other species, please use `effLength` argument")
    }
    ref <- ref[, c(grep(idType, colnames(ref), value = T), "length")]
  } else {
    ref <- effLength
  }

  ## format the data
  colnames(ref) <- c('id', 'length')

  ## get order and match
  inter_gene <- intersect(ref$id, gene_id)
  if (length(inter_gene) != length(gene_id)) {
    cli::cli_alert_warning(col_red(paste(length(gene_id)-length(inter_gene), "gene can't be matched, filtering it ...")))
  }

  exprMat <- exprMat %>% filter(is.element(genes_id, inter_gene))
  ref <- ref[match(exprMat$genes_id, ref$id),]
  expr <- exprMat %>% select(!genes_id) %>% as.matrix()

  ## do the convert
  message(paste("Converting", from, "to", to, "..."))
  if (from == "count") {
    if (to == "tpm") {
      res <- apply(expr, 2, function(x) counts_to_tpm(x, ref$length))
    } else if (is.element(to, c("fpkm", "rpkm"))) {
      res <- apply(expr, 2, function(x) counts_to_fpkm(x, ref$length))
    } else {
      stop("We can't do this convert!")
    }
  } else if (is.element(from, c("fpkm", "rpkm"))) {
    if (to == "tpm") {
      res <- apply(expr, 2, function(x) fpkm_to_tpm(x))
    } else {
      stop("We can't do this convert!")
    }
  } else {
    stop("We can't do this convert!")
  }

  ## keep result
  res <- cbind(exprMat$genes_id, as.data.frame(res))
  colnames(res)[1] <- "genes"
  return(res)
}



##--------------- method to use -----------------------##

#---convert count to TPM---
counts_to_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

#---convert Fpkm/Rpkm to TPM---
fpkm_to_tpm <- function(fpkm) {
  #TPM = ( FPKM / sum of FPKM over all genes/transcripts) * 10^6
  #TPM = ( RPKM / sum of RPKM over all genes/transcripts) * 10^6
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#---convert count to FPKM---
counts_to_fpkm <- function(counts, lengths) {
  exp(log(counts) + log(1e9) - log(lengths) - log(sum(counts)) )
}

