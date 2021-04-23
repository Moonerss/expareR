#' @title Convert read counts to transcripts per million (TPM)
#' @description Convert the count matrix to TPM matrix
#' @name convert_expr
#' @param exprMat A expression matrix, with geneid as rownames and sample as columns.
#' @param idType Type of gene id. "ensembl" "entrez","symbol"
#' @param from The type of matrix, support: "count", "tpm", "fpkm", "rpkm"
#' @param to Which type will be converted to, support: "tpm", "fpkm", "rpkm"
#' @param org Organism type, support human and mouse. We used the gen length information get from ensembl website (www.ensembl.org),
#'    if you use the custom gene length information, you can ignore it.
#' @param effLength The gene length information which contains gene id `id` and gene length `length`. Dafult is NULL, we use the
#'    information stored in this package.
#' @export
#'
#' @examples
#'
convet_expr <- function(exprMat, idType = c("ensembl", "entrez", "symbol"), from = "count", to = "tpm",
                        org = c("human", "mouse"), effLength = NULL) {

  ## check the data format
  if (class(exprMat)!="matrix") exprMat<-as.matrix(exprMat)

  ## get the Organism
  org <- match.arg(org)
  idType <- match.arg(idType)

  ## check whther have gene length information and format the ref data
  if (is.null(effLength)) {
    ref <- ifelse(org == "human", expareR::human, expareR::mouse)
    ref <- ref[, c(grep(idType, colnames(ref), value = T), "length")]
  } else {
    ref <- effLength
  }

  ## format the data
  colnames(ref) <- c('id', 'length')

  ## get order and match
  inter_gene <- intersect(ref$id, rownames(exprMat))
  if (length(inter_gene) != nrow(exprMat)) {
    cli::cli_alert_warning(col_red(paste(nrow(exprMat)-length(inter_gene), "gene can't be matched, filtering it ...")))
  }
  exprMat <- exprMat[inter_gene,]
  ref <- ref[match(inter_gene, ref$id),]

  ## do the convert
  message(paste("Converting", from, "to", to, "..."))
  if (from == "count") {
    if (to == "tpm") {
      res <- apply(exprMat, 2, function(x) counts_to_tpm(x, ref$length))
    } else if (is.element(to, c("fpkm", "rpkm"))) {
      res <- apply(exprMat, 2, function(x) counts_to_fpkm(x, ref$length))
    } else {
      stop("We can't do this convert!")
    }
  } else if (is.element(from, c("fpkm", "rpkm"))) {
    if (to == "tpm") {
      res <- apply(exprMat, 2, function(x) fpkm_to_tpm(x))
    } else {
      stop("We can't do this convert!")
    }
  } else {
    stop("We can't do this convert!")
  }

  ## keep result
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
count_to_fpkm <- function(counts, lengths) {
  exp(log(counts) + log(1e9) - log(lengths) - log(sum(counts)) )
}

