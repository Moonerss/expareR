## code to prepare `gene_length` dataset goes here

#' prepare the gene length information for count convert
prepare_gene_len <- function(org = c("human", "mouse")) {
  org <- match.arg(org)
  datasets <- ifelse(org == "human", "hsapiens_gene_ensembl", "mmusculus_gene_ensembl")
  type <- c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "start_position", "end_position")
  if(org == "mouse") type[3] <- "mgi_symbol"


  mart <- biomaRt::useMart(host = "www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = datasets)
  ensembl <- biomaRt::getBM(attributes=type, mart = mart)
  ensembl$length <- abs(ensembl$end_position - ensembl$start_position)

  return(ensembl)
}

human <- prepare_gene_len(org = "human")
mouse <- prepare_gene_len(org = "mouse")
usethis::use_data(human, overwrite = TRUE)
usethis::use_data(mouse, overwrite = TRUE)
