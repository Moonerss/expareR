
<!-- README.md is generated from README.Rmd. Please edit that file -->

# expareR

<!-- badges: start -->
<!-- badges: end -->

The goal of expareR is to quickly prepare the expression matrix for
transcriptomics analysis.

## Installation

You can install the released version of `expareR` from
[github](https://github.com/Moonerss/expareR) with:

``` r
library(remotes)
remotes::install_github("Moonerss/expareR")
```

## Example

### For the RNA-seq data

``` r
library(expareR)
library(readr)
library(magrittr)

## read count matrix
count_expr <- read_delim(system.file("extdata", "count_expr.txt", package = "expareR"), delim = "\t")
#> 
#> -- Column specification --------------------------------------------------------
#> cols(
#>   .default = col_double(),
#>   symbol = col_character()
#> )
#> i<U+00A0>Use `spec()` for the full column specifications.
head(count_expr[,1:5])
#> # A tibble: 6 x 5
#>   symbol       GSM1523727 GSM1523728 GSM1523729 GSM1523744
#>   <chr>             <dbl>      <dbl>      <dbl>      <dbl>
#> 1 SH3KBP1            14.4       14.3       14.5       14.3
#> 2 RPL41              14.1       14.1       14.2       14.3
#> 3 CD24               14.3       14.3       14.2       14.2
#> 4 COX2               14.1       14.2       14.2       13.9
#> 5 LOC101928826       14.0       14.0       14.0       13.9
#> 6 HUWE1              14.0       14.0       13.9       14.0
```

first, you can convert the `count` to `tpm`:

``` r
# count to tpm
tpm_mat <- convert_expr(exprMat = count_expr, gene_col = "symbol", idType = "symbol", from = "count", to = "tpm", org = "human", effLength = NULL)
#> !<U+00A0>11 gene can't be matched, filtering it ...
#> Converting count to tpm ...
```

second, check the data quality after convert

``` r
## add bad data
tpm_mat[2,2] = NA
tpm_mat[3,3] = Inf

## check the data quality of count matrix
check_expr(tpm_mat, verbose = F)
#> `id_col = NULL`, choose first column as id.
#> !<U+00A0>Have some problem:
#> 1. Have `NA` in matrix
#> 2. Have `Inf` or `-Inf` in matrix
#> 3. Have duplicated genes
#> [1] FALSE
```

if the expression matrix have `NA` or `Inf` value, you can change them:

``` r
transed_mat <- transfer_data(tpm_mat, gene_col = "genes", data_type = c("all", "NA", "Inf"), data_to = 0)
```

and then, you can remove dupliacted genes by the expression, or combine
the gene expression

``` r
# remove duplicated genes by choose the one with bigger mean value
filter_dat <- remove_duplicate(transed_mat, gene_col = "genes", method = "order", value = "mean")

# remove duplicated genes by using the mean value of multiple gene instead of raw value
filter_dat <- remove_duplicate(transed_mat, gene_col = "genes", method = "combin", value = "mean")
```

finally, you log2 transformed the tpm matrix

``` r
# log2 transformed
log2_mat <- log2expr(filter_dat)
#> `gene_col = NULL`, choose first column as gene id.
#> i<U+00A0>log2 transformation finished!
head(log2_mat[,1:4])
#>          GSM1523727 GSM1523728 GSM1523729 GSM1523744
#> SH3KBP1    4.040651   4.037018   4.046019   4.042879
#> RPL41      0.000000  11.930430  11.932492  11.955328
#> CD24       9.854272   8.752582   9.837680   9.840981
#> HUWE1      5.147661   5.145952   5.137949   5.158486
#> SNORD42A  16.351119  16.347461  16.352296  16.362869
#> RPL37A     6.058324   6.062137   6.063762   6.079527
```

instead, you can use pipe connect the command

``` r
log2_mat <- count_expr %>% 
  convert_expr(gene_col = "symbol", idType = "symbol", from = "count", to = "tpm", org = "human", effLength = NULL) %>% 
  remove_duplicate(method = "order", value = "mean") %>% 
  transfer_data(data_type = "all", data_to = 0) %>% 
  log2expr()
#> `gene_col = NULL`, choose first column as gene id.
#> `gene_col = NULL`, choose first column as gene id.
#> !<U+00A0>11 gene can't be matched, filtering it ...
#> Converting count to tpm ...
#> `gene_col = NULL`, choose first column as gene id.
#> i<U+00A0>log2 transformation finished!
```

``` r
head(log2_mat[,1:4])
#>          GSM1523727 GSM1523728 GSM1523729 GSM1523744
#> SH3KBP1    4.040651   4.037018   4.046019   4.042879
#> RPL41     11.929083  11.930430  11.932492  11.955328
#> CD24       9.871570   9.874269   9.860082   9.877379
#> HUWE1      5.147661   5.145952   5.137949   5.158486
#> SNORD42A  16.351119  16.347461  16.352296  16.362869
#> RPL37A     6.058324   6.062137   6.063762   6.079527
```

### For the array data

If you need to process array data, you can use `anno_expr` to transform
probe to genes

``` r
## probe expression
array_expr <- read_delim(system.file("extdata", "array_expr.txt", package = "expareR"), delim = "\t")
#> 
#> -- Column specification --------------------------------------------------------
#> cols(
#>   .default = col_double(),
#>   probe = col_character()
#> )
#> i<U+00A0>Use `spec()` for the full column specifications.

## load annotation file 
data("anno_hug133plus2")

annotated_expr <- anno_expr(arrayMat = array_expr, anno = anno_hug133plus2, probe_col = "probe", symbol = "symbol", probe = "probe_id")
#> i<U+00A0>>>> 100.00% of probe in expression set was annotated
```

and you can process the data like RNA-seq analysis

``` r
log2_array <- array_expr %>% 
  anno_expr(anno = anno_hug133plus2, probe_col = "probe", symbol = "symbol", probe = "probe_id") %>% 
  remove_duplicate(method = "order", value = "mean") %>% 
  transfer_data(data_type = "all", data_to = 0) %>% 
  log2expr()
#> `gene_col = NULL`, choose first column as gene id.
#> `gene_col = NULL`, choose first column as gene id.
#> i<U+00A0>>>> 100.00% of probe in expression set was annotated
#> `gene_col = NULL`, choose first column as gene id.
#> i<U+00A0>log2 transformation is not necessary.

head(log2_array[,1:4])
#>         GSM2696792 GSM2696793 GSM2696794 GSM2696795
#> MIR4640  10.261650  10.740100  10.694980  10.692940
#> RFC2      8.365080   8.463314   8.535923   8.091345
#> HSPA6     7.212391   5.832635   5.616512   8.090054
#> PAX8      7.377803   7.549826   8.186557   7.459644
#> GUCA1A    6.699024   2.621356   2.422049   2.423371
#> MIR5193   7.487896   7.255971   6.111238   6.941331
```
