#' Top-ranked genes that show differential variability
#'
#' @param clrDV_result An output from clrDV() function.
#' @param top large/small/abs. Large for SD ratio > 1; small for SD ratio < 1;
#'            abs. for ranking the DV genes by the absolute value of log2 SD ratio.
#' @param n The threshold rank.
#' @param full.table Default is FALSE. If full.table = TRUE, the maximum likelihood
#' estimates of parameters of skew-normal distribution for group 1 and group 2 will be included.
#'
#'
#' @description Extract top-ranked differential variability (DV) genes. These genes are ranked
#'              using the SD ratio (treatment vs. control) of the CLR-transformed count.
#'
#' @return
#' \item{DV}{The difference of the estimated standard deviation between group 2 and group 1 (sigma2 - sigma1).}
#' \item{se}{The standard error of DV.}
#' \item{z-value}{The observed Wald statistic value.}
#' \item{p-value}{The unadjusted p-value of the Wald test.}
#' \item{adjusted p-value}{The p-value of the Wald test adjusted using the Benjamini-Yekutieli procedure.}
#' \item{sd_ratio}{The ratio of sigma2 and sigma1 (sigma2/sigma1).}
#' \item{LFC}{log2 of the sd_ratio.}
#'
#' @export
#'
#' @import utils
#'
#' @examples
#'    library(clrDV)
#'    data("clrCounts2")
#'    group0 <- c(rep(0, 200), rep(1, 200))
#'    simu_clrDV_test <- clrDV(clrCounts2, group0)
#'    top.DV.genes(simu_clrDV_test, top = "large", n = 5)
#'    top.DV.genes(simu_clrDV_test, top = "small", n = 5)
#'    top.DV.genes(simu_clrDV_test, top = "abs.", n = 5)
#'    top.DV.genes(simu_clrDV_test, top = "abs.", n = 5, full.table = TRUE)
#'    dim(top.DV.genes(simu_clrDV_test, top = "abs.", n = Inf))
#'
#'
top.DV.genes <- function(clrDV_result, top = "abs.", n = 10, full.table = F){
  # top = "abs.", top large absolute value of log2 SD ratio
  # top = "large", top large SD ratio
  # top = "small", top small SD ratio
  # n <= number of DV genes; the default is top 10 DV genes
  dv_table <- clrDV_result[clrDV_result[, 5] < 0.05, ]
  SD_ratio <- dv_table[, 10]/dv_table[, 6] # sigma2/sigma1
  LFC <- log2(SD_ratio)
  dv_table <- cbind(dv_table, SD_ratio, LFC)


  if (n > dim(dv_table)[1] | n == Inf) {n <- dim(dv_table)[1]}

  if (top == "small"){
    top.dv.genes <- head(dv_table[order(dv_table[, 14]), ], n = n)
  } else if (top == "large"){
    top.dv.genes <- head(dv_table[order(dv_table[, 14], decreasing = T), ], n = n)
  } else if (top == "abs.") {
    top.dv.genes <- head(dv_table[order(abs(log2(dv_table[, 14])), decreasing = T), ], n = n)
  }

  if (full.table == T){
    return(top.dv.genes)
  } else {
    return(top.dv.genes[, c(1:5,14,15)])
  }
}
