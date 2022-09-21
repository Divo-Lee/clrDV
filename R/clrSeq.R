#' Fitting the skew-normal distribution to CLR-transformed RNA-Seq data for 2 groups
#'
#' @param data A table of clr-transformed count data, which genes/transcripts on the rows and
#'  samples on columns for 2 groups.
#' @param group A vector specifying the group labels of the data.
#'
#' @return
#'  \item{mu1}{The maximum likelihood estimate of mean parameter for group 1.}
#'  \item{se.mu1}{The standard error of the maximum likelihood estimate of mean parameter for group 1.}
#'  \item{z.mu1}{The Wald statistic for mu1.}
#'  \item{p.mu1}{The p-value of the Wald test for mu1.}
#'  \item{sigma1}{The maximum likelihood estimate of the standard deviation parameter for group 1.}
#'  \item{se.sigma1}{The standard error of the maximum likelihood estimate of the standard deviation parameter for group 1.}
#'  \item{z.sigma1}{The Wald statistic for sigma1.}
#'  \item{p.sigma1}{The p-value of the Wald test for sigma1.}
#'  \item{gamma1}{The maximum likelihood estimate of the skewness parameter for group 1.}
#'  \item{se.gamma1}{The standard error of the maximum likelihood estimate of the skewness parameter for group 1.}
#'  \item{z.gamma1}{The Wald statistic for gamma1.}
#'  \item{p.gamma1}{The p-value of the Wald test for gamma1.}
#'  \item{mu2}{The maximum likelihood estimate of mean parameter for group 2.}
#'  \item{se.mu2}{The standard error of the maximum likelihood estimate of mean parameter for group 2.}
#'  \item{z.mu2}{The Wald statistic for mu2.}
#'  \item{p.mu2}{The p-value of the Wald test for mu2.}
#'  \item{sigma2}{The maximum likelihood estimate of the standard deviation parameter for group 2.}
#'  \item{se.sigma2}{The standard error of the maximum likelihood estimate of the standard deviation parameter for group 2.}
#'  \item{z.sigma2}{The Wald statistic for sigma2.}
#'  \item{p.sigma2}{The p-value of the Wald test for sigma2.}
#'  \item{gamma2}{The maximum likelihood estimate of the skewness parameter for group 2.}
#'  \item{se.gamma2}{The standard error of the maximum likelihood estimate of the skewness parameter for group 2.}
#'  \item{z.gamma2}{The Wald statistic for gamma2.}
#'  \item{p.gamma2}{The p-value of the Wald test for gamma2.}
#'
#' @export
#'
#' @examples
#'    library(clrDV)
#'    data("clrCounts2")
#'    group0 <- c(rep(0, 200), rep(1, 200))
#'    clrSeq(clrCounts2[c(1:5),], group0)
#'
#'    clr_Seq <- clrSeq(clrCounts2, group0)
#'    tail(clr_Seq, 5)
#'
#'
#'
 clrSeq <- function(data = NULL, group = NULL){
  # only for two group
  # data must be clr-transformed counts
  if (is.factor(group) == F){group = as.factor(group)}
  colnames(data) <- as.factor(group)

  clr_count_group1 <- data[, group == levels(group)[1]]
  clr_count_group2 <- data[, group == levels(group)[2]]

  group1_seq <- clr.SN.fit(clr_count_group1)
  group2_seq <- clr.SN.fit(clr_count_group2)

  result <- cbind(group1_seq, group2_seq)
  row.names(result) <- row.names(data)
  colnames(result) <-  c("mu1", "se.mu1", "z-value", "p-value",
                         "sigma1", "se.sigma1", "z-value", "p-value",
                         "gamma1", "se.gamma1", "z-value", "p-value",
                         "mu2", "se.mu2", "z-value", "p-value",
                         "sigma2", "se.sigma2", "z-value", "p-value",
                         "gamma2", "se.gamma2", "z-value", "p-value")
  return(result)
}
