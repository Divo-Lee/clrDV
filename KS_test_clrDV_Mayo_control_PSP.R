###############################################################
#clrDV: A differential variability test for
#RNA-Seq data based on the skew-normal distribution
#Authors: Hongxiang Li and Tsung Fei Khang
#Email: chelsea.divo@hotmail.com
#Latest update: 10 April 2023
#R Codes for DV test on Mayo RNA-Seq dataset
#Part 5: Goodness-of-fit Assessment of the skew-normal model on 
#        the distribution of the CLR-transformed Mayo RNA-Seq
#        data (control vs. PSP) 
###############################################################

## R packages downloaded
# install.packages("BiocManager")
# install.packages("devtools")
# install.packages("compositions")
# BiocManager::install("edgeR", force = T)
# devtools::install_github("Divo-Lee/clrDV")
# install.packages("sn")
library(compositions); library(edgeR)
library(clrDV); library(sn)


###################################
### KS test for control vs. PSP ###
###################################

## Read data
# The use of the data files MayoRNAseq_individual_metadata_031422.csv
# and MayoRNAseq_RNAseq_TCX_geneCounts.csv requires permission from the data owners.
# Request permission from https://adknowledgeportal.synapse.org/ (look for Mayo RNAseq study)

## raw count table
ad_counts <- read.csv('MayoRNAseq_RNAseq_TCX_geneCounts.csv', row.names = 1)
dim(ad_counts)
## meta-data
ad_meta <- read.csv('MayoRNAseq_individual_metadata_031422.csv')
# remove samples give NA in disease name column
ad_meta <- ad_meta[-which(is.na(ad_meta$diagnosis)), ]


# samples id in meta_data
control_id <- (ad_meta$individualID)[ad_meta$diagnosis == "control"]
psp_id <-(ad_meta$individualID)[ad_meta$diagnosis == "progressive supranuclear palsy"]
# samples id in raw count table
ad_counts_id <-   colnames(ad_counts)
#
control_vector <- sapply(ad_counts_id, function(k) k %in% control_id)
control_counts <-  ad_counts[, control_vector]
psp_vector <- sapply(ad_counts_id, function(k) k %in% psp_id)
psp_counts <- ad_counts[, psp_vector]
N_psp <- length(control_counts) + length(psp_counts)
mayo_counts2 <- as.matrix(cbind(control_counts, psp_counts),
                          nrows=dim(ad_counts)[1], ncol = N_psp)

dim(control_counts); dim(psp_counts); dim(mayo_counts2)

# Filtering
CPM_psp <- cpm(mayo_counts2)
keep_psp <- rowMeans(CPM_psp[,1:length(control_counts)]) > 0.5 &
  rowMeans(CPM_psp[,(length(control_counts)+1):N_psp]) > 0.5 &
  apply(mayo_counts2[,1:length(control_counts)], 1, function(k) length(k[k == 0])/length(k)) < 0.85 &
  apply(mayo_counts2[,(length(control_counts)+1):N_psp], 1, function(k) length(k[k == 0])/length(k)) < 0.85

mayo_counts_filter2 <- mayo_counts2[keep_psp, ]
dim(mayo_counts_filter2)


############
group3 = c(rep(0, length(control_counts)), rep(1, length(psp_counts)))

# CLR-transformation
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}

# MLE/MPLE
t2 <- proc.time()
clr.counts3 <- clr.transform(data = mayo_counts_filter2)
clrseq_PSP <- clrSeq(clr.counts3, group = group3)
as.numeric(proc.time() - t2)[3] # computing time, in seconds

control2_sn_CP <- clrseq_PSP[, c("mu1", "sigma1", "gamma1")]
PSP_sn_CP <- clrseq_PSP[, c("mu2", "sigma2", "gamma2")]

## map centered parameters (CP) to direct parameters (DP)
cp_to_dp <- function(mean=NULL, sd=NULL, skewness=NULL){
  b <- sqrt(2/pi)
  
  if(skewness >= 0){
    r <- (2*skewness/(4-pi))^(1/3)
  } else {
    r <- -(2*(- skewness)/(4-pi))^(1/3)
  }
  
  alpha <- r/sqrt(2/pi - (1-2/pi)*r^2)
  delta <- alpha/sqrt(1 + alpha^2)
  omega <- sd/sqrt(1 - (b^2)*delta^2)
  xi <- mean - b*omega*delta
  return(c(xi, omega, alpha))
}


###
control2_sn_DP <- matrix(NA, nrow = dim(control2_sn_CP)[1], ncol = 3)
PSP_sn_DP <- matrix(NA, nrow = dim(PSP_sn_CP)[1], ncol = 3)

for (i in 1:dim(control2_sn_CP)[1]) {
  control2_sn_DP[i,] <- c(cp_to_dp(control2_sn_CP[i,1],
                                   control2_sn_CP[i,2],
                                   control2_sn_CP[i,3]))
}

colnames(control2_sn_DP) <- c("xi1", "omega1", "alpha1")
control2_sn_DP <- as.data.frame(control2_sn_DP)


for (i in 1:dim(PSP_sn_CP)[1]) {
  PSP_sn_DP[i,] <- c(cp_to_dp(PSP_sn_CP[i,1],
                              PSP_sn_CP[i,2],
                              PSP_sn_CP[i,3]))
}

colnames(PSP_sn_DP) <- c("xi2", "omega2", "alpha2")
PSP_sn_DP <- as.data.frame(PSP_sn_DP)


# KS test for the control group (control vs. PSP)
control2_KS_test_pvalue <- vector()
for (i in 1:dim(control2_sn_CP)[1]) {
  ks <- ks.test(clr.counts3[i, 1:78],
                "psn",
                xi = control2_sn_DP$xi1[i],
                omega = control2_sn_DP$omega1[i],
                alpha = control2_sn_DP$alpha1[i])
  control2_KS_test_pvalue[i] <- ks$p.value
}
sum(control2_KS_test_pvalue < 0.05)
# 281/18636 = 0.015
# 98.5% of genes in the control group fitted skew-normal distribution well

# Fig. 4 (c)
hist(sqrt(-log10(control2_KS_test_pvalue)),
     freq = F, breaks = 25,
     xlab =expression(sqrt(-log[10](p))),
     xlim = c(0,4), ylim = c(0,2.6),
     main = NULL, border = "grey83")
abline(v = sqrt(-log10(0.05)), lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(c)")


## KS test for the PSP group (control vs. PSP)
PSP_KS_test_pvalue <- vector()
for (i in 1:dim(control2_sn_CP)[1]) {
  ks <- ks.test(clr.counts3[i, 79:162],
                "psn",
                xi = PSP_sn_DP$xi2[i],
                omega = PSP_sn_DP$omega2[i],
                alpha = PSP_sn_DP$alpha2[i])
  PSP_KS_test_pvalue[i] <- ks$p.value
}
sum(PSP_KS_test_pvalue < 0.05)
# 147/18636 = 0.008
# 99.2% of genes in the AD group fitted skew-normal distribution well

# Fig. 4 (d)
hist(sqrt(-log10(PSP_KS_test_pvalue)), freq = F, breaks = 50,
     main = NULL, border = "grey83",
     xlim = c(0, 4), ylim = c(0, 2.6),
     xlab = expression(sqrt(-log[10](p))))
abline(v = sqrt(-log10(0.05)), lty = 2, lwd = 1.25, col = "red")
title(adj=0, "(d)")


# union of the genes, fit skew-normal well
ks_pvalue_mat <- cbind.data.frame(control2_KS_test_pvalue,
                                  PSP_KS_test_pvalue)
colnames(ks_pvalue_mat) <- c("control", "PSP")
dim(ks_pvalue_mat[(ks_pvalue_mat$control < 0.05 | ks_pvalue_mat$PSP < 0.05), ])
# 375/18636 = 0.02
# overall, the skew-normal distribution fits 98% of the genes well in both the control and PSP groups

###END###
