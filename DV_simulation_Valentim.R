#####################################################
#clrDV: A differential variability test for
#RNA-Seq data based on the skew-normal distribution
#Author: Hongxiang Li
#Email: chelsea.divo@hotmail.com
#Date: 25 September 2022
#R Codes for comparison of methods
#Part 2: simulation study using the Valentim dataset 
######################################################

## Required R packages
# install.packages("compositions")
# install.packages("BiocManager")
# install.packages("devtools")
# BiocManager::install("polyester")
# BiocManager::install("edgeR", force = T)
# install.packages("gamlss")
# devtools::install_github("Divo-Lee/clrDV")
# BiocManager::install("GenomicAlignments")
# BiocManager::install("bumphunter", force = T)
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
# BiocManager::install("missMethyl", force = T)
# BiocManager::install("cqn")
# devtools::install_github("zjdaye/MDSeq")


library(polyester); library(compositions); library(edgeR)
library(MDSeq); library(missMethyl); library(gamlss)

###############################################
### GSE123658_read_counts.gene_level.txt.gz 
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123658
### Valentim Dataset                        
###############################################
Data = read.table(gzfile('GSE123658_read_counts.gene_level.txt.gz'), sep="\t", header = T, check.names = TRUE)
row.names(Data) <- Data[, 1]; Data <- Data[, -1]; dim(Data)
 ## Filtering
CPM <- cpm(Data)
keep <- (rowMeans(CPM[,1:43]) > 0.5  & rowMeans(CPM[,44:82]) > 0.5 & apply(Data[,1:43], 1, function(x) length(x[x==0])/length(x)) < 0.85 & apply(Data[,44:82], 1, function(x) length(x[x==0])/length(x)) < 0.85)
Data <- Data[keep, ]; dim(Data)
  # data_h <- Data[,1:43] # h = healthy group
data_d1 <- Data[,44:82] # d1 = type 1 diabetic group


### Functions
# Modified create_read_numbers() from the polyester R package
 ## simulate DV count table
create_dv_counts <- function(mu, fit, p0 = NULL,
                             N.genes = NULL,
                             N.samples = NULL, # each group
                             size.fc = NULL,
                             seed = NULL){

  if(!is.null(seed)){set.seed(seed)}
  if(is.null(N.genes) | is.null(N.samples)){
    stop(.makepretty("Please specify N.genes and N.spamples.\n"))
  }

  index = sample(1:length(mu), size = N.genes)

  means <- mu[index]
  p0s <- p0[index]
  mean.Mat <- log(means + 0.001) %*% t(rep(1, N.samples))
  mean.vec <- as.vector(mean.Mat)
  size.vec <- predict(fit, mean.vec)$y
  size.Mat <- matrix(size.vec, nrow = N.genes)
  counts.Mat0 <- size.Mat*NA
  counts.Mat1 <- size.Mat*NA

  if (is.null(size.fc)){size.fc <- rep(1, N.genes)}

  if(is.null(p0)){
    for (i in 1:N.genes) {
      counts.Mat0[i, ] <- rnbinom(N.samples,
                                  mu = exp(mean.Mat[i, ]),
                                  size = exp(size.Mat[i, ]))
    }
  } else {
    for (i in 1:N.genes) {
      counts.Mat0[i, ] <- rbinom(N.samples, prob = (1 - p0s[i]), size = 1)*
        rnbinom(N.samples, mu = exp(mean.Mat[i, ]), size = exp(size.Mat[i, ]))
    }
  }

  if(is.null(p0)){
    for (i in 1:N.genes) {
      counts.Mat1[i, ] <- rnbinom(N.samples,
                                  mu = exp(mean.Mat[i, ]),
                                  size = size.fc[i]*exp(size.Mat[i, ]))
    }
  } else {
    for (i in 1:N.genes) {
      counts.Mat1[i, ] <- rbinom(N.samples, prob = (1 - p0s[i]), size = 1)*
        rnbinom(N.samples, mu = exp(mean.Mat[i, ]), size = size.fc[i]*exp(size.Mat[i, ]))
    }
  }
  return(cbind(counts.Mat0, counts.Mat1))
}


 ## CLR-transformation
clr.transform <- function(data = NULL){
  data[data == 0] <- 1/2
  clr.count <- t(clr(t(data)))
  clr.count <- matrix(as.numeric(clr.count),
                      nrow = dim(data)[1],
                      ncol = dim(data)[2])
  row.names(clr.count) <- row.names(data)
  return(clr.count)
}


 ## Function for Comparisons of Methods
DV.test.comparison <- function(data = NULL,
                               N.genes = 2000,
                               N.samples = NULL,
                               prob.dv = 0.1,
                               N.simulations = 30,
                               seed = NULL){
  set.seed(seed)
  seeds <- sample(1:100000, size = 2*N.simulations)
  seeds1 <- seeds[1:N.simulations]
  seeds2 <- seeds[(1 + N.simulations):(2*N.simulations)]

  FC <- c(seq(0.25, 1/2, 0.0001), seq(2, 4, 0.0008))
  N.G <- N.genes
  N.dv <- round(prob.dv*N.G)
  params <- get_params(data)
  n = N.samples # Number of samples for each group

  names <- c("clrDV", "MDSeq", "diffVar", "GAMLSS_BH", "GAMLSS_BY")
  FDR <- matrix(NA, nrow = length(names), ncol = N.simulations)
  Type_II_error <- matrix(NA, nrow = length(names), ncol = N.simulations)
  Time <- matrix(NA, nrow = length(names), ncol = N.simulations)
  row.names(FDR) <- names
  row.names(Type_II_error)<- names
  row.names(Time)<- names

  ## Methods comparison
  i <- 1
  while (i <= N.simulations) {
    set.seed(seeds1[i])
    d.size<- sample(FC, size = N.dv)
    d <- N.G
    s.fc <- c(d.size, rep(1, d-N.dv))

    dv_data <- create_dv_counts(params$mu, params$fit,
                                N.genes = N.G,
                                N.samples = n,
                                size.fc = s.fc,
                                seed = seeds2[i])
    row.names(dv_data) <- paste('gene', 1:N.genes, sep='')
    dv_gene <- row.names(dv_data)[1:N.dv]
    data2 <- dv_data

    ### Data Filter
    CPM2 <- cpm(data2)
    keep2 <- (rowMeans(CPM2[,1:n]) > 0.5 & rowMeans(CPM2[,(n+1):(2*n)]) > 0.5 & apply(data2[,1:n], 1, function(x) length(x[x == 0])/length(x)) < 0.85   & apply(data2[,(n+1):(2*n)], 1, function(x) length(x[x == 0])/length(x)) < 0.85 )
    # sum(keep2); sum(keep2[1:N.dv])
    N.genes_dv_filter <- sum(keep2)
    N.dv_filter <- sum(keep2[1:N.dv])
    dv_gene <- dv_gene[keep2[1:N.dv]]

    data2 <- data2[keep2, ]


    ######################
    ### DV Test Comparison
    group2 = rep(c(0,1), each = n)

    ### clrDV
    t1 <- proc.time()
    clr.counts2 <- clr.transform(data = data2)
    s.d._test <- clrDV(clr.counts2, group = group2)
    s.d._test <- na.omit(s.d._test)

    dv_table <- s.d._test[s.d._test[, 5] < 0.05, ]
    dv_genes_clrDV <- row.names(dv_table) # DV genes called

    FDR_clrDV <- (sum(s.d._test[,5] < 0.05) -  length(intersect(dv_genes_clrDV, dv_gene)))/sum(s.d._test[, 5] < 0.05)
    type_II_error_clrDV <- (N.dv_filter - length(intersect(dv_genes_clrDV, dv_gene)))/N.dv_filter

    FDR[1, i] <- FDR_clrDV
    Type_II_error[1, i] <-type_II_error_clrDV
    Time[1, i] <- as.numeric(proc.time() - t1)[3]


    ### MDSeq
    t2 <- proc.time()
    libsizes2 <- colSums(data2)
    nf2 <- calcNormFactors(data2, method="TMM")
    els2 <- nf2 * libsizes2
    sf2 <- els2 / exp(mean(log(libsizes2)))
    norm.data2 <- t(t(data2) / sf2)

    contrasts2 <- get.model.matrix(as.factor(group2))
    fit.MDSeq.dv <- MDSeq(norm.data2,
                          contrast = contrasts2)
    res.MDSeq.dv <- extract.ZIMD(fit.MDSeq.dv,
                                 get='contrast',
                                 compare=list(A="0",B="1"),
                                 log2FC.threshold = 0)

    res.MDSeq.dv <- na.omit(res.MDSeq.dv)

    dv_table <- res.MDSeq.dv[res.MDSeq.dv$FDR.dispersion < 0.05, ]
    dv_genes_MDSeq <- row.names(dv_table) # DV genes called

    FDR_MDSeq <- (sum(res.MDSeq.dv$FDR.dispersion < 0.05) - length(intersect(dv_genes_MDSeq, dv_gene)))/sum(res.MDSeq.dv$FDR.dispersion < 0.05)
    type_II_error_var_MDSeq <- (N.dv_filter -  length(intersect(dv_genes_MDSeq, dv_gene)))/N.dv_filter

    FDR[2, i] <- FDR_MDSeq
    Type_II_error[2, i] <- type_II_error_var_MDSeq
    Time[2, i] <- as.numeric(proc.time() - t2)[3]


    ### diffVar
    t3 <- proc.time()
    design2 = model.matrix(~group2)

    libsizes2 <- colSums(data2)
    nf2 <- calcNormFactors(data2, method="TMM")
    els2 <- nf2 * libsizes2
    sf2 <- els2 / exp(mean(log(libsizes2)))
    norm.data2 <- t(t(data2) / sf2)


    fit.diffVar <- varFit(norm.data2, design2, coef=c(1,2))
    res.diffVar <- topVar(fit.diffVar, coef=2, number=nrow(norm.data2), sort=F)

    FDR_diffVar <- (sum(res.diffVar$Adj.P.Value < 0.05) - sum(res.diffVar$Adj.P.Value[1:N.dv_filter] < 0.05))/sum(res.diffVar$Adj.P.Value < 0.05)
    type_II_error_diffVar <- (N.dv_filter - sum(res.diffVar$Adj.P.Value[1:N.dv_filter] < 0.05))/N.dv_filter

    FDR[3, i] <- FDR_diffVar
    Type_II_error[3, i] <- type_II_error_diffVar
    Time[3, i] <- as.numeric(proc.time() - t3)[3]


    ### GAMLSS
    t4 <- proc.time()
    design2 = model.matrix(~group2)
    libsizes2 <- colSums(data2)
    nf2 <- calcNormFactors(data2, method="TMM")

    dat2.edgeR <- DGEList(counts=data2, norm.factors=nf2, group=group2)
    dat2.edgeR <- estimateDisp(dat2.edgeR, design2)

    dat2.edgeR$CPM <- cpm.DGEList(dat2.edgeR) # uses same data object as edgeR
    ofs <- log(dat2.edgeR$samples$lib.size * dat2.edgeR$samples$norm.factors)
    dat2.edgeR$samples$offset <- ofs
    gene_i <- seq_along(dat2.edgeR$counts[,1])

    gamlss_NB <- lapply(gene_i, function(i) {
      dat <- data.frame(x = dat2.edgeR$samples$group,
                        y = dat2.edgeR$counts[i,],
                        ofs = dat2.edgeR$samples$offset)
      dat$x <- relevel(dat$x, ref = c("1"))
      m0 <- tryCatch(gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
                            family = NBI(), sigma.start = 0.1, n.cyc = 100),
                     warning= function(w) NULL, error= function(e) NULL)
      m1 <- tryCatch(gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
                            family = NBI(), sigma.start = 0.1, n.cyc = 100),
                     warning= function(w) NULL, error= function(e) NULL)
      res <- data.frame(CV.1 = NA, CV.2 = NA, LR.cv = NA, p.cv = NA)

      if(!any(sapply(list(m0,m1), is.null)))
      {
        res$CV.1 = sqrt(exp(m0$sigma.coefficients)[[1]])
        res$CV.2 = sqrt(exp(m0$sigma.coefficients)[[2]])
        res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]]) /
                           sqrt(exp(m0$sigma.coefficients)[[1]]))
        res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
      }
      res
    })

    gamlss_NB <- do.call(rbind, gamlss_NB)
    rownames(gamlss_NB) <- rownames(dat2.edgeR$counts)[gene_i]
    gamlss_NB <- na.omit(gamlss_NB)
    t5 <- proc.time()

    ## BH FDR
    gamlss_NB$padj.cv <- p.adjust(gamlss_NB$p.cv, "fdr")

    dv_table <- gamlss_NB[gamlss_NB$padj.cv < 0.05, ]
    dv_genes_gamlss <- row.names(dv_table) # DV genes called

    FDR_gamlss <- (sum(gamlss_NB$padj.cv < 0.05) - length(intersect(dv_genes_gamlss, dv_gene)))/sum(gamlss_NB$padj.cv < 0.05)
    type_II_error_gamlss <- (N.dv_filter - length(intersect(dv_genes_gamlss, dv_gene)))/N.dv_filter

    FDR[4, i] <- FDR_gamlss
    Type_II_error[4, i] <- type_II_error_gamlss
    Time[4, i] <- as.numeric(proc.time() - t4)[3]

    ## BY FDR
    t6 <- proc.time()
    gamlss_NB$padj.cv <- p.adjust(gamlss_NB$p.cv, "BY")
    dv_table <- gamlss_NB[gamlss_NB$padj.cv < 0.05, ]
    dv_genes_gamlss <- row.names(dv_table)

    FDR_gamlss <- (sum(gamlss_NB$padj.cv < 0.05) - length(intersect(dv_genes_gamlss, dv_gene)))/sum(gamlss_NB$padj.cv < 0.05)
    type_II_error_gamlss <- (N.dv_filter - length(intersect(dv_genes_gamlss, dv_gene)))/N.dv_filter

    FDR[5, i] <- FDR_gamlss
    Type_II_error[5, i] <- type_II_error_gamlss
    Time[5, i] <- as.numeric(proc.time() - t6)[3] + as.numeric(t5 - t4)[3]

    i <- i + 1
  }
  return(list("FDR" =  FDR, "Type_II_error" = Type_II_error, "Time" = Time))
}



### Comparison of methods
### Be patient, this will take some time  
## Comparison One, 2*50 samples
DV_test1 <- DV.test.comparison(data = data_d1, N.genes = 2000,
                               N.samples = 50, prob.dv = 0.1,
                               N.simulations = 30, seed = 1)
DV_test1$FDR; DV_test1$Type_II_error; DV_test1$Time

write.csv(round(DV_test1$FDR, 4), "DV_FDR_diabetes1_50.csv")
write.csv(round(DV_test1$Type_II_error, 4), "DV_Type_II_error_diabetes1_50.csv")
write.csv(round(DV_test1$Time, 2), "DV_Time_diabetes1_50.csv")


## Comparison Two, 2*100 samples
DV_test2 <- DV.test.comparison(data = data_d1, N.genes = 2000,
                               N.samples = 100, prob.dv = 0.1,
                               N.simulations = 30, seed = 2)
DV_test2$FDR; DV_test2$Type_II_error; DV_test2$Time

write.csv(round(DV_test2$FDR, 4), "DV_FDR_diabetes1_100.csv")
write.csv(round(DV_test2$Type_II_error, 4), "DV_Type_II_error_diabetes1_100.csv")
write.csv(round(DV_test2$Time, 2), "DV_Time_diabetes1_100.csv")


## Comparison Three, 2*125 samples
DV_test3 <- DV.test.comparison(data = data_d1, N.genes = 2000,
                               N.samples = 125, prob.dv = 0.1,
                               N.simulations = 30, seed = 3)
DV_test3$FDR; DV_test3$Type_II_error; DV_test3$Time

write.csv(round(DV_test3$FDR, 4), "DV_FDR_diabetes1_125.csv")
write.csv(round(DV_test3$Type_II_error, 4), "DV_Type_II_error_diabetes1_125.csv")
write.csv(round(DV_test3$Time, 2), "DV_Time_diabetes1_125.csv")


## Comparison Four, 2*150 samples
DV_test4 <- DV.test.comparison(data = data_d1, N.genes = 2000,
                               N.samples = 150, prob.dv = 0.1,
                               N.simulations = 30, seed = 4)
DV_test4$FDR; DV_test4$Type_II_error; DV_test4$Time

write.csv(round(DV_test4$FDR, 4), "DV_FDR_diabetes1_150.csv")
write.csv(round(DV_test4$Type_II_error, 4) , "DV_Type_II_error_diabetes1_150.csv")
write.csv(round(DV_test4$Time, 2), "DV_Time_diabetes1_150.csv")


## Comparison Five, 2*200 samples
DV_test5 <- DV.test.comparison(data = data_d1, N.genes = 2000,
                               N.samples = 200, prob.dv = 0.1,
                               N.simulations = 30, seed = 5)
DV_test5$FDR; DV_test5$Type_II_error; DV_test5$Time

write.csv(round(DV_test5$FDR, 4), "DV_FDR_diabetes1_200.csv")
write.csv(round(DV_test5$Type_II_error, 4), "DV_Type_II_error_diabetes1_200.csv")
write.csv(round(DV_test5$Time, 2), "DV_Time_diabetes1_200.csv")


## Comparison Six, 2*250 samples
DV_test6 <- DV.test.comparison(data = data_d1, N.genes = 2000,
                               N.samples = 250, prob.dv = 0.1,
                               N.simulations = 30, seed = 6)
DV_test6$FDR; DV_test6$Type_II_error; DV_test6$Time

write.csv(round(DV_test6$FDR, 4), "DV_FDR_diabetes1_250.csv")
write.csv(round(DV_test6$Type_II_error, 4), "DV_Type_II_error_diabetes1_250.csv")
write.csv(round(DV_test6$Time, 2), "DV_Time_diabetes1_250.csv")



########################
### DV Test Analysis ###
########################
dv.FDR.d1_new.50samples <- read.csv('DV_FDR_diabetes1_50.csv')
dv.FDR.d1_new.100samples <- read.csv('DV_FDR_diabetes1_100.csv')
dv.FDR.d1_new.125samples <- read.csv('DV_FDR_diabetes1_125.csv')
dv.FDR.d1_new.150samples <- read.csv('DV_FDR_diabetes1_150.csv')
dv.FDR.d1_new.200samples <- read.csv('DV_FDR_diabetes1_200.csv')
dv.FDR.d1_new.250samples <- read.csv('DV_FDR_diabetes1_250.csv')

row.names(dv.FDR.d1_new.50samples) <- as.factor(dv.FDR.d1_new.50samples[, 1])
dv.FDR.d1_new.50samples <- dv.FDR.d1_new.50samples[, -1]
row.names(dv.FDR.d1_new.100samples) <- as.factor(dv.FDR.d1_new.100samples[, 1])
dv.FDR.d1_new.100samples <- dv.FDR.d1_new.100samples[, -1]
row.names(dv.FDR.d1_new.125samples) <- as.factor(dv.FDR.d1_new.125samples[, 1])
dv.FDR.d1_new.125samples <- dv.FDR.d1_new.125samples[, -1]
row.names(dv.FDR.d1_new.150samples) <- as.factor(dv.FDR.d1_new.150samples[, 1])
dv.FDR.d1_new.150samples <- dv.FDR.d1_new.150samples[, -1]
row.names(dv.FDR.d1_new.200samples) <- as.factor(dv.FDR.d1_new.200samples[, 1])
dv.FDR.d1_new.200samples <- dv.FDR.d1_new.200samples[, -1]
row.names(dv.FDR.d1_new.250samples) <- as.factor(dv.FDR.d1_new.250samples[, 1])
dv.FDR.d1_new.250samples <- dv.FDR.d1_new.250samples[, -1]


dv.type2.d1_new.50samples <- read.csv('DV_Type_II_error_diabetes1_50.csv')
dv.type2.d1_new.100samples <- read.csv('DV_Type_II_error_diabetes1_100.csv')
dv.type2.d1_new.125samples <- read.csv('DV_Type_II_error_diabetes1_125.csv')
dv.type2.d1_new.150samples <- read.csv('DV_Type_II_error_diabetes1_150.csv')
dv.type2.d1_new.200samples <- read.csv('DV_Type_II_error_diabetes1_200.csv')
dv.type2.d1_new.250samples <- read.csv('DV_Type_II_error_diabetes1_250.csv')

row.names(dv.type2.d1_new.50samples) <- as.factor(dv.type2.d1_new.50samples[, 1])
dv.type2.d1_new.50samples <- dv.type2.d1_new.50samples[, -1]
row.names(dv.type2.d1_new.100samples) <- as.factor(dv.type2.d1_new.100samples[, 1])
dv.type2.d1_new.100samples <- dv.type2.d1_new.100samples[, -1]
row.names(dv.type2.d1_new.125samples) <- as.factor(dv.type2.d1_new.125samples[, 1])
dv.type2.d1_new.125samples <- dv.type2.d1_new.125samples[, -1]
row.names(dv.type2.d1_new.150samples) <- as.factor(dv.type2.d1_new.150samples[, 1])
dv.type2.d1_new.150samples <- dv.type2.d1_new.150samples[, -1]
row.names(dv.type2.d1_new.200samples) <- as.factor(dv.type2.d1_new.200samples[, 1])
dv.type2.d1_new.200samples <- dv.type2.d1_new.200samples[, -1]
row.names(dv.type2.d1_new.250samples) <- as.factor(dv.type2.d1_new.250samples[, 1])
dv.type2.d1_new.250samples <- dv.type2.d1_new.250samples[, -1]



# FDR and average type II error of diffVar
c(rowMeans(dv.FDR.d1_new.50samples[3, ]),
  rowMeans(dv.FDR.d1_new.100samples[3, ]),
  rowMeans(dv.FDR.d1_new.125samples[3, ]),
  rowMeans(dv.FDR.d1_new.150samples[3, ]),
  rowMeans(dv.FDR.d1_new.200samples[3, ]),
  rowMeans(dv.FDR.d1_new.250samples[3, ]))

c(rowMeans(dv.type2.d1_new.50samples[3, ]),
  rowMeans(dv.type2.d1_new.100samples[3, ]),
  rowMeans(dv.type2.d1_new.125samples[3, ]),
  rowMeans(dv.type2.d1_new.150samples[3, ]),
  rowMeans(dv.type2.d1_new.200samples[3, ]),
  rowMeans(dv.type2.d1_new.250samples[3, ]))



### Computing times in seconds
dv.time.d1_new.50samples <- read.csv('DV_Time_diabetes1_50.csv')
dv.time.d1_new.100samples <- read.csv('DV_Time_diabetes1_100.csv')
dv.time.d1_new.125samples <- read.csv('DV_Time_diabetes1_125.csv')
dv.time.d1_new.150samples <- read.csv('DV_Time_diabetes1_150.csv')
dv.time.d1_new.200samples <- read.csv('DV_Time_diabetes1_200.csv')
dv.time.d1_new.250samples <- read.csv('DV_Time_diabetes1_250.csv')

dv.time.d1_new.50samples <- dv.time.d1_new.50samples[, -1]
dv.time.d1_new.100samples <- dv.time.d1_new.100samples[, -1]
dv.time.d1_new.125samples <- dv.time.d1_new.125samples[, -1]
dv.time.d1_new.150samples <- dv.time.d1_new.150samples[, -1]
dv.time.d1_new.200samples <- dv.time.d1_new.200samples[, -1]
dv.time.d1_new.250samples <- dv.time.d1_new.250samples[, -1]
 # Mean computing time
Time_dv_d1_new_2000genes <- matrix(NA, nrow = 6, ncol = 5)
Time_dv_d1_new_2000genes <- rbind(rowMeans(dv.time.d1_new.50samples),
                                  rowMeans(dv.time.d1_new.100samples),
                                  rowMeans(dv.time.d1_new.125samples),
                                  rowMeans(dv.time.d1_new.150samples),
                                  rowMeans(dv.time.d1_new.200samples),
                                  rowMeans(dv.time.d1_new.250samples))
Time_dv_d1_new_2000genes[ , 3] <- round(Time_dv_d1_new_2000genes[, 3], 1)
Time_dv_d1_new_2000genes[, c(1,2,4,5)] <- round(Time_dv_d1_new_2000genes[, c(1,2,4,5)])
row.names(Time_dv_d1_new_2000genes) <- c("n=50", "n=100", "n=125", "n=150", "n=200", "n=250")
colnames(Time_dv_d1_new_2000genes) <- c("clrDV","MDSeq","diffVar","GAMLSS-BH", "GAMLSS-BY")
Time_dv_d1_new_2000genes   # Average computing Times, in seconds

 # standard deviation of computing time
Time_sd_dv_d1_new <- matrix(NA, nrow = 6, ncol = 5)
Time_sd_dv_d1_new <- rbind(apply(dv.time.d1_new.50samples, 1, sd),
                           apply(dv.time.d1_new.100samples, 1, sd),
                           apply(dv.time.d1_new.125samples, 1, sd),
                           apply(dv.time.d1_new.150samples, 1, sd),
                           apply(dv.time.d1_new.200samples, 1, sd),
                           apply(dv.time.d1_new.250samples, 1, sd))
row.names(Time_sd_dv_d1_new) <- c("n=50", "n=100", "n=125", "n=150", "n=200", "n=250")
colnames(Time_sd_dv_d1_new) <- c("clrDV","MDSeq","diffVar","GAMLSS-BH", "GAMLSS-BY")
Time_sd_dv_d1_new[, c(1,2,4,5)] <- round(Time_sd_dv_d1_new[, c(1,2,4,5)])
Time_sd_dv_d1_new[, 3] <- round(Time_sd_dv_d1_new[,3], 2)
Time_sd_dv_d1_new



### Scatter plots
### 2*50 samples
plot(jitter(as.numeric(dv.FDR.d1_new.50samples[4, ]), amount = 0.002),
     jitter(as.numeric(dv.type2.d1_new.50samples[4, ]), amount = 0.002),
     lwd = 1, cex = 1,
     col = rainbow(7)[2], pch = 17,
     xlim = c(0, 0.35), ylim = c(0, 0.9),
     xlab = "FDR", ylab = "Type II Error")

points(jitter(as.numeric(dv.FDR.d1_new.50samples[5, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.50samples[5, ]), amount = 0.002),
       lwd =1, cex=1,
       col = "green", pch =0)

points(jitter(as.numeric(dv.FDR.d1_new.50samples[2, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.50samples[2, ]), amount = 0.002),
       col = "blue", lwd =1, cex=1, pch = 6 )

points(jitter(as.numeric(dv.FDR.d1_new.50samples[3, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.50samples[3, ]), amount = 0.002),
       lwd =1, cex=1,
       col = rainbow(7)[7], pch =5)

points(jitter(as.numeric(dv.FDR.d1_new.50samples[1, ]), amount = 0.002),
       jitter(as.numeric(dv.type2.d1_new.50samples[1, ]), amount = 0.002),
       lwd =1.25, cex=1.25, pch = 1,
       col = "red")

abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("bottomright", c("clrDV", "MDSeq", "diffVar", "GAMLSS-BH", "GAMLSS-BY"),
       pch = c(1, 6, 5 , 17, 0), pt.cex = c(1.25,1,1,1,1),
       cex= 0.7, pt.lwd = c(1.25, 1, 1, 1, 1),
       col = c("red", "blue", rainbow(7)[7],  rainbow(7)[2], "green"))
title(adj=0, "(a)")


### Scatter plots
### 2*100 samples
plot(jitter(as.numeric(dv.FDR.d1_new.100samples[4, ]), amount = 0.002),
     jitter(as.numeric(dv.type2.d1_new.100samples[4, ]), amount = 0.002),
     lwd = 1, cex = 1,
     col = rainbow(7)[2], pch = 17,
     xlim = c(0, 0.25), ylim = c(0, 0.45),
     xlab = "FDR", ylab = "Type II Error")

points(jitter(as.numeric(dv.FDR.d1_new.100samples[5, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.100samples[5, ]), amount = 0.002),
       lwd =1, cex=1,
       col = "green", pch =0)

points(jitter(as.numeric(dv.FDR.d1_new.100samples[2, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.100samples[2, ]), amount = 0.002),
       col = "blue", lwd =1, cex=1, pch = 6 )

points(jitter(as.numeric(dv.FDR.d1_new.100samples[3, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.100samples[3, ]), amount = 0.002),
       lwd =1, cex=1,
       col = rainbow(7)[7], pch =5)

points(jitter(as.numeric(dv.FDR.d1_new.100samples[1, ]), amount = 0.002),
       jitter(as.numeric(dv.type2.d1_new.100samples[1, ]), amount = 0.002),
       lwd =1.25, cex=1.25, pch = 1,
       col = "red")
abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)
legend("bottomright", c("clrDV", "MDSeq", "diffVar", "GAMLSS-BH", "GAMLSS-BY"),
       pch = c(1, 6, 5 , 17, 0), pt.cex = c(1.25,1,1,1,1),
       cex= 0.7, pt.lwd = c(1.25, 1, 1, 1, 1),
       col = c("red", "blue", rainbow(7)[7],  rainbow(7)[2], "green"))
title(adj=0, "(b)")


### Scatter plots
### 2*125 samples
plot(jitter(as.numeric(dv.FDR.d1_new.125samples[4, ]), amount = 0.002),
     jitter(as.numeric(dv.type2.d1_new.125samples[4, ]), amount = 0.002),
     lwd = 1, cex = 1,
     col = rainbow(7)[2], pch = 17,
     xlim = c(0, 0.25), ylim = c(0, 0.35),
     xlab = "FDR", ylab = "Type II Error")

points(jitter(as.numeric(dv.FDR.d1_new.125samples[5, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.125samples[5, ]), amount = 0.002),
       lwd =1, cex=1,
       col = "green", pch =0)

points(jitter(as.numeric(dv.FDR.d1_new.125samples[2, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.125samples[2, ]), amount = 0.002),
       col = "blue", lwd =1, cex=1, pch = 6 )

points(jitter(as.numeric(dv.FDR.d1_new.125samples[3, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.125samples[3, ]), amount = 0.002),
       lwd =1, cex=1,
       col = rainbow(7)[7], pch =5)

points(jitter(as.numeric(dv.FDR.d1_new.125samples[1, ]), amount = 0.002),
       jitter(as.numeric(dv.type2.d1_new.125samples[1, ]), amount = 0.002),
       lwd =1.25, cex=1.25, pch = 1,
       col = "red")
abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)

legend("bottomright", c("clrDV", "MDSeq", "diffVar", "GAMLSS-BH", "GAMLSS-BY"),
       pch = c(1, 6, 5 , 17, 0), pt.cex = c(1.25,1,1,1,1),
       cex= 0.7, pt.lwd = c(1.25, 1, 1, 1, 1),
       col = c("red", "blue", rainbow(7)[7],  rainbow(7)[2], "green"))
title(adj=0, "(c)")


### Scatter plots
### 2*150 samples
plot(jitter(as.numeric(dv.FDR.d1_new.150samples[4, ]), amount = 0.002),
     jitter(as.numeric(dv.type2.d1_new.150samples[4, ]), amount = 0.002),
     lwd = 1, cex = 1,
     col = rainbow(7)[2], pch = 17,
     xlim = c(0, 0.26), ylim = c(0, 0.26),
     xlab = "FDR", ylab = "Type II Error")

points(jitter(as.numeric(dv.FDR.d1_new.150samples[5, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.150samples[5, ]), amount = 0.002),
       lwd =1, cex=1,
       col = "green", pch =0)

points(jitter(as.numeric(dv.FDR.d1_new.150samples[2, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.150samples[2, ]), amount = 0.002),
       col = "blue", lwd =1, cex=1, pch = 6 )

points(jitter(as.numeric(dv.FDR.d1_new.150samples[3, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.150samples[3, ]), amount = 0.002),
       lwd =1, cex=1,
       col = rainbow(7)[7], pch =5)

points(jitter(as.numeric(dv.FDR.d1_new.150samples[1, ]), amount = 0.002),
       jitter(as.numeric(dv.type2.d1_new.150samples[1, ]), amount = 0.002),
       lwd =1.25, cex=1.25, pch = 1,
       col = "red")
abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)

legend("bottomright", c("clrDV", "MDSeq", "diffVar", "GAMLSS-BH", "GAMLSS-BY"),
       pch = c(1, 6, 5 , 17, 0), pt.cex = c(1.25,1,1,1,1),
       cex= 0.7, pt.lwd = c(1.25, 1, 1, 1, 1),
       col = c("red", "blue", rainbow(7)[7],  rainbow(7)[2], "green"))
title(adj=0, "(d)")


### Scatter plots
### 2*200 samples
plot(jitter(as.numeric(dv.FDR.d1_new.200samples[4, ]), amount = 0.002),
     jitter(as.numeric(dv.type2.d1_new.200samples[4, ]), amount = 0.002),
     lwd = 1, cex = 1,
     col = rainbow(7)[2], pch = 17,
     xlim = c(0, 0.22), ylim = c(0, 0.22),
     xlab = "FDR", ylab = "Type II Error")

points(jitter(as.numeric(dv.FDR.d1_new.200samples[5, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.200samples[5, ]), amount = 0.002),
       lwd =1, cex=1,
       col = "green", pch =0)

points(jitter(as.numeric(dv.FDR.d1_new.200samples[2, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.200samples[2, ]), amount = 0.002),
       col = "blue", lwd =1, cex=1, pch = 6 )

points(jitter(as.numeric(dv.FDR.d1_new.200samples[3, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.200samples[3, ]), amount = 0.002),
       lwd =1, cex=1,
       col = rainbow(7)[7], pch =5)

points(jitter(as.numeric(dv.FDR.d1_new.200samples[1, ]), amount = 0.002),
       jitter(as.numeric(dv.type2.d1_new.200samples[1, ]), amount = 0.002),
       lwd =1.25, cex=1.25, pch = 1,
       col = "red")
abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)

legend("topleft", c("clrDV", "MDSeq", "diffVar", "GAMLSS-BH", "GAMLSS-BY"),
       pch = c(1, 6, 5 , 17, 0), pt.cex = c(1.25,1,1,1,1),
       cex= 0.7, pt.lwd = c(1.25, 1, 1, 1, 1),
       col = c("red", "blue", rainbow(7)[7],  rainbow(7)[2], "green"))
title(adj=0, "(e)")


### Scatter plots
### 2*250 samples
plot(jitter(as.numeric(dv.FDR.d1_new.250samples[4, ]), amount = 0.002),
     jitter(as.numeric(dv.type2.d1_new.250samples[4, ]), amount = 0.002),
     lwd = 1, cex = 1,
     col = rainbow(7)[2], pch = 17,
     xlim = c(0, 0.24), ylim = c(0, 0.24),
     xlab = "FDR", ylab = "Type II Error")

points(jitter(as.numeric(dv.FDR.d1_new.250samples[5, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.250samples[5, ]), amount = 0.002),
       lwd =1, cex=1,
       col = "green", pch =0)

points(jitter(as.numeric(dv.FDR.d1_new.250samples[2, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.250samples[2, ]), amount = 0.002),
       col = "blue", lwd =1, cex=1, pch = 6 )

points(jitter(as.numeric(dv.FDR.d1_new.250samples[3, ]), amount = 0.002),
       jitter(as.numeric( dv.type2.d1_new.250samples[3, ]), amount = 0.002),
       lwd =1, cex=1,
       col = rainbow(7)[7], pch =5)

points(jitter(as.numeric(dv.FDR.d1_new.250samples[1, ]), amount = 0.002),
       jitter(as.numeric(dv.type2.d1_new.250samples[1, ]), amount = 0.002),
       lwd =1.25, cex=1.25, pch = 1,
       col = "red")
abline(h = 0.05, lty = 2, lwd = 1); abline(v = 0.05, lty = 2, lwd = 1)

legend("topleft", c("clrDV", "MDSeq", "diffVar", "GAMLSS-BH", "GAMLSS-BY"),
       pch = c(1, 6, 5 , 17, 0), pt.cex = c(1.25,1,1,1,1),
       cex= 0.7, pt.lwd = c(1.25, 1, 1, 1, 1),
       col = c("red", "blue", rainbow(7)[7],  rainbow(7)[2], "green"))
title(adj=0, "(f)")

###END###

