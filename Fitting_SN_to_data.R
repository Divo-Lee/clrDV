#####################################################
#clrDV: A differential variability test for
#RNA-Seq data based on the skew-normal distribution
#Author: Hongxiang Li
#Email: chelsea.divo@hotmail.com
#Date: 25 September 2022
#R Codes for fitting skew-normal model onto real-world data
#Part 1: Motivational Figure (Fig. 1)
#####################################################

  ## R packages downloaded
# install.packages("compositions")
# install.packages("devtools")
# install.packages("BiocManager")
# BiocManager::install("polyester")
# BiocManager::install("edgeR", force = T)
# devtools::install_github("Divo-Lee/clrDV")
# install.packages("sn")

library(polyester); library(compositions)
library(edgeR); library(sn); library(clrDV)

##############################################
### GSE123658_read_counts.gene_level.txt.gz  
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123658
### Valentim Data                            
##############################################
Data = read.table(gzfile('GSE123658_read_counts.gene_level.txt.gz'), sep="\t", header = T, check.names = TRUE, row.names = 1)
dim(Data)

 # filter
CPM <- cpm(Data)
keep <- (rowMeans(CPM[,1:43]) > 0.5  & rowMeans(CPM[,44:82]) > 0.5 & apply(Data[,1:43], 1, function(x) length(x[x==0])/length(x)) < 0.85 & apply(Data[,44:82], 1, function(x) length(x[x==0])/length(x)) < 0.85)
Data <- Data[keep, ]; dim(Data)

# without replicates
# h = healthy volunteers
# d1 = type 1 diabetic patients
# data_h <- Data[,1:43]
data_d1 <- Data[,44:82]


################################
### GSE150318_counts.csv.gz    
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150318
### Kelmer Data                
################################
Data2 = read.csv(gzfile('GSE150318_counts.csv.gz'), header = T, check.names = TRUE, row.names = 1)
dim(Data2)

data_10weeks <- Data2[ , seq(1,228,2)]
data_20weeks <- Data2[ , seq(2,228,2)]
Data2 <- cbind(data_10weeks, data_20weeks)
 # filter
CPM <- cpm(Data2)
keep <- (rowMeans(CPM[,1:114]) > 0.5 & rowMeans(CPM[,115:228]) > 0.5 & apply(Data2[,1:114], 1, function(x) length(x[x==0])/length(x)) < 0.85 & apply(Data2[,115:228], 1, function(x) length(x[x==0])/length(x)) < 0.85)
Data2 <- Data2[keep, ]

data_10weeks <- Data2[, 1:114]
# data_20weeks <- Data2[, 115:228]

                                                                                                                                                     

############################################################
# modified plot.selm() in "sn" package, for the name of xlab
plot.sn <- function(x, param.type="CP", which = c(1:4), caption,
                    panel = if (add.smooth) panel.smooth else points, main = "",
                    # sub.caption = NULL,
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
                    id.n = 3, labels.id = names(x@residuals.dp),
                    cex.id = 0.75, identline = TRUE, add.smooth = getOption("add.smooth"),
                    label.pos = c(4, 2), cex.caption = 1)
{
  if(!(is(x, "selm"))) stop("object not of class 'selm'")
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  dots <- list(...)
  nmdots <- names(dots)
  p <- slot(x, "size")["p"]
  if(missing(caption))  { caption <-  if(p> 1)
    c("Residuals vs Fitted Values",
      "Residual values and fitted error distribution",
      "Q-Q plot of (scaled DP residuals)^2",
      "P-P plot of (scaled DP residuals)^2") else
        c("Boxplot of observed values",
          "Empirical values and fitted distribution",
          "Q-Q plot of (scaled DP residuals)^2",
          "P-P plot of (scaled DP residuals)^2")}
  all.par <- slot(x, "param")
  param.type <- tolower(param.type)
  param <- all.par[[param.type]]
  if(is.null(param)) { message(paste(
    "Requested param.type='", param.type, "' evaluates to NULL.", sep=""))
    if(param.type == "pseudo-cp" & x@family== "SN")
      message("Pseudo-CP makes no sense for SN family")
    if(param.type == "cp" & x@family== "SC")
      message("CP makes no sense for SC family")
    if(param.type == "cp" & x@family== "ST")
      message("CP of ST family requires nu>4")
    stop("Consider another choice of param.type (DP or pseudo-CP)")
  }
  r <- residuals(x, param.type)
  r.lab <- paste(toupper(param.type), "residuals")
  dp <- if(length(all.par$fixed) > 0) all.par$dp.complete else all.par$dp
  nu. <- switch(x@family, ST = dp[p+3], SN = Inf, SC=1)
  rs <- slot(x,"residuals.dp")/dp[p+1]
  rs2 <- rs^2
  n <- slot(x, "size")["n.obs"]
  yh <- fitted(x, param.type)
  w <- weights(x)
  if (!is.null(w)) {
    wind <- (w != 0)
    r <- r[wind]
    yh <- yh[wind]
    w <- w[wind]
    labels.id <- labels.id[wind]
  }
  else w <- rep(1,n)
  rw <- n*w/slot(x,"size")["nw.obs"]
  cex.pts <- rw * if("cex" %in% nmdots) dots$cex else par("cex")
  if (is.null(id.n))
    id.n <- 0
  else {
    id.n <- as.integer(id.n)
    if (id.n < 0 || id.n > n)
      stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
  }
  if (id.n > 0) {
    if (is.null(labels.id))
      labels.id <- paste(1:n)
    iid <- 1:id.n
    # show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    show.rs <- sort.list(rs2, decreasing = TRUE)[iid]
    # rs2.lab <- paste("(scaled DP residuals)^2")
    text.id <- function(x, y, ind, adj.x = TRUE) {
      labpos <- if (adj.x)
        label.pos[1 + as.numeric(x > mean(range(x)))]
      else 3
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
           pos = labpos, offset = 0.25)
    }
  }
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (show[1]) {
    if(all(is.na(r)) & p>1)  message(paste("CP residuals not available;",
                                           "consider param.type='DP' or 'pseudo-CP'"))
    else {
      if(p == 1){
        y <-  (x@residuals.dp + x@fitted.values.dp)
        boxplot(y, plot=TRUE, col="gray85", border="gray60")
      }
      else { # p>1
        # if (id.n > 0)
        #    ylim <- extendrange(r = ylim, f = 0.08)
        ylim <- range(r, na.rm = TRUE)
        plot(yh, r, xlab = "Fitted values", ylab = r.lab, main = main,
             ylim = ylim, type = "n")
        panel(yh, r, ...)  # previously it included 'cex=cex.pts'
        # if (one.fig) title(sub = sub.caption, ...)
        if (id.n > 0) {
          y.id <- r[show.rs]
          y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
          text.id(yh[show.rs], y.id, show.rs)
        }
        abline(h = 0, lty = 2, col = "gray")
      } }
    mtext(caption[1], 3, 0.5, cex = cex.caption) }
  if (show[2]) {
    if(all(is.na(r)) & p>1) message(
      "CP residuals not available; consider param.type='DP' or 'pseudo-CP'")
    else {
      if (p == 1){
        y <-  (x@residuals.dp + x@fitted.values.dp)
        dp0 <- dp
        xlab="CLR-transformed count"}
      else {
        y <- r
        dp0 <- as.numeric(c(dp[1]-param[1], dp[-(1:p)]))
        xlab=r.lab
      }
      h <- hist(rep(y, w), plot=FALSE)
      extr <- extendrange(x=h$breaks)
      x.pts <- seq(max(extr), min(extr), length=501)
      d.fn <- get(paste("d", tolower(x@family), sep=""), inherits = TRUE)
      pdf <- d.fn(x.pts, dp=dp0)
      plot(c(h$mids, x.pts), c(h$density, pdf), type="n", main=main,
           xlab=xlab,  ylab="probability density")
      hist(rep(y, w), col="gray95", border="gray60", probability=TRUE,
           freq=FALSE, add=TRUE)
      lines(x.pts, pdf, ...)
      rug(y, ticksize=0.02, ...)
      # if (id.n > 0) {     rug(y, ticksize=0.015, ...)
      #   text(y[show.rs], 0, labels.id[show.rs], srt=90, cex=0.5, pos=1,
      #   offset=0.2) }
      mtext(caption[2], 3, 0.25, cex = cex.caption)
    }}
  if (show[3]) {
    ylim <- c(0, max(pretty(rs2)))
    q <- qf((1:n)/(n+1), 1, nu.)
    plot(q, sort(rs2), xlab="Theoretical values", ylab="Empirical values",
         ylim=ylim, type="p", main=main, ...)   # cex=cex.pts
    if(identline) abline(0, 1, lty = 2, col = "gray50")
    # if (one.fig) title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25, cex = cex.caption)
    if (id.n > 0) text.id(q[n+1-iid], rs2[show.rs], show.rs)
  }
  if (show[4]) {
    p <- (1:n)/(n+1)
    pr <- pf(sort(rs2), 1, nu.)
    plot(p, pr, xlab="Theoretical values", ylab="Empirical values",
         xlim=c(0,1), ylim=c(0,1), main=main, ...) # cex=cex.pts,
    if(identline) abline(0, 1, lty = 2, col = "gray50")
    # if (one.fig)  title(sub = sub.caption, ...)
    mtext(caption[4], 3, 0.25, cex = cex.caption)
    if(identline) abline(0, 1, lty = 2, col = "gray50")
    if (id.n > 0)  text.id(p[n+1-iid], pr[n+1-iid], show.rs)
  }
  # if (!one.fig && par("oma")[3] >= 1)
  #     mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}



 ### Fig. 1

##########################################
### Simulation One, Valentim Data, T1D ###
##########################################
 #
params <- get_params(data_d1)
N.genes <- 10000
N.samples <- 500
dat0 = create_read_numbers(params$mu,params$fit, params$p0,
                           m=N.genes, n = N.samples,
                           seed=123456)
row.names(dat0) <- paste('gene', 1:N.genes, sep='')

 # filter
CPM <- cpm(dat0)
keep <- (rowMeans(CPM[,1:500]) > 0.5  &  apply(dat0[,1:500], 1, function(x) length(x[x==0])/length(x)) < 0.85)
dat0 <- dat0[keep, ]; dim(dat0)

 # CLR-transformation
dat0[dat0 == 0] <- 1/2
clr_simu.1 <- t(clr(t(dat0)))
clr_simu.1 <- matrix(as.numeric(clr_simu.1), nrow = dim(dat0)[1], ncol = dim(dat0)[2])


 # plot clr-transformed counts of gene1 of simulated data
 # Fig. 1 (a)
fit_sn1 <- selm(clr_simu.1[1, ] ~ 1, family="SN")
plot.sn(fit_sn1, which=2, caption = NULL, main =NULL)
title(adj=0, "(a)")

 # MLE
clr.SN.fit(clr_simu.1[1, ])

 # we can check any other genes in simulated data,
 # almost all CLR-transformed simulated reads fit skew-normal model well
for (i in 2:100) {
  fit_sn <- selm(clr_simu.1[i, ] ~ 1, family="SN") 
  plot.sn(fit_sn, which=2, 
          caption = NULL, 
          main =paste('gene', i, sep=''))
}                                               


#######################################################
### Simulation Two, Kelmer data, at 10 weeks of age ###
#######################################################
 # simulate dataset
params.2 <- get_params(data_10weeks)
N.genes.2 <- 10000
N.samples.2 <- 500
dat2 = create_read_numbers(params.2$mu, params.2$fit, params.2$p0,
                           m=N.genes.2, n = N.samples.2,
                           seed=103)
row.names(dat2) <- paste('gene', 1:N.genes.2, sep='')

 # filter
CPM <- cpm(dat2)
keep <- (rowMeans(CPM[,1:N.samples.2]) > 0.5  &  apply(dat2[,1:N.samples.2], 1, function(x) length(x[x==0])/length(x)) < 0.85)
dat2 <- dat2[keep, ]; dim(dat2)


 # clr-transformation
dat2[dat2 == 0] <- 1/2
clr_simu.2 <- t(clr(t(dat2)))
clr_simu.2 <- matrix(as.numeric(clr_simu.2), nrow = dim(dat2)[1], ncol = dim(dat2)[2])

 # plot clr-transformed counts of gene1 of simulated data
 # Fig. 1 (b)
fit_sn2 <- selm(clr_simu.2[1, ] ~ 1, family="SN")
plot.sn(fit_sn2, which=2, caption = NULL,
        main =NULL, cex.id = 0)
title(adj=0, "(b)")

 # MLE
clr.SN.fit(clr_simu.2[1, ])

                                                       
 # we can check any other genes in simulated data,
 # almost all CLR-transformed simulated reads fit skew-normal model well                                                      
for (i in 2:100) {
  fit_sn2 <- selm(clr_simu.2[i, ] ~ 1, family="SN") 
  plot.sn(fit_sn2, which=2, 
          caption = NULL, 
          main =paste('gene', i, sep=''))
}
###END###


