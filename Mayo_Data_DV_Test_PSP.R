###################################################
#clrDV: A differential variability test for
#RNA-Seq data based on the skew-normal distribution
#Author: Hongxiang Li
#Email: chelsea.divo@hotmail.com
#Latest update: 23 Feb. 2023
#R Codes for DV test on Mayo RNA-Seq dataset
#Part 5: control vs. PSP comparison
###################################################

## R packages downloaded
# install.packages("readr")
# install.packages("devtools")
# install.packages("BiocManager")
# install.packages("compositions")
# BiocManager::install("polyester")
# BiocManager::install("edgeR", force = T)
# devtools::install_github("Divo-Lee/clrDV")
# install.packages("gamlss")
# BiocManager::install("cqn")
# devtools::install_github("zjdaye/MDSeq")
# install.packages("vioplot")
# install.packages("VennDiagram")
# install.packages("gridExtra")
# install.packages("httr")
# install.packages("jsonlite")
library(readr); library(compositions); library(edgeR)
library(MDSeq); library(gamlss); library(clrDV)
library(vioplot); library(VennDiagram); library(gridExtra)
library(httr); library(jsonlite)

### Read data
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


################
### DV Tests ###
################
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

### clrDV Test
t1 <- proc.time()
clr.counts2 <- clr.transform(data = mayo_counts_filter2)
s.d._test_results <- clrDV(clr.counts2, group = group3)
sum(is.na(s.d._test_results))
#s.d._test_results <- na.omit(s.d._test_results)
s.d._test_results <- as.data.frame(s.d._test_results)
dv_table_clrDV <- s.d._test_results[s.d._test_results$adj_pval < 0.05, ]
dv_genes_clrDV <- row.names(dv_table_clrDV)
as.numeric(proc.time() - t1)[3]
# additional information for downstream analysis
neg_log10_q <- -log10(s.d._test_results$adj_pval) 
SD_ratio <- s.d._test_results$sigma2/s.d._test_results$sigma1 # sigma.2/sigam.1
s.d._test <- cbind(s.d._test_results[1:5], SD_ratio, neg_log10_q)
dv_table_clrDV <- s.d._test[s.d._test$adj_pval < 0.05, ]
dim(s.d._test); dim(dv_table_clrDV); length(dv_genes_clrDV)



### MDSeq
data2 <- mayo_counts_filter2 # used by both MDSeq and GAMLSS
t2 <- proc.time()
libsizes2 <- colSums(data2)
nf2 <- calcNormFactors(data2, method="TMM")
els2 <- nf2 * libsizes2
sf2 <- els2 / exp(mean(log(libsizes2)))

contrasts2 <- get.model.matrix(as.factor(group3))
fit.MDSeq.dv <- MDSeq(data2, offsets = sf2,
                      contrast = contrasts2)
res.MDSeq.dv <- extract.ZIMD(fit.MDSeq.dv,
                             get='contrast',
                             compare=list(A="0",B="1"),
                             log2FC.threshold = 0)
res.MDSeq.dv <- na.omit(res.MDSeq.dv)
# sum(res.MDSeq.dv$FDR.dispersion < 0.05) # DV genes called
dv_table_MDSeq <- res.MDSeq.dv[res.MDSeq.dv$FDR.dispersion < 0.05, ]
dv_genes_MDSeq <- row.names(dv_table_MDSeq) # True positive DV genes
as.numeric(proc.time() - t2)[3]
dim(res.MDSeq.dv)
dim(mayo_counts_filter2)[1] - dim(res.MDSeq.dv)[1]  # the number of genes which got NA values
# there are 45 genes return NA values
length(dv_genes_MDSeq)



### GAMLSS
# Code modified from https://github.com/Vityay/ExpVarQuant/blob/master/ExpVarQuant.R
t4 <- proc.time()
design2 = model.matrix(~group3)
libsizes2 <- colSums(data2)
nf2 <- calcNormFactors(data2, method="TMM")
dat2.edgeR <- DGEList(counts=data2, norm.factors=nf2, group=group3)
dat2.edgeR <- estimateDisp(dat2.edgeR, design2)
dat2.edgeR$CPM <- cpm.DGEList(dat2.edgeR)
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
# BH FDR procedure
gamlss_NB$padj.cv <- p.adjust(gamlss_NB$p.cv, "BH")
dv_table_gamlss <- gamlss_NB[gamlss_NB$padj.cv < 0.05, ]

dv_genes_gamlss <- row.names(dv_table_gamlss)
as.numeric(proc.time() - t4)[3]
dim(gamlss_NB); length(dv_genes_gamlss)



################
### Analysis ###
################
## clrDV PSP DV genes list (Supplementary Table 4)
psp_clrDv_ratio_up_table <- dv_table_clrDV[dv_table_clrDV$SD_ratio > 1, ]
psp_clrDV_ratio_down_table <- dv_table_clrDV[dv_table_clrDV$SD_ratio < 1, ]
dim(psp_clrDv_ratio_up_table); dim(psp_clrDV_ratio_down_table)
# ranked
psp_clrDv_ratio_up_table <- psp_clrDv_ratio_up_table[order(psp_clrDv_ratio_up_table$SD_ratio),]
psp_clrDV_ratio_down_table <- psp_clrDV_ratio_down_table[order(psp_clrDV_ratio_down_table$SD_ratio),]
psp_clrDV_ratio_up_genes <- row.names(psp_clrDv_ratio_up_table)
psp_clrDV_ratio_down_genes <- row.names(psp_clrDV_ratio_down_table)
## Gene id convert to gene symbol
url = "https://www.biotools.fr/human/ensembl_symbol_converter/"
# up ratio
ids_up = psp_clrDV_ratio_up_genes
ids_json_up <- toJSON(ids_up)
body_up <- list(api=1, ids=ids_json_up)
r_up <- POST(url, body = body_up)
output_up <- fromJSON(content(r_up, "text"), flatten=TRUE)
output_up[sapply(output_up, is.null)] <- NA
gene_symbol <- unlist(output_up, use.names = FALSE)
psp_clrDv_ratio_up_table <- cbind.data.frame(gene_symbol, psp_clrDv_ratio_up_table)
psp_clrDv_ratio_up_table <- psp_clrDv_ratio_up_table[, c("gene_symbol", "SD_ratio", "adj_pval")]
psp_clrDv_ratio_up_table
# write.csv(psp_clrDv_ratio_up_table, "psp_DV_genes_up_ratio.csv", quote = F, row.names = T)

# down ratio
ids_down = psp_clrDV_ratio_down_genes
ids_json_down <- toJSON(ids_down)
body_down <- list(api=1, ids=ids_json_down)
r_down <- POST(url, body = body_down)
output_down <- fromJSON(content(r_down, "text"), flatten=TRUE)
output_down[sapply(output_down, is.null)] <- NA # some genes do not have gene symbol
gene_symbol <- unlist(output_down, use.names = FALSE)
for (i in 1:length(ids_down)) {
  if(is.na(gene_symbol[i])){
    gene_symbol[i] <- ids_down[i]
  }
} # replace the null gene symbol with gene id
psp_clrDV_ratio_down_table <- cbind.data.frame(gene_symbol, psp_clrDV_ratio_down_table)
psp_clrDV_ratio_down_table <- psp_clrDV_ratio_down_table[, c("gene_symbol", "SD_ratio", "adj_pval")]
# write.csv(psp_clrDV_ratio_down_table, "psp_DVgenes_down_ratio.csv", quote = F, row.names = T)
# head(psp_clrDV_ratio_down_table)


## Volcano plot, black and white with dashed line (q=0.05)
# Fig. S1 (b)
par(mar= c(5, 4.6, 4 ,1))
plot(log(SD_ratio,2), sqrt(neg_log10_q),
     xlim = c(-2, 2), ylim = c(0, 3.5),
     xlab = expression(paste(log[2]("SD ratio"))),
     ylab = expression(sqrt(-log[10](q))),
     col = c(rgb(0,0,0,0.3)), pch=1, lwd = 1.5)
abline(h = sqrt(-log10(0.05)), lty = 2, lwd = 1.25)
title(adj=0, "(b)")



### Comparison One, GAMLSS uses BH FDR
## Venn Diagram (BH)
# Fig. S2 (a)
clrDV = as.factor(dv_genes_clrDV)
GAMLSS = as.factor(dv_genes_gamlss)
MDSeq = as.factor(dv_genes_MDSeq)
Length_A<-length(clrDV)
Length_B<-length(GAMLSS)
Length_C<-length(MDSeq)
Length_AB<-length(intersect(clrDV,GAMLSS))
Length_BC<-length(intersect(GAMLSS,MDSeq))
Length_AC<-length(intersect(clrDV,MDSeq))
Length_ABC<-length(intersect(intersect(clrDV,GAMLSS),MDSeq))

vd <- venn.diagram(list(clrDV=clrDV,MDSeq=MDSeq,"GAMLSS-BH"=GAMLSS),
                   filename=NULL, lwd=1, lty=2,
                   col="transparent",
                   fill=c('red','green','cornflowerblue'),
                   cat.col="black",
                   height = 400, width = 400)
# grid.arrange(gTree(children=vd), top=textGrob(expression(bold("(a)")), x = 0.1,  hjust = 0))
grid.draw(vd)


## violin plots of the distribution of log2 SD ratio
# Fig. S2 (c)
length(setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq))

diff_clrDV1 <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
diff_gamlss1 <- setdiff(setdiff(dv_genes_gamlss, dv_genes_clrDV), dv_genes_MDSeq)
diff_MDSeq1 <- setdiff(setdiff(dv_genes_MDSeq, dv_genes_clrDV), dv_genes_gamlss)

cgm_genes1 <- intersect(intersect(dv_genes_clrDV,dv_genes_MDSeq), dv_genes_gamlss)
length(cgm_genes1)

vioplot(log2(dv_table_clrDV[cgm_genes1, 6]),
        log2(dv_table_clrDV[diff_clrDV1, 6]),
        log2(s.d._test[diff_gamlss1,6]),
        log2(s.d._test[diff_MDSeq1,6]),
        names=c("cgm", "c-g-m", "g-c-m","m-c-g"),
        col = 0, border = "black", pchMed="",
        rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
        ylab = expression(paste(log[2]("SD ratio"))),
        horizontal = F)
stripchart(list(log2(dv_table_clrDV[cgm_genes1, 6]),
                log2(dv_table_clrDV[diff_clrDV1, 6]),
                log2(s.d._test[diff_gamlss1,6]),
                log2(s.d._test[diff_MDSeq1,6]) ),
           vertical = T, pch = "." , cex = 3,
           method = "jitter", add = T, lwd = 3,
           col = c(rgb(0,0,0,0.1)))
title(adj=0, "(c)")


### ### ### ### ### ###
## Supplementary Table 6
### union_DV_table_BH (GAMLSS use BH FDR)
union_dv_genes_psp_BH <- union(union(dv_genes_clrDV, dv_genes_MDSeq), dv_genes_gamlss)
length(union_dv_genes_psp_BH)

union_dv_table_psp_BH <- s.d._test[union_dv_genes_psp_BH, c("SD_ratio", "adj_pval")]
indicator_clrDV <- as.integer(union_dv_genes_psp_BH %in% dv_genes_clrDV)
indicator_MDSeq <- as.integer(union_dv_genes_psp_BH %in% dv_genes_MDSeq)
indicator_gamlss_BH <- as.integer(union_dv_genes_psp_BH %in% dv_genes_gamlss)

union_dv_table_psp_BH <- cbind(union_dv_table_psp_BH, indicator_clrDV,
                               indicator_MDSeq, indicator_gamlss_BH)
# gene id convert to gene symbol
ids = union_dv_genes_psp_BH
ids_json <- toJSON(ids)
body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)
output <- fromJSON(content(r, "text"), flatten=TRUE)
output[sapply(output, is.null)] <- NA # some genes do not have gene symbol
gene_symbol <- unlist(output, use.names = FALSE)
for (i in 1:length(ids)) {
  if(is.na(gene_symbol[i])){
    gene_symbol[i] <- ids[i]
  }
}
union_dv_table_psp_BH <- cbind.data.frame(gene_symbol, union_dv_table_psp_BH)
union_dv_table_psp_BH <- union_dv_table_psp_BH[order(union_dv_table_psp_BH$SD_ratio) , ]
# write.csv(union_dv_table_psp_BH, "union_DV_table_psp_BH.csv", quote = F, row.names = T)
# head(union_dv_table_psp_BH)

## unique dv genes detected by clrDV (GAMLSS use BH FDR)
unique_clrDV_dv_genes_psp_BH <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
unique_clrDV_dv_table_psp_BH <- union_dv_table_psp_BH[unique_clrDV_dv_genes_psp_BH, c(1,2,3)]
unique_clrDV_dv_table_psp_BH <- unique_clrDV_dv_table_psp_BH[order(unique_clrDV_dv_table_psp_BH[,2]),]
# write.csv(unique_clrDV_dv_table_psp_BH, "unique_clrDV_DV_table_psp_BH.csv", quote = F, row.names = T)
# head(unique_clrDV_dv_table_psp_BH)


######################################
### Comparison Two, GAMLSS uses BY FDR
gamlss_NB$padj.cv <- p.adjust(gamlss_NB$p.cv, "BY")
dv_table_gamlss <- gamlss_NB[gamlss_NB$padj.cv < 0.05, ]
dim(dv_table_gamlss)
dv_genes_gamlss <- row.names(dv_table_gamlss)

## Venn Diagram
# Fig. S2 (b)
clrDV = as.factor(dv_genes_clrDV)
GAMLSS = as.factor(dv_genes_gamlss)
MDSeq = as.factor(dv_genes_MDSeq)
Length_A<-length(clrDV)
Length_B<-length(GAMLSS)
Length_C<-length(MDSeq)
Length_AB<-length(intersect(clrDV,GAMLSS))
Length_BC<-length(intersect(GAMLSS,MDSeq))
Length_AC<-length(intersect(clrDV,MDSeq))
Length_ABC<-length(intersect(intersect(clrDV,GAMLSS),MDSeq))

vd2 <- venn.diagram(list(clrDV=clrDV,MDSeq=MDSeq,"GAMLSS-BY"=GAMLSS),
                    filename=NULL, lwd=1, lty=2,
                    col="transparent",
                    fill=c('red','green','cornflowerblue'),
                    cat.col="black",
                    height = 400, width = 400)
# grid.arrange(gTree(children=vd2), top=textGrob(expression(bold("(b)")), x = 0.1, hjust = 0))
grid.draw(vd2)


## violin plots of the distribution of estimated log2 SD ratio
# Fig. S2 (d)
cgm_genes2 <- intersect(intersect(dv_genes_clrDV,dv_genes_MDSeq), dv_genes_gamlss)
length(cgm_genes2)

diff_clrDV2 <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
diff_gamlss2 <- setdiff(setdiff(dv_genes_gamlss, dv_genes_clrDV), dv_genes_MDSeq)
diff_MDSeq2 <- setdiff(setdiff(dv_genes_MDSeq, dv_genes_clrDV), dv_genes_gamlss)

vioplot(log2(dv_table_clrDV[cgm_genes2, 6]),
        log2(dv_table_clrDV[diff_clrDV2, 6]),
        log2(s.d._test[diff_gamlss2,6]),
        log2(s.d._test[diff_MDSeq2,6]),
        names=c("cgm", "c-g-m", "g-c-m","m-c-g"),
        col = 0, border = "black", pchMed="",
        rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
        ylab = expression(paste(log[2]("SD ratio"))),
        horizontal = F)
stripchart(list(log2(dv_table_clrDV[cgm_genes2, 6]),
                log2(dv_table_clrDV[diff_clrDV2, 6]),
                log2(s.d._test[diff_gamlss2,6]),
                log2(s.d._test[diff_MDSeq2,6]) ),
           vertical = T, pch ="." , cex = 3,
           method = "jitter", add = T, lwd =1,
           col = c(rgb(0,0,0,0.1)))
title(adj=0, "(d)")


########################
## Supplementary Table 6
### union_DV_table_BY (GAMLSS use BY FDR)
union_dv_genes_psp_BY <- union(union(dv_genes_clrDV, dv_genes_MDSeq), dv_genes_gamlss)
length(union_dv_genes_psp_BY)

union_dv_table_psp_BY <- s.d._test[union_dv_genes_psp_BY,  c("SD_ratio", "adj_pval")]
indicator_clrDV <- as.integer(union_dv_genes_psp_BY %in% dv_genes_clrDV)
indicator_MDSeq <- as.integer(union_dv_genes_psp_BY %in% dv_genes_MDSeq)
indicator_gamlss_BY <- as.integer(union_dv_genes_psp_BY %in% dv_genes_gamlss)

union_dv_table_psp_BY <- cbind(union_dv_table_psp_BY, indicator_clrDV,
                               indicator_MDSeq, indicator_gamlss_BY)
# convert gene id to gene symbol
ids = union_dv_genes_psp_BY
ids_json <- toJSON(ids)
body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)
output <- fromJSON(content(r, "text"), flatten=TRUE)
output[sapply(output, is.null)] <- NA # some genes do not have gene symbol
gene_symbol <- unlist(output, use.names = FALSE)
for (i in 1:length(ids)) {
  if(is.na(gene_symbol[i])){
    gene_symbol[i] <- ids[i]
  }
}
union_dv_table_psp_BY <- cbind.data.frame(gene_symbol, union_dv_table_psp_BY)
union_dv_table_psp_BY <- union_dv_table_psp_BY[order(union_dv_table_psp_BY$SD_ratio), ]
# write.csv(union_dv_table_psp_BY, "union_DV_table_psp_BY.csv", quote = F, row.names = T)
# head(union_dv_table_psp_BY)

# unique dv genes detected by clrDV
unique_clrDV_dv_genes_psp_BY <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
unique_clrDV_dv_table_psp_BY <- union_dv_table_psp_BY[unique_clrDV_dv_genes_psp_BY, c(1,2,3)]
dim(unique_clrDV_dv_table_psp_BY)
unique_clrDV_dv_table_psp_BY <- unique_clrDV_dv_table_psp_BY[order(unique_clrDV_dv_table_psp_BY$SD_ratio), ]
# write.csv(unique_clrDV_dv_table_psp_BY, "unique_clrDV_DV_table_psp_BY.csv", quote = F, row.names = T)
# head(unique_clrDV_dv_table_psp_BY)

###END###
