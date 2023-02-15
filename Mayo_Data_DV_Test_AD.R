#######################################################
#clrDV: A differential variability test for
#RNA-Seq data based on the skew-normal distribution
#Author: Hongxiang Li
#Email: chelsea.divo@hotmail.com
#Latest update: 25 September 2022
#R Codes for DV test on Mayo RNA-Seq dataset
#Part 4: control vs. AD comparison
#######################################################

 ## R packages downloaded
# install.packages("readr")
# install.packages("BiocManager")
# install.packages("devtools")
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

# raw count table
ad_counts <- read.csv('MayoRNAseq_RNAseq_TCX_geneCounts.csv')
row.names(ad_counts) <- ad_counts$ensembl_id
ad_counts <- ad_counts[, -1]; dim(ad_counts)
# meta-data
ad_meta <- read.csv('MayoRNAseq_individual_metadata_031422.csv')
sum(is.na(ad_meta$diagnosis)) # check NA in disease name column
# remove samples which give NA in disease name column, in meta-data
ad_meta <- ad_meta[-which(is.na(ad_meta$diagnosis)), ]
 
# samples id in meta-data
control_id <- (ad_meta$individualID)[ad_meta$diagnosis == "control"]
ad_id <- (ad_meta$individualID)[ad_meta$diagnosis == "Alzheimer Disease"]
# samples id in raw count table
ad_counts_id <-  parse_number(colnames(ad_counts))

control_vector <- sapply(ad_counts_id, function(k) k %in% control_id)
control_counts <-  ad_counts[, control_vector]
ad_vector <- sapply(ad_counts_id, function(k) k %in% ad_id)
ad_counts1 <- ad_counts[, ad_vector]
dim(control_counts); dim(ad_counts1)
N_ad <- length(control_counts) + length(ad_counts1)
mayo_counts1 <- as.matrix(cbind(control_counts, ad_counts1), 
nrows=dim(ad_counts)[1], ncol = N_ad)

## Filter
CPM2 <- cpm(mayo_counts1)

keep <- rowMeans(CPM2[,1:length(control_counts)]) > 0.5 & 
rowMeans(CPM2[,(length(control_counts)+1):N_ad]) > 0.5 &
apply(mayo_counts1[,1:length(control_counts)], 1, function(k) length(k[k == 0])/length(k)) < 0.85 &
apply(mayo_counts1[,(length(control_counts)+1):N_ad], 1, function(k) length(k[k == 0])/length(k)) < 0.85

mayo_counts_filter <- mayo_counts1[keep, ]
dim(mayo_counts_filter)

################
### DV Tests ###
################
group2 = c(rep(0, length(control_counts)), rep(1, length(ad_counts1))) # 78 control, 82 AD

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

##############
### clrDV Test
t1 <- proc.time()
clr.counts2 <- clr.transform(data = mayo_counts_filter)
s.d._test_results <- clrDV(clr.counts2, group = group2)
sum(is.na(s.d._test_results)) # check NA value
s.d._test_results <- na.omit(s.d._test_results)
dv_table_clrDV <- s.d._test_results[s.d._test_results[, 5] < 0.05, ]
dv_genes_clrDV <- row.names(dv_table_clrDV)
as.numeric(proc.time() - t1)[3]
  # additional information for downstream analysis
neg_log10_q <- -log10(s.d._test_results[, 5]) # q-value is the adjusted p-value
SD_ratio <- s.d._test_results[,10]/s.d._test_results[, 6] # sigma.2/sigam.1
s.d._test_results <- cbind(s.d._test_results, SD_ratio, neg_log10_q)
s.d._test <- s.d._test_results[, c(1,2,3,4,5,14,15)]
dv_table_clrDV <- s.d._test[s.d._test[, 5] < 0.05, ]
dim(s.d._test); dim(dv_table_clrDV); length(dv_genes_clrDV)



#########
### MDSeq
data2 <- mayo_counts_filter # used by both MDSeq and GAMLSS
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
dv_table_MDSeq <- res.MDSeq.dv[res.MDSeq.dv$FDR.dispersion < 0.05, ]
dv_genes_MDSeq <- row.names(dv_table_MDSeq) # DV genes called
as.numeric(proc.time() - t2)[3]
dim(res.MDSeq.dv); length(dv_genes_MDSeq)



##########
### GAMLSS
# Code modified from https://github.com/Vityay/ExpVarQuant/blob/master/ExpVarQuant.R
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
gamlss_NB$padj.cv <- p.adjust(gamlss_NB$p.cv, "BH") # BH FDR procedure
dv_table_gamlss <- gamlss_NB[gamlss_NB$padj.cv < 0.05, ]
dv_genes_gamlss <- row.names(dv_table_gamlss)
as.numeric(proc.time() - t4)[3]
dim(gamlss_NB); dim(dv_table_gamlss); length(dv_genes_gamlss)



########################
### DV Test Analysis ###
########################

## clrDV AD DV genes list (Supplementary Table 3)
clrDv_ratio_up_table <- dv_table_clrDV[dv_table_clrDV[, 6] > 1, ]
clrDV_ratio_down_table <- dv_table_clrDV[dv_table_clrDV[, 6] < 1, ]
dim(clrDv_ratio_up_table); dim(clrDV_ratio_down_table)

clrDv_ratio_up_table <- clrDv_ratio_up_table[order(clrDv_ratio_up_table[,6]),]
clrDV_ratio_down_table <- clrDV_ratio_down_table[order(clrDV_ratio_down_table[,6]),]
clrDV_ratio_up_genes <- row.names(clrDv_ratio_up_table)
clrDV_ratio_down_genes <- row.names(clrDV_ratio_down_table)
 ## Gene id convert to gene symbol
url = "https://www.biotools.fr/human/ensembl_symbol_converter/"
 # up ratio
ids_up = clrDV_ratio_up_genes
ids_json_up <- toJSON(ids_up)
body_up <- list(api=1, ids=ids_json_up)
r_up <- POST(url, body = body_up)
output_up <- fromJSON(content(r_up, "text"), flatten=TRUE)
output_up[sapply(output_up, is.null)] <- NA
gene_symbol <- unlist(output_up, use.names = FALSE)
clrDv_ratio_up_table <- cbind.data.frame(gene_symbol, clrDv_ratio_up_table)
clrDv_ratio_up_table <- clrDv_ratio_up_table[, c(1,7,6)]
# clrDv_ratio_up_table
# write.csv(clrDv_ratio_up_table, "AD_DVgenes_up_ratio.csv")

 # down ratio
ids_down = clrDV_ratio_down_genes
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
clrDV_ratio_down_table <- cbind.data.frame(gene_symbol, clrDV_ratio_down_table)
clrDV_ratio_down_table <- clrDV_ratio_down_table[, c(1,7,6)]
 # write.csv(clrDV_ratio_down_table, "AD_DVgenes_down_ratio.csv")
 # head(clrDV_ratio_down_table)


#########################
## volcano plot for clrDV
 # Fig. S1 (a)
par(mar= c(5, 4.6, 4 ,1))
plot(log(SD_ratio,2), sqrt(neg_log10_q),
     xlim = c(-2, 2), ylim = c(0, 3.5),
     xlab = expression(paste(log[2]("SD ratio"))),
     ylab = expression(sqrt(-log[10](q))),
     col = c(rgb(0,0,0,0.3)), pch=1, lwd = 1.5)
abline(h = sqrt(-log10(0.05)), lty = 2, lwd = 1.25)
title(adj=0, "(a)")



######################################
### Comparison One, GAMLSS uses BH FDR
  # DV genes uniquely detected
diff_clrDV1 <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
diff_gamlss1 <- setdiff(setdiff(dv_genes_gamlss, dv_genes_clrDV), dv_genes_MDSeq)
diff_MDSeq1 <- setdiff(setdiff(dv_genes_MDSeq, dv_genes_clrDV), dv_genes_gamlss)
  # overlap
cgm_genes1 <- intersect(intersect(dv_genes_clrDV,dv_genes_MDSeq), dv_genes_gamlss)

## Venn Diagram (BH)
 # Fig. 4 (a)
clrDV = as.factor(dv_genes_clrDV)
GAMLSS_BH = as.factor(dv_genes_gamlss)
MDSeq = as.factor(dv_genes_MDSeq)
Length_A <- length(clrDV)
Length_B <- length(GAMLSS_BH)
Length_C <- length(MDSeq)
Length_AB <- length(intersect(clrDV,GAMLSS_BH))
Length_BC <- length(intersect(GAMLSS_BH,MDSeq))
Length_AC <- length(intersect(clrDV,MDSeq))
Length_ABC <- length(intersect(intersect(clrDV,GAMLSS_BH),MDSeq))

vd <- venn.diagram(list(clrDV=clrDV,MDSeq=MDSeq,"GAMLSS-BH"=GAMLSS_BH),
                  filename=NULL, lwd=1, lty=2,
                  col="transparent",
                  fill=c('red','green','cornflowerblue'),
                  cat.col="black",
                  height = 400, width = 400)
 # grid.arrange(gTree(children=vd), top=grid::textGrob(expression(bold("(a)")), x = 0.1,  hjust = 0))
grid.draw(vd)


 # Fig. 4 (c)
vioplot(log2(s.d._test[cgm_genes1, 6]),
        log2(s.d._test[diff_clrDV1, 6]),
        log2(s.d._test[diff_gamlss1, 6]),
        log2(s.d._test[diff_MDSeq1, 6]),
        names=c("cgm", "c-g-m", "g-c-m","m-c-g"),
        col = c(0,0,0,0), border ="black", pchMed="",
        rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
        ylab = expression(paste(log[2]("SD ratio"))),
        horizontal = F)
stripchart(list(log2(s.d._test[cgm_genes1, 6]),
                log2(s.d._test[diff_clrDV1, 6]),
                log2(s.d._test[diff_gamlss1, 6]),
                log2(s.d._test[diff_MDSeq1, 6]) ),
           vertical = T, pch = "." , cex = 3.5,
           method = "jitter", add = T, lwd = 2.5,
           col = c(rgb(0,0,0,0.1)))
title(adj=0, "(c)")



###############
## (Supplementary Table 5)
## union_DV_table_BH (GAMLSS use BH FDR) 
union_dv_genes_AD_BH <- union(union(dv_genes_clrDV, dv_genes_MDSeq),
                              dv_genes_gamlss)
length(union_dv_genes_AD_BH)

union_dv_table_AD_BH <- s.d._test_results[union_dv_genes_AD_BH, c(14,5)]
indicator_clrDV <- as.integer(union_dv_genes_AD_BH %in% dv_genes_clrDV)
indicator_MDSeq <- as.integer(union_dv_genes_AD_BH %in% dv_genes_MDSeq)
indicator_gamlss_BH <- as.integer(union_dv_genes_AD_BH %in% dv_genes_gamlss)

union_dv_table_AD_BH <- cbind(union_dv_table_AD_BH, indicator_clrDV,
                              indicator_MDSeq, indicator_gamlss_BH)
# gene id convert to gene symbol
ids = union_dv_genes_AD_BH
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
union_dv_table_AD_BH <- cbind.data.frame(gene_symbol, union_dv_table_AD_BH)
 # ranked by SD ratio
union_dv_table_AD_BH <- union_dv_table_AD_BH[order(union_dv_table_AD_BH[ ,2]) , ]
 # write.csv(union_dv_table_AD_BH, "union_dv_table_AD_BH.csv")
 # head(union_dv_table_AD_BH)

# unique dv genes detected by clrDV
 # BH means that GAMLSS used BH FDR
unique_clrDV_dv_genes_AD_BH <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
unique_clrDV_dv_table_AD_BH <- union_dv_table_AD_BH[unique_clrDV_dv_genes_AD_BH, 1:3]
unique_clrDV_dv_table_AD_BH <- unique_clrDV_dv_table_AD_BH[order(unique_clrDV_dv_table_AD_BH[,2]),]
unique_clrDV_dv_table_AD_BH
 # write.csv(unique_clrDV_dv_table_AD_BH, "unique_clrDV_dv_table_AD_BH.csv")
 


################################################
### Comparison Two, GAMLSS uses BY FDR procedure
gamlss_NB$padj.cv <- p.adjust(gamlss_NB$p.cv, "BY")
dv_table_gamlss <- gamlss_NB[gamlss_NB$padj.cv < 0.05, ]
dim(dv_table_gamlss)
dv_genes_gamlss <- row.names(dv_table_gamlss)
length(setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq))

## Venn Diagram (BY)
 # Fig. 4 (b)
clrDV = as.factor(dv_genes_clrDV)
GAMLSS_BY = as.factor(dv_genes_gamlss)
MDSeq = as.factor(dv_genes_MDSeq)
Length_A<-length(clrDV)
Length_B<-length(GAMLSS_BY)
Length_C<-length(MDSeq)
Length_AB<-length(intersect(clrDV,GAMLSS_BY))
Length_BC<-length(intersect(GAMLSS_BY,MDSeq))
Length_AC<-length(intersect(clrDV,MDSeq))
Length_ABC<-length(intersect(intersect(clrDV,GAMLSS_BY),MDSeq))

vd2 <- venn.diagram(list(clrDV=clrDV,MDSeq=MDSeq,"GAMLSS-BY"=GAMLSS_BY),
                  filename=NULL, lwd=1, lty=2,
                  col="transparent",
                  fill=c('red','green','cornflowerblue'),
                  cat.col="black",
                  #                  reverse=TRUE,
                  height = 400, width = 400)
 # grid.arrange(gTree(children=vd2), top=grid::textGrob(expression(bold("(b)")), x = 0.1,  hjust = 0))
grid.draw(vd2)


## violin plots
# Fig. 4 (d)
diff_clrDV2 <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss), dv_genes_MDSeq)
diff_gamlss2 <- setdiff(setdiff(dv_genes_gamlss, dv_genes_clrDV), dv_genes_MDSeq)
diff_MDSeq2 <- setdiff(setdiff(dv_genes_MDSeq, dv_genes_clrDV), dv_genes_gamlss)
 # overlap
cgm_genes2 <- intersect(intersect(dv_genes_clrDV,dv_genes_MDSeq), dv_genes_gamlss)
length(cgm_genes2)

vioplot(log2(s.d._test[cgm_genes2, 6]),
        log2(s.d._test[diff_clrDV2, 6]),
        log2(s.d._test[diff_gamlss2, 6]),
        log2(s.d._test[diff_MDSeq2, 6]),
        names=c("cgm", "c-g-m", "g-c-m","m-c-g"),
        col = c(0,0,0,0), border ="black", pchMed="",
        rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
        ylab = expression(paste(log[2]("SD ratio"))),
        horizontal = F)
stripchart(list(log2(s.d._test[cgm_genes2, 6]),
                log2(s.d._test[diff_clrDV2, 6]),
                log2(s.d._test[diff_gamlss2,6]),
                log2(s.d._test[diff_MDSeq2,6]) ),
           vertical = T, pch ="." , cex = 3.5,
           method = "jitter", add = T, lwd =3,
           col = c(rgb(0,0,0,0.1)))
title(adj=0, "(d)")



### ### ### ### ###
### union_DV_table_BY (GAMLSS use BY FDR) 
## (Supplementary Table 5)
dv_genes_gamlss_BY <- row.names(dv_table_gamlss)
union_dv_genes_AD_BY <- union(union(dv_genes_clrDV, dv_genes_MDSeq),
                              dv_genes_gamlss_BY)
length(union_dv_genes_AD_BY)

union_dv_table_AD_BY <- s.d._test_results[union_dv_genes_AD_BY, c(14,5)]
indicator_clrDV <- as.integer(union_dv_genes_AD_BY %in% dv_genes_clrDV)
indicator_MDSeq <- as.integer(union_dv_genes_AD_BY %in% dv_genes_MDSeq)
indicator_gamlss_BY <- as.integer(union_dv_genes_AD_BY %in% dv_genes_gamlss_BY)

union_dv_table_AD_BY <- cbind(union_dv_table_AD_BY, indicator_clrDV,
                              indicator_MDSeq, indicator_gamlss_BY)
  # convert gene id to gene symbol
ids = union_dv_genes_AD_BY
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
union_dv_table_AD_BY <- cbind.data.frame(gene_symbol, union_dv_table_AD_BY)
union_dv_table_AD_BY <- union_dv_table_AD_BY[order(union_dv_table_AD_BY[, 2]), ]
 # write.csv(union_dv_table_AD_BY, "union_dv_table_AD_BY.csv")
 # head(union_dv_table_AD_BY)

# unique dv genes detected by clrDV
unique_clrDV_dv_genes_AD_BY <- setdiff(setdiff(dv_genes_clrDV, dv_genes_gamlss_BY), dv_genes_MDSeq)
unique_clrDV_dv_table_AD_BY <- union_dv_table_AD_BY[unique_clrDV_dv_genes_AD_BY, c(1,2,3)]
dim(unique_clrDV_dv_table_AD_BY)
unique_clrDV_dv_table_AD_BY <- unique_clrDV_dv_table_AD_BY[order(unique_clrDV_dv_table_AD_BY[,2]),]
 # write.csv(unique_clrDV_dv_table_AD_BY, "unique_clrDV_dv_table_AD_BY.csv")
 # head(unique_clrDV_dv_table_AD_BY)


########################################
## Violin plots for large-magnitude SD ratio DV genes
## Fig. S3
top_dv_genes_to_plot <- c("ENSG00000142494", "ENSG00000119147",
                          "ENSG00000124107", "ENSG00000119681",
                          "ENSG00000203618", "ENSG00000141456")
# ENSG00000142494	SLC47A1
# ENSG00000119147	C2orf40
# ENSG00000124107	SLPI
# ENSG00000119681	LTBP2
# ENSG00000203618	GP1BB
# ENSG00000141456	PELP1
labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")

par(mfrow=c(3,2))
for (i in 1:6) {
  vioplot(clr.counts2[top_dv_genes_to_plot[i], ][1:78],
          clr.counts2[top_dv_genes_to_plot[i], ][79:160],
          names=c("control", "AD"),  pchMed="",
          col = c(0,0), border ="black", horizontal = T,
          rectCol=rgb(0,0,0,0), lineCol=rgb(0,0,0,0),
          xlab="CLR-transformed count", ylab="Condition")

  stripchart(list(clr.counts2[top_dv_genes_to_plot[i], ][1:78],
                  clr.counts2[top_dv_genes_to_plot[i], ][79:160]),
             method="jitter", vertical=F,  add=TRUE,
             pch=1, cex=0.5, col="black")
  title(adj=0, labels[i])
}


###END###
