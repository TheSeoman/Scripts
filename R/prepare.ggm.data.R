source('Scripts/R/paths.R')
source('Scripts/R/residuals.R')

require('GenomicRanges')
require('illuminaHumanv3.db')

library(Rsamtools)

load(PATHS$HERV.MEQTL.OVERLAP.DATA)
load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$METH.COV.MATRIX.DATA)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$SNP.TFBS.OVERLAP.DATA)
load(PATHS$EXPR.TFBS.OVERLAP.DATA)
load(PATHS$METH.TFBS.OVERLAP.DATA)

load(PATHS$HERV.METH.OVERLAP.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)


load(PATHS$SNP.RANGES.DATA)

snp.samples <- scan(PATHS$F.SNP.SAMPLES)

get.snp.data <- function(snp.range) {
  data = scanTabix(PATHS$F.SNP, param=snp.range)
  snp.data.list <- lapply(data, function(x) strsplit(x, '\t'))
  snp.data.table <- data.frame(matrix(unlist(snp.data.list), ncol=length(snp.data.list), byrow=FALSE))
  snp.data.table <- snp.data.table[-(1:5), , drop = FALSE]
  colnames(snp.data.table) <- names(snp.range)
  rownames(snp.data.table) <- snp.samples
  return(snp.data.table)
}

covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)
id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                         & covariates.all$axio_s4f4 %in% snp.samples 
                         & covariates.all$meth_f4 %in% rownames(meth.matrix), c('expr_s4f4ogtt', 'axio_s4f4', 'meth_f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]
id.map$expr_s4f4ogtt <- as.character(id.map$expr_s4f4ogtt)
id.map$axio_s4f4 <- as.character(id.map$axio_s4f4)
id.map$meth_f4 <- as.character(id.map$meth_f4)

genes <- unlist(as.list(illuminaHumanv3SYMBOL))
genes <- genes[!is.na(genes)]

set <- 'hervS1'
filter <- 'snp'
snp.count.threshold <- 5
meqtl.pairs <- get(paste0(set, '.meqtl.trans.overlap'))[[filter]]
meqtl.count <- table(meqtl.pairs$snp)[table(meqtl.pairs$snp) > 0]
snps <- names(meqtl.count[meqtl.count >= snp.count.threshold])

total.meth.ids <- unique(meqtl.pairs[meqtl.pairs$snp %in% snps, 'cpg'])

eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('cis.', filter)]]

# snps in meQTLs with snp and meth probe in herv
snps <- intersect(as.character(unique(meqtl.pairs$snp)), as.character(unique(eqtl.pairs$snps)))


for(snp in snps[1:5]) {
  cat(paste0('Processing snp: ', snp), fill = T)
  snp.tfbs.ids <- snp.tfbs.overlap$pairs[snp.tfbs.overlap$pairs$snp.id == snp, 'tfbs.id']
  snp.tfbs.genes <- unique(snp.tfbs.overlap$tfbs.ranges[snp.tfbs.ids]$TF)
  cat(paste0('tfs overlapping with snp: ', paste(snp.tfbs.genes, collapse=', ') ), fill = T)
  meth.ids <- as.character(meqtl.pairs[meqtl.pairs$snp == snp, 'cpg'])
  cat(paste0('cpgs associated with snp: ', paste(meth.ids, sep=', ') ), fill = T)
  meth.tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% meth.ids, 'tfbs.id'])
  meth.tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[meth.tfbs.ids]$TF)
  cat(paste0('tfs overlapping with cpgs: ', paste(meth.tfbs.genes, collapse=', ') ), fill = T)
  expr.ids <- as.character(eqtl.pairs[eqtl.pairs$snps == snp, 'gene'])
  cat(paste0('expr probes associated with snp: ', paste(expr.ids, sep=', ') ), fill = T)
  expr.values <- expr.residuals[, expr.ids, drop = F]
  expr.tfbs.ids <- unique(expr.tfbs.overlap$pairs[expr.tfbs.overlap$pairs$expr.id %in% expr.ids, 'tfbs.id'])
  expr.tfbs.genes <- unique(expr.tfbs.overlap$tfbs.ranges[expr.tfbs.ids]$TF)
  cat(paste0('tfs overlapping with expr probes: ', paste(expr.tfbs.genes, collapse=', ') ), fill = T)
  
  tfbs.genes <- unique(c(snp.tfbs.genes, meth.tfbs.genes, expr.tfbs.genes))
  
  tfbs.expr.values <- lapply(tfbs.genes, function(gene) {
    probe.ids <- unique(names(genes[genes == gene]))
    if (length(probe.ids) == 1) {
      expr <- expr.residuals[, probe.ids[1]]
    } else {
      expr <- apply(expr.residuals[, probe.ids], 1, function(x) mean(x))
    }
    return(expr)
  } )
  tfbs.expr.values <- matrix(unlist(tfbs.expr.values), byrow=FALSE, ncol = length(tfbs.expr.values))
  colnames(tfbs.expr.values) <- tfbs.genes
  rownames(tfbs.expr.values) <- rownames(expr.residuals)
  
  meth.values <- get.residuals(meth.matrix, 'meth', meth.ids)
  
  snp.values <- get.snp.data(snp.ranges[snp])
  
  data.matrix <- cbind.data.frame(snp.values[id.map$axio_s4f4, , drop=F], meth.values[id.map$meth_f4, , drop=F], expr.values[id.map$expr_s4f4ogtt, , drop=F], tfbs.expr.values[id.map$expr_s4f4ogtt,])
  rownames(data.matrix) <- id.map$expr_s4f4ogtt
  
  save(data.matrix, file = paste0(PATHS$DATA.DIR, 'ggm/', set, '/', snp, '.RData'))
}
