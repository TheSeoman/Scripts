source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')
source('Scripts/R/util.R')

require(ggplot2)

load(PATHS$MAF001.RES.ME.DATA)
load(PATHS$EXPR.GENE.ANNOT.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

cis.pos.pairs <- eqtl.me$cis$ntest
trans.pos.pairs <- eqtl.me$trans$ntests

cis.pairs <- eqtl.me$cis$eqtls
cis.pairs$snps <- as.character(cis.pairs$snps)
cis.pairs$gene <- as.character(cis.pairs$gene)
cis.snps <- unique(cis.pairs$snps)
cis.probes <- unique(cis.pairs$gene)
cis.genes <- unique(na.omit(probe2gene[cis.probes]))

probes.distances <- distanceToNearest(expr.ranges, snp.ranges)
pos.cis.probes <- names(expr.ranges[mcols(probes.distances)$distance < 5e5])
pos.cis.genes <- unique(probe2gene[pos.cis.probes[pos.cis.probes %in% names(probe2gene)]])
cis.gene.enrichment <- go.enrichment(cis.genes, pos.cis.genes, gsc, c('BP'))

trans.pairs <- eqtl.me$trans$eqtls
trans.pairs$snps <- as.character(trans.pairs$snps)
trans.pairs$gene <- as.character(trans.pairs$gene)
trans.snps <- unique(as.character(trans.pairs$snps))
trans.probes <- unique(as.character(trans.pairs$gene))
trans.genes <- unique(na.omit(probe2gene[trans.probes]))

# schramm eqtl
require(gdata)
schramm.cis.eqtl <- read.xls(PATHS$F.CIS.EQTL, sheet = 1, header = TRUE)
schramm.cis.eqtl <- schramm.cis.eqtl[,c(1, 5, 6)]
colnames(schramm.cis.eqtl) <- c('snpid', 'probeid', 'gene')
for(i in 1:3) {
  cis.eqtl[,i] <- as.character(cis.eqtl[,i])
}

schramm.trans.eqtl <- read.xls(PATHS$F.TRANS.EQTL, sheet = 1, header = TRUE)
schramm.trans.eqtl <- schramm.trans.eqtl[,c(2,4, 5)]
colnames(schramm.trans.eqtl) <- c('snpid', 'probeid', 'gene')

rep.cis <- apply(schramm.cis.eqtl, 1, function(x) {
  return(any(cis.pairs$snps==x[1] & cis.pairs$gene ==x[2]))
})

save(rep.cis, file = paste0(PATHS$DATA.DIR, 'eQTL/schramm.rep.cis.RData'))

load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$EXPR.DATA)

load(PATHS$SNP.SAMPLES.DATA)
covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)
id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                         & covariates.all$axio_s4f4 %in% snp.samples,
                         c('expr_s4f4ogtt', 'axio_s4f4')] 
id.map <- id.map[order(id.map$expr_s4f4ogtt),]
id.map$expr_s4f4ogtt <- as.character(id.map$expr_s4f4ogtt)
id.map$axio_s4f4 <- as.character(id.map$axio_s4f4)

get.snp.data <- function(snp.range, snp.samples) {
  require(Rsamtools)
  data = scanTabix(PATHS$F.SNP, param=snp.range)
  snp.data.list <- lapply(data, function(x) strsplit(x, '\t'))
  snp.data.table <- data.frame(matrix(unlist(snp.data.list), nrow=length(snp.samples)+5, byrow=F), stringsAsFactors = FALSE)
  colnames(snp.data.table) <- snp.data.table[2, ]
  snp.data.table <- snp.data.table[-(1:5), names(snp.range), drop = FALSE]
  snp.data.table[, names(snp.range)] <- as.numeric(snp.data.table[, names(snp.range)])
  rownames(snp.data.table) <- snp.samples
  return(snp.data.table)
}

plot.single.eqtl <- function(eqtl.pair) {
  snp.id <- eqtl.pair$snps
  expr.id <- eqtl.pair$gene
  snp.range <- snp.ranges[snp.id]
  cat(paste0('Generating boxplot for ', snp.id, '-', expr.id), fill=T)
  snp.data <- get.snp.data(snp.range, snp.samples)[id.map$axio_s4f4, ]
  snp.data <- as.factor(round(snp.data))
  levels(snp.data) <- c(paste0(snp.range$from, snp.range$from), paste0(snp.range$from, snp.range$to), paste0(snp.range$to, snp.range$to))
  expr.data <- expr.residuals[id.map$expr_s4f4ogtt, expr.id]
  p <- ggplot(data.frame(snp.data=snp.data, expr.data=expr.data), aes(snp.data,expr.data))
  p <- p + geom_boxplot() + geom_jitter(width=0.2)
  p <- p + scale_x_discrete(name = 'Genotype') + scale_y_continuous(name = 'Expression residual')
  return(p)
}
