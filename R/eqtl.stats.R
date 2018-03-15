source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')
source('Scripts/R/util.R')

require(ggplot2)
require(gridExtra)

load(PATHS$MAF001.RES.ME.DATA)
load(PATHS$EXPR.GENE.ANNOT.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$METH.RESIDUALS.DATA)

load(PATHS$SNP.SAMPLES.DATA)

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

covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)
eqtl.id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                         & covariates.all$axio_s4f4 %in% snp.samples,
                         c('expr_s4f4ogtt', 'axio_s4f4')] 
eqtl.id.map <- eqtl.id.map[order(eqtl.id.map$expr_s4f4ogtt),]
eqtl.id.map$expr_s4f4ogtt <- as.character(eqtl.id.map$expr_s4f4ogtt)
eqtl.id.map$axio_s4f4 <- as.character(eqtl.id.map$axio_s4f4)

eqtm.id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                              & covariates.all$meth_f4 %in% rownames(meth.residuals),
                              c('expr_s4f4ogtt', 'meth_f4')] 
eqtm.id.map <- eqtm.id.map[order(eqtm.id.map$expr_s4f4ogtt),]
eqtm.id.map$expr_s4f4ogtt <- as.character(eqtm.id.map$expr_s4f4ogtt)
eqtm.id.map$meth_f4 <- as.character(eqtm.id.map$meth_f4)

get.snp.data <- function(snp.range, snp.samples) {
  require(Rsamtools)
  data = scanTabix(PATHS$F.SNP, param=snp.range)
  snp.data.list <- lapply(data, function(x) strsplit(x, '\t'))
  snp.data.table <- data.frame(matrix(unlist(snp.data.list), nrow=length(snp.samples)+5, byrow=F), stringsAsFactors = FALSE)
  colnames(snp.data.table) <- snp.data.table[2, ]
  snp.data.table <- snp.data.table[-(1:5), names(snp.range), drop = FALSE]
  snp.data.table <- data.frame(data.matrix(snp.data.table))
  rownames(snp.data.table) <- snp.samples
  return(snp.data.table)
}

plot.single.eqtl <- function(eqtl.pair, title) {
  snp.id <- as.character(eqtl.pair$snps)
  expr.id <- as.character(eqtl.pair$gene)
  snp.range <- snp.ranges[snp.id]
  cat(paste0('Generating boxplot for ', snp.id, '-', expr.id), fill=T)
  snp.data <- get.snp.data(snp.range, snp.samples)[eqtl.id.map$axio_s4f4, ]
  snp.data <- as.factor(round(snp.data))
  levels(snp.data) <- c(paste0(snp.range$from, snp.range$from), paste0(snp.range$from, snp.range$to), paste0(snp.range$to, snp.range$to))
  expr.data <- expr.residuals[eqtl.id.map$expr_s4f4ogtt, expr.id]
  p <- ggplot(data.frame(snp.data=snp.data, expr.data=expr.data), aes(snp.data,expr.data))
  p <- p + geom_boxplot() + geom_jitter(width=0.2)
  p <- p + scale_x_discrete(name = 'Genotype') + scale_y_continuous(name = 'Expression residual') + labs(title=paste0(title, ' ', paste(snp.id, expr.id, sep = '-')))
  return(p)
}

plot.single.eqtm <- function(eqtm.pair, title) {
  meth.id <- as.character(eqtm.pair$snps)
  expr.id <- as.character(eqtm.pair$gene)
  meth.data <- meth.residuals[eqtm.id.map$meth_f4, meth.id]
  expr.data <- expr.residuals[eqtm.id.map$expr_s4f4ogtt, expr.id]
  p <- ggplot(data.frame(meth.data=meth.data, expr.data=expr.data), aes(meth.data, expr.data)) + geom_point() 
  p <- p + labs(title=paste0(title, ' ', paste(meth.id, expr.id, sep = '-')), x='Methylation residual', y='Expression residual')
  return(p)
}

eqtl.pair <- eqtl.me$cis$eqtls[1, ]

best.cis.box <- plot.single.eqtl(cis.pairs[1,], 'A')
worst.cis.box <- plot.single.eqtl(cis.pairs[dim(cis.pairs)[1],], 'B')

best.trans.box <- plot.single.eqtl(trans.pairs[1,], 'C')
worst.trans.box <- plot.single.eqtl(trans.pairs[dim(trans.pairs)[1]-5,], 'D')

pdf(file = paste0(PATHS$PLOT.DIR, 'best_worst_eqtl_boxplots.pdf'), width = 6.3, height = 4.5)
grid.arrange(best.cis.box, worst.cis.box, best.trans.box, worst.trans.box, ncol = 2, widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()

best.eqtm.cis.scatter <- plot.single.eqtm(eqtm.me$cis$eqtls[1,], 'A')
worst.eqtm.cis.scatter <- plot.single.eqtm(eqtm.me$cis$eqtls[dim(eqtm.me$cis$eqtls)[1],], 'B')

best.eqtm.trans.scatter <- plot.single.eqtm(eqtm.me$trans$eqtls[1,], 'C')
worst.eqtm.trans.scatter <- plot.single.eqtm(eqtm.me$trans$eqtls[dim(eqtm.me$trans$eqtls)[1],], 'D')

pdf(file = paste0(PATHS$PLOT.DIR, 'best_worst_eqtm_scatter.pdf'), width = 6.3, height = 4.5)
grid.arrange(best.eqtm.cis.scatter, worst.eqtm.cis.scatter, best.eqtm.trans.scatter, worst.eqtm.trans.scatter, ncol = 2, widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()
