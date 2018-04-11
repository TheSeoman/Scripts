source('Scripts/R/paths.R')

require(MatrixEQTL)

load(PATHS$SNP.RANGES.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$EXPR.GENE.ANNOT.DATA)

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/hervS2.meqtl.snp.250kb.string/')
snp <- 'rs11079148'

load(paste0(GGM.DIR, 'data.meta.RData'))


calc.ggm.eqtls <- function(GGM.DIR, id, type) {
  load(paste0(GGM.DIR, 'data/', id, '.RData'))
  load(paste0(GGM.DIR, 'data.meta.RData'))
  
  meta <- data.meta[[id]]
  
  snp.data <- SlicedData$new()
  if(type == 'snp') {
    snp.data$CreateFromMatrix(t(ggm.data[, id, drop = F]))
    snp.range <- snp.ranges[id]
    snp.pos <- data.frame(snp = id, chr = seqnames(snp.range), pos = start(snp.range))
    expr.ids <- unique(c(meta$snp.genes, meta$snp.no.gene.probes, meta$meth.genes, meta$meth.no.gene.probes, meta$tfbs.genes, meta$path.genes))
  } else {
    snps <- meta$best.snps
    snp.range <- snp.ranges[snps]
    snp.pos <- data.frame(snp = snps, chr = seqnames(snp.range), pos = start(snp.range))
    expr.ids <- unique(c(meta$snp.genes, meta$snp.no.gene.probes, meta$seed.meth.genes, meta$seed.meth.no.gene.probes, meta$extra.meth.genes, meta$extra.meth.no.gene.probes, meta$tfbs.genes, meta$path.genes))
  }
  
  expr.data <- SlicedData$new()
  expr.data$CreateFromMatrix(t(ggm.data[, expr.ids, drop = F]))
  
  expr.pos <- data.frame(matrix(sapply(expr.ids, function(id) {
    if(id %in% names(expr.ranges)) {
      range <- expr.ranges[id]
      return(c(id, as.character(seqnames(range)), start(range), end(range)))
    } else {
      probe.ids <- names(probe2gene[probe2gene == id])
      ranges <- expr.ranges[probe.ids]
      return(c(id, as.character(seqnames(ranges)[1]), min(c(start(ranges), end(ranges))), max(c(start(ranges), end(ranges)))))
    }
  }), ncol = 4, byrow = T), stringsAsFactors = F)
  colnames(expr.pos) <- c('geneid', 'chr', 's1', 's2')
  expr.pos$s1 <- as.numeric(expr.pos$s1)
  expr.pos$s2 <- as.numeric(expr.pos$s2)
  
  me = Matrix_eQTL_main(
    snps = snp.data,
    gene = expr.data,
    cvrt = SlicedData$new(),
    output_file_name = 'temp.tsv',
    pvOutputThreshold = 1,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = 'temp2.tsv',
    pvOutputThreshold.cis = 1,
    snpspos = snp.pos,
    genepos = expr.pos,
    cisDist = 5e5,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  save(me, file = paste0(GGM.DIR, 'eqtl/', id, '.RData'))
  return(me)
}

rs11079148.me <- calc.ggm.eqtls(GGM.DIR, 'rs11079148', 'snp')
rs11079148.ILMN_1902558.boxplot <- plot.single.eqtl.from.ggm.data(load.ggm.data(GGM.DIR, 'rs11079148'), 'rs11079148', 'ILMN_1902558', 'A')
rs11079148.SNORD101.boxplot <- plot.single.eqtl.from.ggm.data(load.ggm.data(GGM.DIR, 'rs11079148'), 'rs11079148', 'SNORD101', 'B')


rs1833977.me <- calc.ggm.eqtls(GGM.DIR, 'rs1833977', 'snp')
rs1833977.LDHD.boxplot <- plot.single.eqtl.from.ggm.data(load.ggm.data(GGM.DIR, 'rs1833977'), 'rs1833977', 'LDHD', 'B')
rs1833977.LDHD.boxplot <- plot.single.eqtl.from.ggm.data(load.ggm.data(GGM.DIR, 'rs1833977'), 'rs1833977', 'LDHD', 'B')


load.ggm.data <- function(GGM.DIR, id) {
  load(paste0(GGM.DIR, 'data/', id, '.RData'))
  return(ggm.data)
}

plot.single.eqtl.from.ggm.data <- function(ggm.data, snp.id, expr.id, title) {
  snp.range <- snp.ranges[snp.id]
  snp.data <- as.factor(round(ggm.data[, snp.id]))
  levels(snp.data) <- c(paste0(snp.range$from, snp.range$from), paste0(snp.range$from, snp.range$to), paste0(snp.range$to, snp.range$to))
  expr.data <- ggm.data[, expr.id, drop = T]
  p <- ggplot(data.frame(snp.data=snp.data, expr.data=expr.data), aes(snp.data,expr.data))
  p <- p + geom_boxplot() + geom_jitter(width=0.2)
  p <- p + scale_x_discrete(name = 'Genotype') + scale_y_continuous(name = 'Expression residual') + labs(title=paste0(title, ' ', paste(snp.id, expr.id, sep = '-')))
  return(p)
}



