source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')
source('Scripts/R/util.R')

require(ggplot2)

load(PATHS$EQTM.ME.DATA)
load(PATHS$EXPR.GENE.ANNOT.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$METH.RANGES.DATA)
load(PATHS$HERV.EQTM.OVERLAP.DATA)

cis.pos.pairs <- eqtm.me$cis$ntest
trans.pos.pairs <- eqtm.me$trans$ntests

cis.pairs <- eqtm.me$cis$eqtls
cis.pairs$snps <- as.character(cis.pairs$snps)
cis.pairs$gene <- as.character(cis.pairs$gene)
cis.cpgs <- unique(cis.pairs$snps)
cis.probes <- unique(cis.pairs$gene)
cis.genes <- unique(na.omit(probe2gene[cis.probes]))

probes.distances <- distanceToNearest(expr.ranges, snp.ranges)
pos.cis.probes <- names(expr.ranges[mcols(probes.distances)$distance < 5e5])
pos.cis.genes <- unique(probe2gene[pos.cis.probes[pos.cis.probes %in% names(probe2gene)]])
cis.gene.enrichment <- go.enrichment(cis.genes, pos.cis.genes, gsc, c('BP'))

trans.pairs <- eqtm.me$trans$eqtls
trans.pairs$snps <- as.character(trans.pairs$snps)
trans.pairs$gene <- as.character(trans.pairs$gene)
trans.snps <- unique(as.character(trans.pairs$snps))
trans.probes <- unique(as.character(trans.pairs$gene))
trans.genes <- unique(na.omit(probe2gene[trans.probes]))

hervS2.cis.either.pairs <- hervS2.eqtm.overlap$cis.either
hervS2.cis.either.cpgs <- as.character(unique(hervS2.cis.either.pairs$snps))
hervS2.cis.either.probes <- as.character(unique(hervS2.cis.either.pairs$gene))
hervS2.cis.either.genes <- unique(na.omit(probe2gene[hervS2.cis.either.probes]))


hervS2.trans.either.pairs <- hervS2.eqtm.overlap$trans.either
hervS2.trans.either.cpgs <- as.character(unique(hervS2.trans.either.pairs$snps))
hervS2.trans.either.probes <- as.character(unique(hervS2.trans.either.pairs$gene))
hervS2.trans.either.genes <- unique(na.omit(probe2gene[hervS2.trans.either.probes]))

hervS2.either.genes <- unique(c(hervS2.cis.either.genes, hervS2.trans.either.genes))
