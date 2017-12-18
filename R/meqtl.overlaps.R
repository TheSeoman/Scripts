source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')

library(FDb.InfiniumMethylation.hg19)

load(PATHS$METH.OVERLAP.DATA)
load(PATHS$HERV.SNP.RANGES.DATA)
if(!file.exists(PATHS$MEQTL.PAIRS.DATA)){
  load(PATHS$MEQTL.COSMO.DATA)
  meqtl.pairs <- cosmo[, c('cpg', 'snp', 'p.disco', 'p.comb')]
  save(meqtl.pairs, file = PATHS$MEQTL.PAIRS.DATA)
} else {
  load(PATHS$MEQTL.PAIRS.DATA)
}

if(!file.exists(PATHS$MEQTL.SNP.RANGES.DATA)) {
  load(PATHS$MEQTL.COSMO.DATA)
  first <- !duplicated(cosmo$snp)
  meqtl.snp.ranges <- GRanges(seqnames = paste0('chr', cosmo$snp.chr[first]), ranges = IRanges(start = cosmo$snp.pos[first], width = 1))
  names(meqtl.snp.ranges) <- cosmo$snp[first]
  save(meqtl.snp.ranges, file = PATHS$MEQTL.SNP.RANGES.DATA)
} 

if(!file.exists(PATHS$MEQTL.METH.RANGES.DATA)) {
  load(PATHS$MEQTL.COSMO.DATA)
  first <- !duplicated(cosmo$cpg)
  meqtl.meth.ranges <- GRanges(seqnames = paste0('chr', cosmo$cpg.chr[first]), ranges = IRanges(start = cosmo$cpg.pos[first], width = 1))
  names(meqtl.meth.ranges) <- cosmo$cpg[first]
  save(meqtl.meth.ranges, file = PATHS$MEQTL.METH.RANGES.DATA)
} 


find.meqtl.overlap <- function(meth.ranges, snp.ranges, meqtl.pairs) {
  out <- list()
  out$snp <- meqtl.pairs[meqtl.pairs$snp %in% names(snp.ranges),]
  out$meth <- meqtl.pairs[meqtl.pairs$cpg %in% names(meth.ranges),]
  out$both <- meqtl.pairs[meqtl.pairs$cpg %in% names(meth.ranges) & meqtl.pairs$snp %in% names(snp.ranges),]
  out$either <- both.meqtl <- meqtl.pairs[meqtl.pairs$cpg %in% names(meth.ranges) | meqtl.pairs$snp %in% names(snp.ranges),]
  return(out)
}

export.meth.overlap.genes <- function (herv.meqtl.overlap, meth.annotation, prefix) {
  snp.meqtl.genes <- unique(meth.annotation[unique(herv.meqtl.overlap$snp$cpg), 'nearestGeneSymbol'])
  snp.meqtl.genes <- snp.meqtl.genes[!is.na(snp.meqtl.genes)]
  if (length(snp.meqtl.genes) != 0) {
    write(snp.meqtl.genes, file = paste0(PATHS$DATA.DIR, 'meQTL/', prefix, '.snp.meqtl.genes.txt'))
  }
  meth.meqtl.genes <- unique(meth.annotation[unique(herv.meqtl.overlap$meth$cpg), 'nearestGeneSymbol'])
  meth.meqtl.genes <- meth.meqtl.genes[!is.na(meth.meqtl.genes)]
  if (length(snp.meqtl.genes) != 0) {
    write(meth.meqtl.genes, file = paste0(PATHS$DATA.DIR, 'meQTL/', prefix, '.meth.meqtl.genes.txt'))
  }
}

meqtl.overlap.go.enrichment <- function (herv.meqtl.overlap, meth.annotation) {
  universe <- unique(meqtl.gene.annotation$nearestGeneSymbol)
  out <- list()
  snp.meqtl.genes <- unique(meqtl.gene.annotation[herv.meqtl.overlap$snp$cpg, 'nearestGeneSymbol'])
  snp.meqtl.genes <- snp.meqtl.genes[!is.na(snp.meqtl.genes)]
  snp.meqtl.enrichment <- go.enrichment(snp.meqtl.genes, universe, gsc, c('BP'))
  
  meth.meqtl.genes <- unique(meqtl.gene.annotation[herv.meqtl.overlap$meth$cpg, 'nearestGeneSymbol'])
  meth.meqtl.genes <- meth.meqtl.genes[!is.na(meth.meqtl.genes)]
  meth.meqtl.enrichment <- go.enrichment(meth.meqtl.genes, universe, gsc, c('BP'))

  out$snp <- snp.meqtl.enrichment
  out$meth <- meth.meqtl.enrichment
  
  return(out)    
}

if (!file.exists(PATHS$MEQTL.GENE.ANNOTATION.DATA)) {
  meth.ranges <- getPlatform(platform='HM450', genome='hg19')
  meqtl.nearest.genes <- getNearestGene(meth.ranges[unique(meqtl.pairs$cpg)])
  meqtl.gene.annotation <- meqtl.nearest.genes[meqtl.nearest.genes$distance == 0, ]
  save(meqtl.gene.annotation, file = PATHS$MEQTL.GENE.ANNOTATION.DATA)
} else {
  load(PATHS$MEQTL.GENE.ANNOTATION.DATA)
}


hervS1.meqtl.overlap <- find.meqtl.overlap(meth.S1.overlap$essay.ranges, hervS1.snp.ranges, meqtl.pairs)
hervS2.meqtl.overlap <- find.meqtl.overlap(meth.S2.overlap$essay.ranges, hervS2.snp.ranges, meqtl.pairs)
hervS3.meqtl.overlap <- find.meqtl.overlap(meth.S3.overlap$essay.ranges, hervS3.snp.ranges, meqtl.pairs)

hervS1.1kb.meqtl.overlap <- find.meqtl.overlap(meth.S1.1kb.overlap$essay.ranges, hervS1.1kb.snp.ranges, meqtl.pairs)
hervS2.1kb.meqtl.overlap <- find.meqtl.overlap(meth.S2.1kb.overlap$essay.ranges, hervS2.1kb.snp.ranges, meqtl.pairs)
hervS3.1kb.meqtl.overlap <- find.meqtl.overlap(meth.S3.1kb.overlap$essay.ranges, hervS3.1kb.snp.ranges, meqtl.pairs)

hervS1.2kb.meqtl.overlap <- find.meqtl.overlap(meth.S1.2kb.overlap$essay.ranges, hervS1.2kb.snp.ranges, meqtl.pairs)
hervS2.2kb.meqtl.overlap <- find.meqtl.overlap(meth.S2.2kb.overlap$essay.ranges, hervS2.2kb.snp.ranges, meqtl.pairs)
hervS3.2kb.meqtl.overlap <- find.meqtl.overlap(meth.S3.2kb.overlap$essay.ranges, hervS3.2kb.snp.ranges, meqtl.pairs)

save(hervS1.meqtl.overlap, hervS2.meqtl.overlap, hervS3.meqtl.overlap, hervS1.1kb.meqtl.overlap, hervS2.1kb.meqtl.overlap, hervS3.1kb.meqtl.overlap, 
     hervS1.2kb.meqtl.overlap, hervS2.2kb.meqtl.overlap, hervS3.2kb.meqtl.overlap, file = PATHS$HERV.MEQTL.OVERLAP.DATA)

hervS1.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS1.meqtl.overlap, meqtl.gene.annotation)
hervS2.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS2.meqtl.overlap, meqtl.gene.annotation)
hervS3.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS3.meqtl.overlap, meqtl.gene.annotation)

hervS1.1kb.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS1.1kb.meqtl.overlap, meqtl.gene.annotation)
hervS2.1kb.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS2.1kb.meqtl.overlap, meqtl.gene.annotation)
hervS3.1kb.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS3.1kb.meqtl.overlap, meqtl.gene.annotation)

hervS1.2kb.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS1.2kb.meqtl.overlap, meqtl.gene.annotation)
hervS2.2kb.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS2.2kb.meqtl.overlap, meqtl.gene.annotation)
hervS3.2kb.meqtl.enrichment <- meqtl.overlap.go.enrichment(hervS3.2kb.meqtl.overlap, meqtl.gene.annotation)

save(hervS1.meqtl.enrichment, hervS2.meqtl.enrichment, hervS3.meqtl.enrichment, hervS1.1kb.meqtl.enrichment, hervS2.1kb.meqtl.enrichment, hervS3.1kb.meqtl.enrichment, 
     hervS1.2kb.meqtl.enrichment, hervS2.2kb.meqtl.enrichment, hervS3.2kb.meqtl.enrichment, file = PATHS$HERV.MEQTL.ENRICHMENT.DATA)
