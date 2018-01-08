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
  out$either <- meqtl.pairs[meqtl.pairs$cpg %in% names(meth.ranges) | meqtl.pairs$snp %in% names(snp.ranges),]
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

get.meqtl.overlap.go.enrichment <- function (herv.meqtl.overlap, meth.annotation) {
  universe <- unique(meqtl.gene.annotation$nearestGeneSymbol)
  out <- list()
  snp.meqtl.genes <- unique(meqtl.gene.annotation[herv.meqtl.overlap$snp$cpg, 'nearestGeneSymbol'])
  snp.meqtl.genes <- snp.meqtl.genes[!is.na(snp.meqtl.genes)]
  out$snp <- go.enrichment(snp.meqtl.genes, universe, gsc, c('BP'))
  
  meth.meqtl.genes <- unique(meqtl.gene.annotation[herv.meqtl.overlap$meth$cpg, 'nearestGeneSymbol'])
  meth.meqtl.genes <- meth.meqtl.genes[!is.na(meth.meqtl.genes)]
  out$meth <- go.enrichment(meth.meqtl.genes, universe, gsc, c('BP'))
  return(out)    
}

get.meqtl.overlap.chromhmm.annotation <- function (meqtl.overlap, meth.annotation, snp.annotation) {
  #annotation of meqtl-snps, where the snps themselves lie in herv elements
  out$snp.snp <- snp.annotation[unique(meqtl.overlap$snp$snp),]
  #annotation of meqtl-cpgs, where the associated snps lie in herv elements
  out$snp.cpg <- meth.annotation[unique(meqtl.overlap$snp$cpg),]
  out$cpg.snp <- snp.annotation[unique(meqtl.overlap$meth$snp),]
  out$cpg.cpg <- meth.annotation[unique(meqtl.overlap$meth$cpg),]
  out$both.snp <- snp.annotation[unique(meqtl.overlap$both$snp),]
  out$both.cpg <- meth.annotation[unique(meqtl.overlap$both$cpg),]
  out$either.snp <- snp.annotation[unique(meqtl.overlap$either$snp),]
  out$either.cpg <- meth.annotation[unique(meqtl.overlap$either$cpg),]
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

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    overlap.name <- paste0('herv', set, flanking, '.meqtl.overlap')
    assign(overlap.name, find.meqtl.overlap(get(paste0('meth.', set, '.overlap'))$essay.ranges, get(paste0('herv', set, '.snp.ranges')), meqtl.pairs))
    assign(paste0('herv', set, flanking, '.meqtl.enrichment'), get.meqtl.overlap.go.enrichment(get(overlap.name), meqtl.gene.annotation))
    assign(paste0('herv', set, flanking, '.meqtl.annotation', extract.meqtl.annotation(get(overlap.name), meth.chromhmm.states, snp.chromhmm.states)))
  }
}

save(hervS1.meqtl.overlap, hervS2.meqtl.overlap, hervS3.meqtl.overlap, hervS1.1kb.meqtl.overlap, hervS2.1kb.meqtl.overlap, hervS3.1kb.meqtl.overlap, 
     hervS1.2kb.meqtl.overlap, hervS2.2kb.meqtl.overlap, hervS3.2kb.meqtl.overlap, file = PATHS$HERV.MEQTL.OVERLAP.DATA)

save(hervS1.meqtl.annotation, hervS1.1kb.meqtl.annotation, hervS1.2kb.meqtl.annotation, hervS2.meqtl.annotation, hervS2.1kb.meqtl.annotation, 
     hervS2.2kb.meqtl.annotation, hervS3.meqtl.annotation, hervS3.1kb.meqtl.annotation, hervS3.2kb.meqtl.annotation, file = PATHS$HERV.MEQTL.CHROMHMM.ANNOTATION.DATA)

save(hervS1.meqtl.enrichment, hervS2.meqtl.enrichment, hervS3.meqtl.enrichment, hervS1.1kb.meqtl.enrichment, hervS2.1kb.meqtl.enrichment, hervS3.1kb.meqtl.enrichment, 
     hervS1.2kb.meqtl.enrichment, hervS2.2kb.meqtl.enrichment, hervS3.2kb.meqtl.enrichment, file = PATHS$HERV.MEQTL.ENRICHMENT.DATA)
