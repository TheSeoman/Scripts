source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')

library(FDb.InfiniumMethylation.hg19)

cat('Loading methylation-herv overlaps...', fill = TRUE)
load(PATHS$HERV.METH.OVERLAP.DATA)
cat('Loading herv-snp overlap...', fill = TRUE)
load(PATHS$HERV.SNP.OVERLAP.DATA)
cat('Loading methylation chromHMM annotation...', fill = TRUE)
load(PATHS$METH.CHROMHMM.DATA)
cat('Loading snp chromHMM annotation...', fill = TRUE)
load(PATHS$SNP.CHROMHMM.DATA)


if(!file.exists(PATHS$MEQTL.TRANS.PAIRS.DATA)) {
  load(PATHS$MEQTL.TRANS.DATA)
  meqtl.trans.pairs <- cosmo[, c('cpg', 'snp', 'p.disco', 'p.comb')]
  rm(cosmo)
  save(meqtl.trans.pairs, file = PATHS$MEQTL.TRANS.PAIRS.DATA)
} else {
  load(PATHS$MEQTL.TRANS.PAIRS.DATA)
}

if(!file.exists(PATHS$MEQTL.PAIRS.DATA)){
  cat('Generating cosmo-pairs file...', fill = TRUE)
  load(PATHS$MEQTL.COSMO.DATA)
  meqtl.pairs <- cosmo[, c('cpg', 'snp', 'p.disco', 'p.comb')]
  rm(cosmo)
  save(meqtl.pairs, file = PATHS$MEQTL.PAIRS.DATA)
} else {
  cat('Loading cosmo-pairs file...', fill = TRUE)
  load(PATHS$MEQTL.PAIRS.DATA)
}

get.herv.meqtl.overlap <- function(meqtl.pairs, snp.ids, meth.ids) {
  out <- list()
  out$snp <- meqtl.pairs[meqtl.pairs$snp %in% snp.ids,]
  out$meth <- meqtl.pairs[meqtl.pairs$cpg %in% meth.ids,]
  out$both <- meqtl.pairs[meqtl.pairs$cpg %in% meth.ids & meqtl.pairs$snp %in% snp.ids,]
  out$either <- meqtl.pairs[meqtl.pairs$cpg %in% meth.ids | meqtl.pairs$snp %in% snp.ids,]
  return(out)
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
  
  both.meqtl.genes <- unique(meqtl.gene.annotation[herv.meqtl.overlap$both$cpg, 'nearestGeneSymbol'])
  both.meqtl.genes <- both.meqtl.genes[!is.na(both.meqtl.genes)]
  out$both <- go.enrichment(both.meqtl.genes, universe, gsc, c('BP'))
  
  either.meqtl.genes <- unique(meqtl.gene.annotation[herv.meqtl.overlap$either$cpg, 'nearestGeneSymbol'])
  either.meqtl.genes <- either.meqtl.genes[!is.na(either.meqtl.genes)]
  out$either <- go.enrichment(either.meqtl.genes, universe, gsc, c('BP'))
  return(out)    
}

get.meqtl.overlap.chromhmm.annotation <- function (meqtl.overlap, meth.annotation, snp.annotation) {
  out <- list()
  #annotation of meqtl-snps, where the snps themselves lie in herv elements
  out$snp.snp <- snp.annotation[unique(as.character(meqtl.overlap$snp$snp)),]
  #annotation of meqtl-cpgs, where the associated snps lie in herv elements
  out$snp.cpg <- meth.annotation[unique(as.character(meqtl.overlap$snp$cpg)),]
  out$cpg.snp <- snp.annotation[unique(as.character(meqtl.overlap$meth$snp)),]
  out$cpg.cpg <- meth.annotation[unique(as.character(meqtl.overlap$meth$cpg)),]
  out$both.snp <- snp.annotation[unique(as.character(meqtl.overlap$both$snp)),]
  out$both.cpg <- meth.annotation[unique(as.character(meqtl.overlap$both$cpg)),]
  out$either.snp <- snp.annotation[unique(as.character(meqtl.overlap$either$snp)),]
  out$either.cpg <- meth.annotation[unique(as.character(meqtl.overlap$either$cpg)),]
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
    cat(paste0('Processing: herv', set, flanking), fill = TRUE)
    overlap.name <- paste0('herv', set, flanking, '.meqtl.overlap')
    cat('Getting overlaps...', fill = TRUE)
    assign(overlap.name, get.herv.meqtl.overlap(meqtl.pairs, names(get(paste0('herv', set, flanking, '.snp.overlap'))$snp.ranges), names(get(paste0('herv', set, flanking, '.meth.overlap'))$meth.ranges)))
    cat('Calculating enrichments...', fill = TRUE)
    assign(paste0('herv', set, flanking, '.meqtl.enrichment'), get.meqtl.overlap.go.enrichment(get(overlap.name), meqtl.gene.annotation))
    cat('Getting annotations', fill = TRUE)
    assign(paste0('herv', set, flanking, '.meqtl.annotation.part'), get.meqtl.overlap.chromhmm.annotation(get(overlap.name), meth.chromhmm.states, snp.chromhmm.states))
  }
}

cat('Saving herv-meqtl overlaps...', fill = TRUE)
save(hervS1.meqtl.overlap, hervS2.meqtl.overlap, hervS3.meqtl.overlap, hervS1.1kb.meqtl.overlap, hervS2.1kb.meqtl.overlap, hervS3.1kb.meqtl.overlap,
     hervS1.2kb.meqtl.overlap, hervS2.2kb.meqtl.overlap, hervS3.2kb.meqtl.overlap, file = PATHS$HERV.MEQTL.OVERLAP.DATA)

cat('Saving herv-meqtl enrichments...', fill = TRUE)
save(hervS1.meqtl.enrichment, hervS2.meqtl.enrichment, hervS3.meqtl.enrichment, hervS1.1kb.meqtl.enrichment, hervS2.1kb.meqtl.enrichment, hervS3.1kb.meqtl.enrichment,
     hervS1.2kb.meqtl.enrichment, hervS2.2kb.meqtl.enrichment, hervS3.2kb.meqtl.enrichment, file = PATHS$HERV.MEQTL.ENRICHMENT.DATA)

cat('Saving herv-meqtl chromHMM annotations...', fill = TRUE)
save(hervS1.meqtl.annotation, hervS1.1kb.meqtl.annotation, hervS1.2kb.meqtl.annotation, hervS2.meqtl.annotation, hervS2.1kb.meqtl.annotation,
     hervS2.2kb.meqtl.annotation, hervS3.meqtl.annotation, hervS3.1kb.meqtl.annotation, hervS3.2kb.meqtl.annotation, file = PATHS$HERV.MEQTL.CHROMHMM.ANNOTATION.DATA)
cat('Done.', fill = TRUE)