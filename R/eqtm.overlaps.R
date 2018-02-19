source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')


cat('Loading expression overlap data...', fill = TRUE)
load(PATHS$HERV.EXPR.OVERLAP.DATA)
cat('Loading methylation overlap data...', fill = TRUE)
load(PATHS$HERV.METH.OVERLAP.DATA)

# cat('Loading chromhmm data for expression probes', fill = TRUE)
# load(PATHS$EXPR.CHROMHMM.DATA)
# cat('Loading chromhmm data for expression snps', fill = TRUE)
# load(PATHS$SNP.CHROMHMM.DATA)

cat('Loading matrix-eqtl result', fill = TRUE)
load(PATHS$EQTM.ME.DATA)



get.herv.eqtm.overlap <- function(cis.eqtm, trans.eqtm, expr.ids, meth.ids) {
  out <- list()
  out$cis.expr <- cis.eqtm[cis.eqtm$gene %in% expr.ids,]
  out$trans.expr <- trans.eqtm[trans.eqtm$gene %in% expr.ids,]
  
  out$cis.meth <-  cis.eqtm[cis.eqtm$snps %in% meth.ids,]
  out$trans.meth <- trans.eqtm[trans.eqtm$snps %in% meth.ids,]
  
  out$cis.both <- cis.eqtm[cis.eqtm$gene %in% expr.ids & cis.eqtm$snps %in% meth.ids,]
  out$trans.both <-  trans.eqtm[trans.eqtm$gene %in% expr.ids & trans.eqtm$snps %in% meth.ids,]
  
  out$cis.either <- cis.eqtm[cis.eqtm$gene %in% expr.ids | cis.eqtm$snps %in% meth.ids,]
  out$trans.either <- trans.eqtm[trans.eqtm$gene %in% expr.ids | trans.eqtm$snps %in% meth.ids,]
  return(out)
}

get.eqtm.overlap.go.enrichment <- function (herv.eqtm.overlap, expr.annotation) {
  universe <- unique(expr.annotation[!is.na(expr.annotation)])
  out <- list()
  meth.trans.genes <- unique(expr.annotation[herv.eqtm.overlap$trans.meth$gene])
  meth.trans.genes <- meth.trans.genes[!is.na(meth.trans.genes)]
  if (length(meth.trans.genes) > 0) {
    out$meth.trans <- go.enrichment(meth.trans.genes, universe, gsc, c('BP'))
  }
  
  expr.trans.genes <- unique(expr.annotation[herv.eqtm.overlap$trans.expr$gene])
  expr.trans.genes <- expr.trans.genes[!is.na(expr.trans.genes)]
  if (length(expr.trans.genes) > 0) {
    out$expr.trans <- go.enrichment(expr.trans.genes, universe, gsc, c('BP'))
  }
  
  both.trans.genes <- unique(expr.annotation[herv.eqtm.overlap$trans.both$gene])
  both.trans.genes <- both.trans.genes[!is.na(both.trans.genes)]
  if (length(both.trans.genes) > 0) {
    out$both.trans <- go.enrichment(both.trans.genes, universe, gsc, c('BP'))
  }
  
  either.trans.genes <- unique(expr.annotation[herv.eqtm.overlap$trans.either$gene])
  either.trans.genes <- either.trans.genes[!is.na(either.trans.genes)]
  if (length(either.trans.genes) > 0) {
    out$either.trans <- go.enrichment(either.trans.genes, universe, gsc, c('BP'))
  }
  
  meth.cis.genes <- unique(expr.annotation[herv.eqtm.overlap$cis.meth$gene])
  meth.cis.genes <- meth.cis.genes[!is.na(meth.cis.genes)]
  if (length(meth.cis.genes) > 0) {
    out$meth.cis <- go.enrichment(meth.cis.genes, universe, gsc, c('BP'))
  }
  
  expr.cis.genes <- unique(expr.annotation[herv.eqtm.overlap$cis.expr$gene])
  expr.cis.genes <- expr.cis.genes[!is.na(expr.cis.genes)]
  if (length(expr.cis.genes) > 0) {
    out$expr.cis <- go.enrichment(expr.cis.genes, universe, gsc, c('BP'))
  }
  
  both.cis.genes <- unique(expr.annotation[herv.eqtm.overlap$cis.both$gene])
  both.cis.genes <- both.cis.genes[!is.na(both.cis.genes)]
  if (length(both.cis.genes) > 0) {
    out$both.cis <- go.enrichment(both.cis.genes, universe, gsc, c('BP'))
  }
  
  either.cis.genes <- unique(expr.annotation[herv.eqtm.overlap$cis.either$gene])
  either.cis.genes <- either.cis.genes[!is.na(either.cis.genes)]
  if (length(either.cis.genes) > 0) {
    out$either.cis <- go.enrichment(either.cis.genes, universe, gsc, c('BP'))
  }
  return(out)    
}

get.eqtl.overlap.chromhmm.annotation <- function(eqtm.overlap, meth.annotation, expr.annotation) {
  out <- list()
  #annotation of cis-eqtl-methylation probes, where the probes themselves lie in herv elements
  out$cis.meth.meth <- meth.annotation[eqtl.overlap$meth.cis$snps,]
  #annotation of cis-eqtm-methylation probes, where the associated expression probes lie in herv elements
  out$cis.expr.meth <- meth.annotation[eqtl.overlap$expr.cis$snps,]
  out$cis.both.meth <- meth.annotation[eqtl.overlap$both.cis$snps,]
  out$cis.either.meth <- meth.annotation[eqtl.overlap$either.cis$snps,]
  out$trans.meth.meth <- meth.annotation[eqtl.overlap$meth.trans$snps,]
  out$trans.expr.meth <- meth.annotation[eqtl.overlap$expr.trans$snps,]
  out$trans.both.meth <- meth.annotation[eqtl.overlap$both.trans$snps,]
  out$trans.either.meth <- meth.annotation[eqtl.overlap$either.trans$snps,]
  
  #annotation of trans-eqtm-expression probes, where associated methylation probes themselves lie in herv elements
  out$cis.meth.expr <- expr.annotation[eqtl.overlap$meth.cis$snps,]
  #annotation of trans-eqtm-expression probes, where the probes themselves lie in herv elements
  out$cis.expr.expr <- expr.annotation[eqtl.overlap$expr.cis$snps,]
  out$cis.both.expr <- expr.annotation[eqtl.overlap$both.cis$snps,]
  out$cis.either.expr <- expr.annotation[eqtl.overlap$either.cis$snps,]
  out$trans.meth.expr <- expr.annotation[eqtl.overlap$meth.trans$snps,]
  out$trans.expr.expr <- expr.annotation[eqtl.overlap$expr.trans$snps,]
  out$trans.both.expr <- expr.annotation[eqtl.overlap$both.trans$snps,]
  out$trans.either.expr <- expr.annotation[eqtl.overlap$either.trans$snps,]
  
  return(out)
}

cat('Generating expression probes gene annotations for eqtl-probes...', fill = TRUE)
require(illuminaHumanv3.db)
eqtm.genes <- unlist(as.list(illuminaHumanv3SYMBOL))[unique(c(as.character(eqtm$cis$eqtls$gene), as.character(eqtm$trans$eqtls$gene)))]

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    cat(paste0('Processing: herv', set, flanking), fill = TRUE)
    overlap.name <- paste0('herv', set, flanking, '.eqtm.overlap')
    assign(overlap.name, get.herv.eqtm.overlap(eqtm.me$cis$eqtls, eqtm.me$trans$eqtls, names(get(paste0('herv', set, flanking, '.expr.overlap'))$expr.ranges), names(get(paste0('herv', set, flanking, '.meth.overlap'))$meth.ranges)))
    assign(paste0('herv', set, flanking, '.eqtm.enrichment'), get.eqtm.overlap.go.enrichment(get(overlap.name), eqtm.genes))
    assign(paste0('herv', set, flanking, '.eqtm.annotation'), get.eqtm.chromhmm.annotation(get(overlap.name), snp.chromhmm.states, expr.chromhmm.annotation))
  }
}

cat('Saving results...', fill = TRUE)
save(hervS1.eqtm.overlap, hervS2.eqtm.overlap, hervS3.eqtm.overlap, hervS1.1kb.eqtm.overlap, hervS2.1kb.eqtm.overlap, hervS3.1kb.eqtm.overlap,
     hervS1.2kb.eqtm.overlap, hervS2.2kb.eqtm.overlap, hervS3.2kb.eqtm.overlap, file = PATHS$HERV.EQTM.OVERLAP.DATA)

save(hervS1.eqtm.enrichment, hervS2.eqtm.enrichment, hervS3.eqtm.enrichment, hervS1.1kb.eqtm.enrichment, hervS2.1kb.eqtm.enrichment, hervS3.1kb.eqtm.enrichment,
     hervS1.2kb.eqtm.enrichment, hervS2.2kb.eqtm.enrichment, hervS3.2kb.eqtm.enrichment, file = PATHS$HERV.EQTM.ENRICHMENT.DATA)

save(hervS1.eqtm.annotation, hervS2.eqtm.annotation, hervS3.eqtm.annotation, hervS1.1kb.eqtm.annotation, hervS2.1kb.eqtm.annotation, hervS3.1kb.eqtm.annotation,
     hervS1.2kb.eqtm.annotation, hervS2.2kb.eqtm.annotation, hervS3.2kb.eqtm.annotation, file = PATHS$HERV.EQTM.ANNOTATION.DATA)

cat('Extracting significant enrichment results', fill = TRUE)
total.significant.enrichment <- data.frame(matrix(ncol = 11, nrow = 0))
for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    eqtm.enrichment <- paste0('herv', set, flanking, '.eqtm.enrichment')
    for (condition in ls(get(eqtm.enrichment))) {
      significant <- extract.significant(paste0('herv', set, flanking), condition, get(eqtm.enrichment)[[condition]])
      total.significant.enrichment <- rbind(total.significant.enrichment, significant)
    }
  }
}

write.table(total.significant, file = paste0(PATHS$DATA.DIR, 'eQTM/BP.enrichment.summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)


