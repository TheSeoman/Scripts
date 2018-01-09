source("Scripts/R/paths.R")
source('Scripts/R/go.enrichment.R')

cat('Loading expression overlap data...')
load(PATHS$EXPR.OVERLAP.DATA)
cat('Loading herv snp ranges')
load(PATHS$HERV.SNP.RANGES.DATA)

cat('Loading chromhmm data for expression probes')
load(PATHS$EXPR.CHROMHMM.DATA)
cat('Loading chromhmm data for expression snps')
load(PATHS$SNP.CHROMHMM.DATA)

cat('Loading matrix-eqtl result')
load(PATHS$MAF001.ME)


get.herv.eqtl.overlap <- function (cis.eqtl, trans.eqtl, snp.ids, expr.ids) {
  out <- list()
  out$snp.cis <- cis.eqtl[cis.eqtl$snps %in% snp.ids,]
  out$snp.trans <- trans.eqtl[trans.eqtl$snps %in% snp.ids,]
  
  out$expr.cis <-  cis.eqtl[cis.eqtl$gene %in% expr.ids,]
  out$expr.trans <- trans.eqtl[trans.eqtl$gene %in% expr.ids,]
  
  out$both.cis <- cis.eqtl[cis.eqtl$snps %in% snp.ids & cis.eqtl$gene %in% expr.ids,]
  out$both.trans <-  trans.eqtl[trans.eqtl$snps %in% snp.ids & trans.eqtl$gene %in% expr.ids,]
  
  out$either.cis <- cis.eqtl[cis.eqtl$snps %in% snp.ids | cis.eqtl$gene %in% expr.ids,]
  out$either.trans <- trans.eqtl[trans.eqtl$snps %in% snp.ids | trans.eqtl$gene %in% expr.ids,]
  
  return(out)
}

get.eqtl.overlap.go.enrichment <- function (herv.eqtl.overlap, expr.annotation) {
  universe <- unique(expr.annotation[!is.na(expr.annotation)])
  out <- list()
  snp.trans.genes <- unique(expr.annotation[herv.eqtl.overlap$snp.trans$gene])
  snp.trans.genes <- snp.trans.genes[!is.na(snp.trans.genes)]
  if (length(snp.trans.genes) > 0) {
    out$snp.trans <- go.enrichment(snp.trans.genes, universe, gsc, c('BP'))
  }
  
  expr.trans.genes <- unique(expr.annotation[herv.eqtl.overlap$expr.trans$gene])
  expr.trans.genes <- expr.trans.genes[!is.na(expr.trans.genes)]
  if (length(expr.trans.genes) > 0) {
    out$expr.trans <- go.enrichment(expr.trans.genes, universe, gsc, c('BP'))
  }
  
  both.trans.genes <- unique(expr.annotation[herv.eqtl.overlap$both.trans$gene])
  both.trans.genes <- both.trans.genes[!is.na(both.trans.genes)]
  if (length(both.trans.genes) > 0) {
    out$both.trans <- go.enrichment(both.trans.genes, universe, gsc, c('BP'))
  }
  
  either.trans.genes <- unique(expr.annotation[herv.eqtl.overlap$either.trans$gene])
  either.trans.genes <- either.trans.genes[!is.na(either.trans.genes)]
  if (length(either.trans.genes) > 0) {
    out$either.trans <- go.enrichment(either.trans.genes, universe, gsc, c('BP'))
  }
  
  snp.cis.genes <- unique(expr.annotation[herv.eqtl.overlap$snp.cis$gene])
  snp.cis.genes <- snp.cis.genes[!is.na(snp.cis.genes)]
  if (length(snp.cis.genes) > 0) {
    out$snp.cis <- go.enrichment(snp.cis.genes, universe, gsc, c('BP'))
  }
  
  expr.cis.genes <- unique(expr.annotation[herv.eqtl.overlap$expr.cis$gene])
  expr.cis.genes <- expr.cis.genes[!is.na(expr.cis.genes)]
  if (length(expr.cis.genes) > 0) {
    out$expr.cis <- go.enrichment(expr.cis.genes, universe, gsc, c('BP'))
  }
  
  both.cis.genes <- unique(expr.annotation[herv.eqtl.overlap$both.cis$gene])
  both.cis.genes <- both.cis.genes[!is.na(both.cis.genes)]
  if (length(both.cis.genes) > 0) {
    out$both.cis <- go.enrichment(both.cis.genes, universe, gsc, c('BP'))
  }
  
  either.cis.genes <- unique(expr.annotation[herv.eqtl.overlap$either.cis$gene])
  either.cis.genes <- either.cis.genes[!is.na(either.cis.genes)]
  if (length(either.cis.genes) > 0) {
    out$either.cis <- go.enrichment(either.cis.genes, universe, gsc, c('BP'))
  }
  return(out)    
}

filter.expr.chromhmm.annotation <- function(expr.annotation.list, expr.ids) {
  filtered.annotation.list <- lapply(expr.annotation.list, function(annotation) annotation[[expr.ids]])
}

get.eqtl.overlap.chromhmm.annotation <- function(eqtl.overlap, snp.annotation, expr.annotation.list) {
  out <- list()
  #annotation of cis-eqtl-snps, where the snps themselves lie in herv elements
  out$cis.snp.snp <- snp.annotation[eqtl.overlap$snp.cis$snps,]
  #annotation of cis-eqtl-snps, where the associated expression probes lie in herv elements
  out$cis.expr.snp <- snp.annotation[eqtl.overlap$expr.cis$snps,]
  out$cis.both.snp <- snp.annotation[eqtl.overlap$both.cis$snps,]
  out$cis.either.snp <- snp.annotation[eqtl.overlap$either.cis$snps,]
  out$trans.snp.snp <- snp.annotation[eqtl.overlap$snp.trans$snps,]
  out$trans.expr.snp <- snp.annotation[eqtl.overlap$expr.trans$snps,]
  out$trans.both.snp <- snp.annotation[eqtl.overlap$both.trans$snps,]
  out$trans.either.snp <- snp.annotation[eqtl.overlap$either.trans$snps,]
  
  out$cis.snp.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$snp.cis$expr)
  out$cis.expr.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$expr.cis$expr)
  out$cis.both.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$both.cis$expr)
  out$cis.either.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$either.cis$expr)
  out$trans.snp.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$snp.trans$expr)
  out$trans.expr.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$expr.trans$expr)
  out$trans.both.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$both.trans$expr)
  out$trans.either.expr <- filter.expr.chromhmm.annotation(expr.annotation.list, eqtl.overlap$either.trans$expr)
  
  return(out)
}

# MAFOO1 snps set eqtl
cat('Generating expression probes gene annotations for eqtl-probes...', fill = TRUE)
require(illuminaHumanv3.db)
eqtl.genes <- unlist(as.list(illuminaHumanv3SYMBOL))[unique(c(as.character(me$cis$eqtls$gene), as.character(me$trans$eqtls$gene)))]

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    cat(paste0('Processing: herv', set, flanking), fill = TRUE)
    overlap.name <- paste0('herv', set, flanking, '.eqtl.overlap')
    assign(overlap.name, get.herv.eqtl.overlap(me$cis$eqtls, me$trans$eqtls, names(get(paste0('herv', set, flanking, '.snp.ranges'))), get(paste0('expr.', set, flanking, '.overlap'))$essay.ranges$ids))
    assign(paste0('herv', set, flanking, '.eqtl.enrichment'), get.eqtl.overlap.go.enrichment(get(overlap.name), eqtl.genes))
    assign(paste0('herv', set, flanking, '.eqtl.annotation'), get.eqtl.chromhmm.annotation(get(overlap.name), snp.annotation, expr.annotation))
  }
}

cat('Saving results...', fill = TRUE)
save(hervS1.eqtl.overlap, hervS2.eqtl.overlap, hervS3.eqtl.overlap, hervS1.1kb.eqtl.overlap, hervS2.1kb.eqtl.overlap, hervS3.1kb.eqtl.overlap, 
     hervS1.2kb.eqtl.overlap, hervS2.2kb.eqtl.overlap, hervS3.2kb.eqtl.overlap, file = PATHS$HERV.EQTL.OVERLAP.DATA)

save(hervS1.eqtl.enrichment, hervS2.eqtl.enrichment, hervS3.eqtl.enrichment, hervS1.1kb.eqtl.enrichment, hervS2.1kb.eqtl.enrichment, hervS3.1kb.eqtl.enrichment, 
     hervS1.2kb.eqtl.enrichment, hervS2.2kb.eqtl.enrichment, hervS3.2kb.eqtl.enrichment, file = PATHS$HERV.EQTL.ENRICHMENT.DATA)

save(hervS1.eqtl.annotation, hervS2.eqtl.annotation, hervS3.eqtl.annotation, hervS1.1kb.eqtl.annotation, hervS2.1kb.eqtl.annotation, hervS3.1kb.eqtl.annotation, 
     hervS1.2kb.eqtl.annotation, hervS2.2kb.eqtl.annotation, hervS3.2kb.eqtl.annotation, file = PATHS$HERV.EQTL.ANNOTATION.DATA)

cat('Extracting significant enrichment results', fill = TRUE)
total.significant.enrichment <- data.frame(matrix(ncol = 11, nrow = 0))
for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    eqtl.enrichment <- paste0('herv', set, flanking, '.eqtl.enrichment')
    for (condition in ls(get(eqtl.enrichment))) {
      significant <- extract.significant(paste0('herv', set, flanking), condition, get(eqtl.enrichment)[[condition]])
      total.significant.enrichment <- rbind(total.significant, significant)
    }
  }
}

write.table(total.significant, file = paste0(PATHS$DATA.DIR, 'eQTL/BP.enrichment.summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
