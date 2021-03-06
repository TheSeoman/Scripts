source("Scripts/R/paths.R")
source('Scripts/R/go.enrichment.R')

cat('Loading expression overlap data...', fill = TRUE)
load(PATHS$HERV.EXPR.OVERLAP.DATA)
cat('Loading herv snp overlap data')
load(PATHS$HERV.SNP.OVERLAP.DATA)

# cat('Loading chromhmm data for expression probes', fill = TRUE)
# load(PATHS$EXPR.CHROMHMM.DATA)
# cat('Loading chromhmm data for snps', fill = TRUE)
# load(PATHS$SNP.CHROMHMM.DATA)

cat('Loading matrix-eqtl result', fill = TRUE)

load(PATHS$MAF001.RES.ME.DATA)
load(PATHS$EXPR.GENE.ANNOT.DATA)
# MAFOO1 snps set eqtl
cat('Generating expression probes gene annotations for eqtl-probes...', fill = TRUE)
eqtl.genes <- unique(na.omit(probe2gene[unique(c(as.character(eqtl.me$cis$eqtls$gene), as.character(eqtl.me$trans$eqtls$gene)))]))

get.eqtl.overview <- function(eqtl.table) {
  out <- list()
  out$pairs <- dim(eqtl.table)[1]
  out$snps <- length(unique(eqtl.table$snps))
  out$probes <- length(unique(eqtl.table$gene))
  genes <- unique(eqtl.genes[as.character(unique(eqtl.table$gene))])
  out$genes <- length(genes[!is.na(genes)])
  return(out)
}

get.herv.eqtl.overlap <- function (cis.eqtl, trans.eqtl, snp.ids, expr.ids) {
  out <- list()
  out$cis.snp <- cis.eqtl[cis.eqtl$snps %in% snp.ids,]
  out$trans.snp <- trans.eqtl[trans.eqtl$snps %in% snp.ids,]
  
  out$cis.expr <-  cis.eqtl[cis.eqtl$gene %in% expr.ids,]
  out$trans.expr <- trans.eqtl[trans.eqtl$gene %in% expr.ids,]
  
  out$cis.both <- cis.eqtl[cis.eqtl$snps %in% snp.ids & cis.eqtl$gene %in% expr.ids,]
  out$trans.both <-  trans.eqtl[trans.eqtl$snps %in% snp.ids & trans.eqtl$gene %in% expr.ids,]
  
  out$cis.either <- cis.eqtl[cis.eqtl$snps %in% snp.ids | cis.eqtl$gene %in% expr.ids,]
  out$trans.either <- trans.eqtl[trans.eqtl$snps %in% snp.ids | trans.eqtl$gene %in% expr.ids,]
  
  return(out)
}

get.eqtl.overlap.go.enrichment <- function (herv.eqtl.overlap, eqtl.genes) {
  out <- list()
  for(type in c('snp', 'expr', 'both', 'either')) {
    type.eqtl.genes <- unique(na.omit(probe2gene[unique(c(as.character(herv.eqtl.overlap[[paste0('cis.', type)]]$gene), 
                                                          as.character(herv.eqtl.overlap[[paste0('trans.', type)]]$gene)))]))
    cat(paste0('Calculating GO enrichment for HERV ', type, ' eQTLs on ', length(type.eqtl.genes), ' genes.'), fill = T)
    if(length(type.eqtl.genes) > 0) {
      out[[type]] <- go.enrichment(type.eqtl.genes, eqtl.genes, gsc, c('BP'))
    } else {
      out[[type]] <- NULL
    }
  }
  return(out)
}

print.eqtl.overlap.summary <- function(overlap) {
  for(key in names(overlap)) {
    cat(key, fill = TRUE) 
    e <- overlap[[key]]
    cat(paste(dim(e)[1], length(unique(e$snps)), length(unique(e$gene)), sep = '\t'), fill = TRUE)
  }
}

get.eqtl.overlap.chromhmm.annotation <- function(eqtl.overlap, snp.annotation, expr.annotation) {
  out <- list()
  #annotation of cis-eqtl-snps, where the snps themselves lie in herv elements
  out$cis.snp.snp <- snp.annotation[eqtl.overlap$cis.snp$snps,]
  #annotation of cis-eqtl-snps, where the associated expression probes lie in herv elements
  out$cis.expr.snp <- snp.annotation[eqtl.overlap$cis.expr$snps,]
  out$cis.both.snp <- snp.annotation[eqtl.overlap$cis.both$snps,]
  out$cis.either.snp <- snp.annotation[eqtl.overlap$cis.either$snps,]
  out$trans.snp.snp <- snp.annotation[eqtl.overlap$trans.snp$snps,]
  out$trans.expr.snp <- snp.annotation[eqtl.overlap$trans.expr$snps,]
  out$trans.both.snp <- snp.annotation[eqtl.overlap$trans.both$snps,]
  out$trans.either.snp <- snp.annotation[eqtl.overlap$trans.either$snps,]
  
  #annotation of cis-eqtl-snps, where the snps themselves lie in herv elements
  out$cis.snp.expr <- expr.annotation[eqtl.overlap$cis.snp$snps,]
  #annotation of cis-eqtl-snps, where the associated expression probes lie in herv elements
  out$cis.expr.expr <- expr.annotation[eqtl.overlap$cis.expr$snps,]
  out$cis.both.expr <- expr.annotation[eqtl.overlap$cis.both$snps,]
  out$cis.either.expr <- expr.annotation[eqtl.overlap$cis.either$snps,]
  out$trans.snp.expr <- expr.annotation[eqtl.overlap$trans.snp$snps,]
  out$trans.expr.expr <- expr.annotation[eqtl.overlap$trans.expr$snps,]
  out$trans.both.expr <- expr.annotation[eqtl.overlap$trans.both$snps,]
  out$trans.either.expr <- expr.annotation[eqtl.overlap$trans.either$snps,]
  
  return(out)
}

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    cat(paste0('Processing: herv', set, flanking), fill = TRUE)
    overlap.name <- paste0('herv', set, flanking, '.eqtl.overlap')
    assign(overlap.name, get.herv.eqtl.overlap(eqtl.me$cis$eqtls, eqtl.me$trans$eqtls, names(get(paste0('herv', set, flanking, '.snp.overlap'))$snp.ranges), names(get(paste0('herv', set, flanking, '.expr.overlap'))$expr.ranges)))
    assign(paste0('herv', set, flanking, '.eqtl.enrichment'), get.eqtl.overlap.go.enrichment(get(overlap.name), eqtl.genes))
    # assign(paste0('herv', set, flanking, '.eqtl.annotation'), get.eqtl.overlap.chromhmm.annotation(get(overlap.name), snp.chromhmm.states, expr.chromhmm.annotation))
  }
}

cat('Saving results...', fill = TRUE)
save(hervS1.eqtl.overlap, hervS2.eqtl.overlap, hervS3.eqtl.overlap, hervS1.1kb.eqtl.overlap, hervS2.1kb.eqtl.overlap, hervS3.1kb.eqtl.overlap,
     hervS1.2kb.eqtl.overlap, hervS2.2kb.eqtl.overlap, hervS3.2kb.eqtl.overlap, file = PATHS$HERV.EQTL.OVERLAP.DATA)

save(hervS1.eqtl.enrichment, hervS2.eqtl.enrichment, hervS3.eqtl.enrichment, hervS1.1kb.eqtl.enrichment, hervS2.1kb.eqtl.enrichment, hervS3.1kb.eqtl.enrichment,
     hervS1.2kb.eqtl.enrichment, hervS2.2kb.eqtl.enrichment, hervS3.2kb.eqtl.enrichment, file = PATHS$HERV.EQTL.ENRICHMENT.DATA)

# save(hervS1.eqtl.annotation, hervS2.eqtl.annotation, hervS3.eqtl.annotation, hervS1.1kb.eqtl.annotation, hervS2.1kb.eqtl.annotation, hervS3.1kb.eqtl.annotation,
     # hervS1.2kb.eqtl.annotation, condition, hervS3.2kb.eqtl.annotation, file = PATHS$HERV.EQTL.CHROMHMM.ANNOTATION.DATA)

cat('Extracting significant enrichment results', fill = TRUE)
total.significant.enrichment <- data.frame(matrix(ncol = 11, nrow = 0))
for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    eqtl.enrichment <- paste0('herv', set, flanking, '.eqtl.enrichment')
    for (condition in names(get(eqtl.enrichment))) {
      significant <- extract.significant(paste0('herv', set, flanking), condition, get(eqtl.enrichment)[[condition]])
      total.significant.enrichment <- rbind(total.significant.enrichment, significant)
    }
  }
}

write.table(total.significant.enrichment, file = paste0(PATHS$DATA.DIR, 'eQTL/BP.enrichment.summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
