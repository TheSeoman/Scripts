source("Scripts/R/paths.R")
source('Scripts/R/go.enrichment.R')

load(PATHS$EXPR.OVERLAP.DATA)
load(PATHS$HERV.SNP.INFO.DATA)

load(PATHS$EXPR.CHROMHMM.DATA)
load(PATHS$SNP.CHROMHMM.DATA)

get.herv.eqtl.overlap <- function (cis.eqtl, trans.eqtl, snp, expr) {
  out <- list()
  out$snp.cis <- cis.eqtl[cis.eqtl$snps %in% rownames(snp),]
  out$snp.trans <- trans.eqtl[trans.eqtl$snps %in% rownames(snp),]
  
  out$expr.cis <-  cis.eqtl[cis.eqtl$gene %in% rownames(expr),]
  out$expr.trans <- trans.eqtl[trans.eqtl$gene %in% rownames(expr),]
  
  out$both.cis <- cis.eqtl[cis.eqtl$snps %in% rownames(snp) & cis.eqtl$gene %in% rownames(expr),]
  out$both.trans <-  trans.eqtl[trans.eqtl$snps %in% rownames(snp) & trans.eqtl$gene %in% rownames(expr),]
  
  out$either.cis <- cis.eqtl[cis.eqtl$snps %in% rownames(snp) | cis.eqtl$gene %in% rownames(expr),]
  out$either.trans <- trans.eqtl[trans.eqtl$snps %in% rownames(snp) | trans.eqtl$gene %in% rownames(expr),]
  
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
  return(out)    
}

get.eqtl.overlap.chromhmm.annotation <- function(eqtl.overlap, snp.annotation, expr.annotation) {
  out <- list()
  #annotation of cis-eqtl-snps, where the snps themselves lie in herv elements
  out$cis.snp.snp <- snp.annotation[eqtl.overlap$snp.cis$snps,]
  out$cis.expr.snp <- snp.annotation[eqtl.overlap$expr.cis$snps,]
  out$cis.both.snp <- snp.annotation[eqtl.overlap$both.cis$snps,]
  out$cis.either.snp <- snp.annotation[eqtl.overlap$either.cis$snps,]
  out$trans.snp.snp <- snp.annotation[eqtl.overlap$snp.trans$snps,]
  out$trans.expr.snp <- snp.annotation[eqtl.overlap$expr.trans$snps,]
  out$trans.both.snp <- snp.annotation[eqtl.overlap$both.trans$snps,]
  out$trans.either.snp <- snp.annotation[eqtl.overlap$either.trans$snps,]
  
  
  
  return(out)
}

# schramm eqtl
require(gdata)
cis.eqtl <- read.xls(PATHS$F.CIS.EQTL, sheet = 1, header = TRUE)
cis.eqtl <- cis.eqtl[,c(1, 5, 6)]
colnames(cis.eqtl) <- c('snpid', 'probeid', 'gene')
trans.eqtl <- read.xls(PATHS$F.TRANS.EQTL, sheet = 1, header = TRUE)
trans.eqtl <- trans.eqtl[,c(2,4, 5)]
colnames(trans.eqtl) <- c('snpid', 'probeid', 'gene')
S2 <- find.herv.eqtl(cis.eqtl, NULL, hervS2.snp.info, expr.S2.overlap$essay.data)
export.genes(S2, 'eQTL/S2.schramm.')


# small snps set eqtl
require(illuminaHumanv3.db)
genes <- unlist(as.list(illuminaHumanv3SYMBOL[as.character(me$all$eqtls$gene)]))
load(PATHS$HERV.SMALL.ME)
cis.eqtl <- cbind(as.character(me$all$eqtls$snps), as.character(me$all$eqtls$gene), genes)
colnames(cis.eqtl) <- c('snpid', 'probeid', 'gene')

# MAFOO1 snps set eqtl
require(illuminaHumanv3.db)
load(PATHS$HERV.MAF001.ME)
eqtl.genes <- unlist(as.list(illuminaHumanv3SYMBOL))[unique(c(as.character(me$cis$eqtls$gene), as.character(me$trans$eqtls$gene)))]
hervS1.eqtl.overlap <- find.herv.eqtl(me$cis$eqtls, me$trans$eqtls, hervS1.snp.info, expr.S1.overlap$essay.data)
hervS1.eqtl.overlap.enrichment <- eqtl.overlap.go.enrichment(hervS1.eqtl.overlap, eqtl.genes)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    overlap.name <- paste0('herv', set, flanking, '.eqtl.overlap')
    assign(overlap.name, get.herv.eqtl.overlap(me$cis$eqtls, me$trans$eqtls, get(paste0('herv', set, flanking, '.snp.info')), get(paste0('expr.', set, flanking, 'overlap'))$essay.data))
    assign(paste0('herv', set, flanking, '.eqtl.enrichment'), get.eqtl.overlap.go.enrichment(get(overlap.name), eqtl.genes))
    assign(paste0('herv', set, flanking, '.eqtl.annotation'), get.eqtl.chromhmm.annotation(get(overlap.name), snp.annotation, expr.annotation))
     
  }
}

save(hervS1.eqtl.overlap, hervS2.eqtl.overlap, hervS3.eqtl.overlap, hervS1.1kb.eqtl.overlap, hervS2.1kb.eqtl.overlap, hervS3.1kb.eqtl.overlap, 
     hervS1.2kb.eqtl.overlap, hervS2.2kb.eqtl.overlap, hervS3.2kb.eqtl.overlap, file = PATHS$HERV.EQTL.OVERLAP.DATA)

save(hervS1.eqtl.enrichment, hervS2.eqtl.enrichment, hervS3.eqtl.enrichment, hervS1.1kb.eqtl.enrichment, hervS2.1kb.eqtl.enrichment, hervS3.1kb.eqtl.enrichment, 
     hervS1.2kb.eqtl.enrichment, hervS2.2kb.eqtl.enrichment, hervS3.2kb.eqtl.enrichment, file = PATHS$HERV.EQTL.ENRICHMENT.DATA)

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
