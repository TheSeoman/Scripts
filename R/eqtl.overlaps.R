source("Scripts/R/paths.R")
source('Scripts/R/go.enrichment.R')

load(PATHS$EXPR.OVERLAP.DATA)
load(PATHS$HERV.SNP.INFO.DATA)

find.herv.eqtl <- function (cis.eqtl, trans.eqtl, snp, expr) {
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

export.genes <- function(herv.eqtl.list, prefix) {
  if (dim(herv.eqtl.list$snp.cis.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$snp.cis.eqtl$gene[!is.na(herv.eqtl.list$snp.cis.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'snp.cis.genes.txt'))
  }
  if (dim(herv.eqtl.list$snp.trans.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$snp.trans.eqtl$gene[!is.na(herv.eqtl.list$snp.trans.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'snp.trans.genes.txt'))
  }
  if (dim(herv.eqtl.list$expr.cis.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$expr.cis.eqtl$gene[!is.na(herv.eqtl.list$expr.cis.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'expr.cis.genes.txt'))
  }
  if (dim(herv.eqtl.list$expr.trans.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$expr.trans.eqtl$gene[!is.na(herv.eqtl.list$expr.trans.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'expr.trans.genes.txt'))
  }
  if (dim(herv.eqtl.list$both.cis.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$both.cis.eqtl$gene[!is.na(herv.eqtl.list$both.cis.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'both.cis.genes.txt'))
  }
  if (dim(herv.eqtl.list$both.trans.eqtl)[1] != 0) {
    write(as.character(unique(hervsave(hervS1..eqtl.list$both.trans.eqtl$gene[!is.na(herv.eqtl.list$both.trans.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'both.trans.genes.txt'))
  }
}

eqtl.overlap.go.enrichment <- function (herv.eqtl.overlap, expr.annotation) {
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
  for (flanking in c('.', '.1kb.', '.2kb.')) {
    eqtl.overlap <- paste0('herv', set, flanking, 'eqtl.overlap')
    assign(eqtl.overlap, find.herv.eqtl(me$cis$eqtls, me$trans$eqtls, get(paste0('herv', set, flanking, 'snp.info')), get(paste0('expr.', set, flanking, 'overlap'))$essay.data))
    eqtl.enrichment <- paste0('herv', set, flanking, 'eqtl.enrichment')
    assign(eqtl.enrichment, eqtl.overlap.go.enrichment(get(eqtl.overlap), eqtl.genes))
    saveRDS(get(eqtl.enrichment), file = paste0(PATHS$DATA.DIR, 'eQTL/go.enrichments/herv', set, flanking, 'overlap.enrichment.rds'))
  }
}

save(hervS1.eqtl.overlap, hervS2.eqtl.overlap, hervS3.eqtl.overlap, hervS1.1kb.eqtl.overlap, hervS2.1kb.eqtl.overlap, hervS3.1kb.eqtl.overlap, 
     hervS1.2kb.eqtl.overlap, hervS2.2kb.eqtl.overlap, hervS3.2kb.eqtl.overlap, file = PATHS$HERV.EQTL.OVERLAP.DATA)

save(hervS1.eqtl.enrichment, hervS2.eqtl.enrichment, hervS3.eqtl.enrichment, hervS1.1kb.eqtl.enrichment, hervS2.1kb.eqtl.enrichment, hervS3.1kb.eqtl.enrichment, 
     hervS1.2kb.eqtl.enrichment, hervS2.2kb.eqtl.enrichment, hervS3.2kb.eqtl.enrichment, file = PATHS$HERV.EQTL.ENRICHMENT.DATA)

