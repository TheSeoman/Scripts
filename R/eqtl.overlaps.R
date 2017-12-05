source("Scripts/R/paths.R")

load(PATHS$EXPR.OVERLAP.DATA)
load(PATHS$HERV.SNP.INFO.DATA)

find.herv.eqtl <- function (cis.eqtl, trans.eqtl, snp, expr) {
  out <- list()
  snp.cis.eqtl <- cis.eqtl[cis.eqtl$snpid %in% rownames(snp),]
  snp.trans.eqtl <- trans.eqtl[trans.eqtl$snpid %in% rownames(snp),]
  
  expr.cis.eqtl <- cis.eqtl[cis.eqtl$probeid %in% rownames(expr),]
  expr.trans.eqtl <- trans.eqtl[trans.eqtl$probeid %in% rownames(expr),]
  
  both.cis.eqtl <- cis.eqtl[cis.eqtl$snpid %in% rownames(snp) & cis.eqtl$probeid %in% rownames(expr),]
  both.trans.eqtl <- trans.eqtl[trans.eqtl$snpid %in% rownames(snp) & trans.eqtl$probeid %in% rownames(expr),]
  
  either.cis.eqtl <- cis.eqtl[cis.eqtl$snpid %in% rownames(snp) | cis.eqtl$probeid %in% rownames(expr),]
  either.trans.eqtl <- trans.eqtl[trans.eqtl$snpid %in% rownames(snp) | trans.eqtl$probeid %in% rownames(expr),]
  
  out$snp.cis.eqtl <- snp.cis.eqtl
  out$snp.trans.eqtl <- snp.trans.eqtl
  
  out$expr.cis.eqtl <- expr.cis.eqtl
  out$expr.trans.eqtl <- expr.trans.eqtl
  
  out$both.cis.eqtl <- both.cis.eqtl
  out$both.trans.eqtl <- both.trans.eqtl
  
  out$either.cis.eqtl <- either.cis.eqtl
  out$either.trans.eqtl <- either.trans.eqtl
  
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
    write(as.character(unique(herv.eqtl.list$both.trans.eqtl$gene[!is.na(herv.eqtl.list$both.trans.eqtl$gene)])), file = paste0(PATHS$DATA.DIR, prefix, 'both.trans.genes.txt'))
  }
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
all.genes <- unlist(as.list(illuminaHumanv3SYMBOL))
cis.genes <- unlist(as.list(illuminaHumanv3SYMBOL[as.character(me$cis$eqtls$gene)]))
cis.eqtl <- as.data.frame(cbind(as.character(me$cis$eqtls$snps), as.character(me$cis$eqtls$gene), cis.genes))
colnames(cis.eqtl) <- c('snpid', 'probeid', 'gene')
trans.genes <- unlist(as.list(illuminaHumanv3SYMBOL[as.character(me$trans$eqtls$gene)]))
trans.eqtl <- as.data.frame(cbind(as.character(me$trans$eqtls$snps), as.character(me$trans$eqtls$gene), trans.genes))
colnames(trans.eqtl) <- c('snpid', 'probeid', 'gene')
hervS1.eqtl <- find.herv.eqtl(cis.eqtl, trans.eqtl, hervS1.snp.info, expr.S1.overlap$essay.data)
S2 <- find.herv.eqtl(cis.eqtl, trans.eqtl, hervS2.snp.info, expr.S2.overlap$essay.data)
export.genes(S2, 'eQTL/S2.MAF001.')

library(GSEABase)
library(hgu95av2.db)
library(GO.db)
