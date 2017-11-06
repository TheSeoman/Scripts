require(gdata)
DATA.DIR = '/media/data/Masterarbeit/data/'

F.CIS.EQTL = paste0(DATA.DIR, 'eQTL/journal.pone.0093844.s005.XLS')
F.TRANS.EQTL = paste0(DATA.DIR, 'eQTL/journal.pone.0093844.s004.XLS')

S2.SNP.DATA <- paste0(DATA.DIR, 'SNPs/hervS2.SNP.RData')
EXPR.OVERLAP.DATA <- paste0(DATA.DIR, 'overlaps/expression.RData')

load(EXPR.OVERLAP.DATA)
load(S2.SNP.DATA)


cis.eqtl <- read.xls(F.CIS.EQTL, sheet = 1, header = TRUE)
trans.eqtl <- read.xls(F.TRANS.EQTL, sheet = 1, header = TRUE)

find.herv.eqtl <- function (cis.eqtl, trans.eqtl, snp, expr) {
  out <- list()
  snp.cis.eqtl <- cis.eqtl[cis.eqtl$top.SNP.KORA.F4 %in% rownames(snp),]
  snp.trans.eqtl <- trans.eqtl[trans.eqtl$top.SNP.KORA.F4 %in% rownames(snp),]

  expr.cis.eqtl <- cis.eqtl[cis.eqtl$Probe_Id %in% rownames(expr),]
  expr.trans.eqtl <- trans.eqtl[trans.eqtl$Probe_Id %in% rownames(expr),]

  both.cis.eqtl <- cis.eqtl[cis.eqtl$top.SNP.KORA.F4 %in% rownames(snp) & cis.eqtl$Probe_Id %in% rownames(expr),]
  both.trans.eqtl <- trans.eqtl[trans.eqtl$top.SNP.KORA.F4 %in% rownames(snp) & trans.eqtl$Probe_Id %in% rownames(expr),]
  
  out$snp.cis.eqtl <- snp.cis.eqtl
  out$snp.trans.eqtl <- snp.trans.eqtl
  
  out$expr.cis.eqtl <- expr.cis.eqtl
  out$expr.trans.eqtl <- expr.trans.eqtl
  
  out$both.cis.eqtl <- both.cis.eqtl
  out$both.trans.eqtl <- both.trans.eqtl
  
  return(out)
}

export.genes <- function(herv.eqtl.list, prefix) {
  if (dim(herv.eqtl.list$snp.cis.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$snp.cis.eqtl$Gene)), file = paste0(DATA.DIR, prefix, 'snp.cis.genes.txt'))
  }
  if (dim(herv.eqtl.list$snp.trans.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$snp.trans.eqtl$Gene)), file = paste0(DATA.DIR, prefix, 'snp.trans.genes.txt'))
  }
  if (dim(herv.eqtl.list$expr.cis.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$expr.cis.eqtl$Gene)), file = paste0(DATA.DIR, prefix, 'expr.cis.genes.txt'))
  }
  if (dim(herv.eqtl.list$expr.trans.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$expr.trans.eqtl$Gene)), file = paste0(DATA.DIR, prefix, 'expr.trans.genes.txt'))
  }
  if (dim(herv.eqtl.list$both.cis.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$both.cis.eqtl$Gene)), file = paste0(DATA.DIR, prefix, 'both.cis.genes.txt'))
  }
  if (dim(herv.eqtl.list$both.trans.eqtl)[1] != 0) {
    write(as.character(unique(herv.eqtl.list$both.trans.eqtl$Gene)), file = paste0(DATA.DIR, prefix, 'both.trans.genes.txt'))
  }
}

S2 <- find.herv.eqtl(cis.eqtl, trans.eqtl, hervS2.SNPs$snpInfo, expr.S2.overlap$essay.data)

export.genes(S2, 'eQTL/S2.schramm.')






