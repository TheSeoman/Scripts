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

#check cis eqtl

load(PATHS$MAF001.ME)
distances <- apply(me$cis$eqtls[1:4,], 1, function(eqtl) {
  snp.id <- as.character(eqtl[1])
  expr.id <- as.character(eqgl[2])
  return(abs(snp.pos[as.character(eqtl[1]),'pos'] - expr.pos[as.character(eqtl[2]), 'start']))
})

expr.pos <- cbind.data.frame(as.character(seqnames(expr.ranges)), as.numeric(start(expr.ranges)), as.numeric(end(expr.ranges)))
rownames(expr.pos) <- names(expr.ranges)
colnames(expr.pos) <- c('chr', 'start', 'end')
