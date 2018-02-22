

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
