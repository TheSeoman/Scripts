load(PATHS$SNP.RANGES.DATA)
load(PATHS$EXPR.RANGES.DATA)

load(PATHS$MAF001.ME)

cis.distances <- apply(me$cis$eqtls, 1, function (eqtl) {
  snp.id <- as.character(eqtl[1])
  expr.id <- as.character(eqtl[2])
  if (cis.snp.pos[snp.id, 1] != cis.expr.pos[expr.id, 1]) {
    distance <- -1
  } else {
    distance <- min(abs(cis.snp.pos[snp.id, 2] - cis.expr.pos[expr.id, 2]),
                    abs(cis.snp.pos[snp.id, 2] - cis.expr.pos[expr.id, 3]))
  }
  return(distance)
})

true.cis.eqtl <- me$cis$eqtls[cis.distances <= 5e5,]

cis.snp.ids <- unique(as.character(me$cis$eqtls$snps))
cis.snp.ranges <- snp.ranges[cis.snp.ids]
cis.snp.pos <- cbind.data.frame(as.character(seqnames(cis.snp.ranges)), as.numeric(start(cis.snp.ranges)))
rownames(cis.snp.pos) <- cis.snp.ids

cis.expr.ids <- unique(as.character(me$cis$eqtls$gene))
cis.expr.ranges <- expr.ranges[cis.expr.ids]
cis.expr.pos <- cbind.data.frame(as.character(seqnames(cis.expr.ranges)), as.numeric(start(cis.expr.ranges)), as.numeric(end(cis.expr.ranges)))
rownames(cis.expr.pos) <- cis.expr.ids

trans.snp.ids <- unique(as.character(me$trans$eqtls$snps))
trans.snp.ranges <- snp.ranges[trans.snp.ids]
trans.snp.pos <- cbind.data.frame(as.character(seqnames(trans.snp.ranges)), as.numeric(start(trans.snp.ranges)))
rownames(trans.snp.pos) <- trans.snp.ids

trans.expr.ids <- intersect(unique(as.character(me$trans$eqtls$gene)), names(expr.ranges))
trans.expr.ranges <- expr.ranges[trans.expr.ids]
trans.expr.pos <- cbind.data.frame(as.character(seqnames(trans.expr.ranges)), as.numeric(start(trans.expr.ranges)), as.numeric(end(trans.expr.ranges)))
rownames(trans.expr.pos) <- trans.expr.ids

trans.distances <- apply(me$trans$eqtls, 1, function (eqtl) {
  snp.id <- as.character(eqtl[1])
  expr.id <- as.character(eqtl[2])
  if (as.character(trans.snp.pos[snp.id, 1]) != as.character(trans.expr.pos[expr.id, 1]) | !expr.id %in% rownames(trans.expr.pos)) {
    distance <- -1
  } else {
    distance <- min(abs(trans.snp.pos[snp.id, 2] - trans.expr.pos[expr.id, 2]),
                    abs(trans.snp.pos[snp.id, 2] - trans.expr.pos[expr.id, 3]))
  }
  return(distance)
})

true.trans.eqtl <- me$trans$eqtls[1:1000,][trans.distances == -1 | trans.distances > 5e5,]
