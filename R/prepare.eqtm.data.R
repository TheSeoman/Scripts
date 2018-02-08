source('Scripts/R/paths.R')

require(GenomicRanges)
use.residuals <- T

covariates.all <-
  read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)

load(PATHS$EXPR.RANGES.DATA)
expr.pos <- data.frame(cbind(names(expr.ranges), as.character(seqnames(expr.ranges)), start(expr.ranges), end(expr.ranges)))
rownames(expr.pos) <- names(expr.ranges)
colnames(expr.pos) <- c('geneid', 'chr', 's1', 's2')

if(use.residuals) {
  load(PATHS$EXPR.RESIDUALS.DATA)
  expr.data <- t(expr.residuals)
  load(PATHS$METH.RESIDUALS.DATA)
  meth.data <- t(meth.residuals)
} else {
  load(PATHS$EXPR.DATA)
  expr.data <- f4.norm
  load(PATHS$METH.DATA)
  meth.data <- beta
}
expr.pos <- expr.pos[rownames(expr.pos) %in% rownames(expr.data),]

load(PATHS$METH.RANGES.DATA)
meth.pos <- data.frame(cbind(names(meth.ranges), as.character(seqnames(meth.ranges)), start(meth.ranges)))
rownames(meth.pos) <- names(meth.ranges)

meth.pos <- meth.pos[rownames(meth.pos) %in% rownames(meth.data),]
colnames(meth.pos) <- c('snpid', 'chr', 'pos')

id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(expr.data) &
                         covariates.all$meth_f4 %in% colnames(meth.data) , c('expr_s4f4ogtt', 'meth_f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt), ]

expr.filtered <- expr.data[rownames(expr.pos), as.character(id.map$expr_s4f4ogtt)]
meth.filtered <- meth.data[rownames(meth.pos), as.character(id.map$meth_f4)]
colnames(meth.filtered) <- colnames(expr.filtered)

if(!use.residuals) {
  covariates.filtered <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(expr.data) &
                                      covariates.all$meth_f4 %in% colnames(meth.data) , 
                                      c('expr_s4f4ogtt', 'ucsex', 'utalteru', 'utbmi', 'ul_wbc')]
  write.table(t(covariates.filtered),
              PATHS$F.EXMEQTL.COVARIATES.FILTERED,
              sep = '\t',
              quote = FALSE,
              col.names = FALSE)
}

write.table(
  expr.pos,
  PATHS$F.EQTM.EXPRESSION.POS,
  sep = '\t',
  quote = FALSE,
  row.names = FALSE)


write.table(
  expr.filtered,
  PATHS$F.EQTM.EXPRESSION.FILTERED,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE)

write.table(
  meth.pos,
  PATHS$F.EQTM.METHYLATION.POS,
  sep = '\t',
  quote = FALSE,
  row.names = FALSE)

write.table(
  meth.filtered,
  PATHS$F.EQTM.METHYLATION.FILTERED ,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE)


