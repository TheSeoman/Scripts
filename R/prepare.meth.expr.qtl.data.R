source('Scripts/R/paths.R')

covariates.all <-
  read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)

load(PATHS$EXPR.RANGES.DATA)
expr.pos <- data.frame(cbind(names(expr.ranges), as.character(seqnames(expr.ranges)), start(expr.ranges), end(expr.ranges)))
rownames(expr.pos) <- names(expr.ranges)
colnames(expr.pos) <- c('geneid', 'chr', 's1', 's2')

load(PATHS$EXPR.DATA)
expr.pos <- expr.pos[rownames(expr.pos) %in% rownames(f4.norm),]

load(PATHS$METH.RANGES.DATA)
meth.pos <- data.frame(cbind(names(meth.ranges), as.character(seqnames(meth.ranges)), start(meth.ranges)))
rownames(meth.pos) <- names(meth.ranges)

load(PATHS$METH.DATA)
meth.pos <- meth.pos[rownames(meth.pos) %in% rownames(beta),]

id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(f4.norm) &
                         covariates.all$meth_f4 %in% colnames(beta) , c('expr_s4f4ogtt', 'meth_f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt), ]
colnames(expr.pos) <- c('geneid', 'chr', 'pos')

expr.filtered <- f4.norm[rownames(expr.pos), as.character(id.map$expr_s4f4ogtt)]
meth.filtered <- beta[rownames(meth.pos), as.character(id.map$meth_f4)]
colnames(meth.filtered) <- colnames(expr.filtered)

covariates.filtered <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(f4.norm) &
                                      covariates.all$meth_f4 %in% colnames(beta) , 
                                      c('expr_s4f4ogtt', 'ucsex', 'utalteru', 'utbmi', 'ul_wbc')]


write.table(
  expr.pos,
  PATHS$F.EXMEQTL.EXPRESSION.POS,
  sep = '\t',
  quote = FALSE,
  row.names = FALSE)


write.table(
  expr.filtered,
  PATHS$F.EXMEQTL.EXPRESSION.FILTERED,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE)

write.table(
  meth.pos,
  PATHS$F.EXMEQTL.METHYLATION.POS,
  sep = '\t',
  quote = FALSE,
  row.names = FALSE)

write.table(
  meth.filtered,
  PATHS$F.EXMEQTL.METHYLATION.FILTERED ,
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE)

write.table(t(covariates.filtered),
            PATHS$F.EXMEQTL.COVARIATES.FILTERED,
            sep = '\t',
            quote = FALSE,
            col.names = FALSE)
