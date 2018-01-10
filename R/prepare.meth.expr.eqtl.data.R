source('Scripts/R/paths.R')

covariates.all <-
  read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)

load(PATHS$EXPR.RANGES.DATA)
expr.pos <- data.frame(cbind(names(expr.ranges), as.character(seqnames(expr.ranges)), start(expr.ranges), end(expr.ranges)))
rownames(expr.pos) <- names(expr.ranges)
colnames(expr.pos) <- c('geneid', 'chr', 's1', 's2')

load(PATHS$METH.RANGES.DATA)
meth.pos <- data.frame((cbind(names(meth.ranges), as.character(seqnames(meth.ranges)), start(meth.ranges), end(meth.ranges)))

load(PATHS$EXPR.DATA)
load(PATHS$METH.DATA)
id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(f4.norm) & covariates.all$meth_f4 %in% colnames(beta) , c('expr_s4f4ogtt', 'meth_f4')]
expr.ordered <- f4.norm[rownames(expr.pos), as.character(id.map$expr_s4f4ogtt)]




meth.ordered <- beta[, as.character(id.map$meth_f4)]

get.expr.pos <- function(expr.ranges) {

}