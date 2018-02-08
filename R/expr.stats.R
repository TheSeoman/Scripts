source('Scripts/R/paths.R')

require(GenomicRanges)
require(illuminaHumanv3.db)

load(PATHS$EXPR.DATA)
load(PATHS$EXPR.RANGES)

expr.data <- f4.norm[names(expr.ranges), ]

probe2gene <- unlist(as.list(illuminaHumanv3SYMBOL))
probe2gene <- probe2gene[names(expr.ranges)]
probe2gene <- probe2gene[!is.na(probe2gene)]

expr.mean <- apply(expr.data, 1, mean)
expr.sd <- apply(expr.data, 1, sd)
expr.var <- apply(expr.data, 1, var)

pdf('~/Scripts/thesis/figures/expr.var.pdf', width = 7, height = 4)
layout(matrix(c(1:2), 1, 2, byrow = T))
plot(expr.mean, expr.sd/expr.mean * 100, pch = 3, xlab = 'mean expression', ylab = 'coefficent of variation', main = 'Expression per probe')
hist(expr.sd/expr.mean * 100, breaks = seq(0, 35, 0.25), xlab = 'coefficent of variation', main = '')
dev.off()
 