source('Scripts/R/paths.R')

require(GenomicRanges)
require(illuminaHumanv3.db)
require(reshape2)

load(PATHS$EXPR.DATA)
load(PATHS$EXPR.RANGES)
load(PATHS$EXPR.RESIDUALS.DATA)

expr.data <- f4.norm[names(expr.ranges), ]

expr.flat <- melt(expr.data)
g <- ggplot(expr.flat, aes(value)) + geom_histogram(breaks = seq(5, 15.4, by = 0.1)) + xlab("expression value")
g <- g + theme(text = element_text(size=14))

expr.residuals.flat <- melt(expr.residuals)
g <- ggplot(expr.residuals.flat, aes(value)) + geom_histogram(breaks = seq(-1, 1, by = 0.1)) + xlab("expression residual")

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



load(PATHS$METH.RANGES)
