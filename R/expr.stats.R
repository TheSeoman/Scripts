source('Scripts/R/paths.R')

require(GenomicRanges)
require(illuminaHumanv3.db)
require(reshape2)
require(grid)
require(gridExtra)

load(PATHS$EXPR.DATA)
load(PATHS$EXPR.RANGES)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)

expr.data <- f4.norm[names(expr.ranges), ]

expr.flat <- melt(expr.data)
g <- ggplot(expr.flat, aes(value)) + geom_histogram(breaks = seq(5, 14, by = 0.1)) + xlab("expression value")
g <- g + theme(text = element_text(size=14))

expr.residuals.flat <- melt(expr.residuals)
expr.residuals.flat[expr.residuals.flat$value <= 1 & expr.esiduals.flat$value >= -1, ]
g2 <- ggplot(expr.residuals.flat, aes(value)) + geom_histogram(breaks = seq(-1.5, 1.5, by = 0.05)) + xlab("expression residual")

pdf(file = paste0(PATHS$PLOT.DIR, 'expr_raw_res.pdf'), width = 7, height = 4)
expr.raw.res.dist <- grid.arrange(g, g2, ncol = 2)
dev.off()

hervS2.expr.data <- expr.data[names(hervS2.expr.overlap$expr.ranges),]
probe2gene <- unlist(as.list(illuminaHumanv3SYMBOL))
probe2gene <- probe2gene[names(expr.ranges)]
probe2gene <- probe2gene[!is.na(probe2gene)]

plot.coefficient.variance(data, type) {
  mean <- apply(data, 1, mean)
  sd <- apply(data, 1, sd)
  var <- apply(data, 1, var)
  
  
}



expr.mean <- apply(hervS2.expr.data, 1, mean)
expr.sd <- apply(hervS2.expr.data, 1, sd)
expr.var <- apply(hervS2.expr.data, 1, var)

pdf(paste0(PATHS$PLOT.DIR, 'expr_var_pdf'), width = 7, height = 4)
layout(matrix(c(1:2), 1, 2, byrow = T))
plot(expr.mean, expr.sd/expr.mean * 100, pch = 3, xlab = 'mean expression', ylab = 'coefficent of variation')
hist(expr.sd/expr.mean * 100, breaks = seq(0, 35, 0.25), xlab = 'coefficent of variation', main = '')
dev.off()



load(PATHS$METH.RANGES)

load()
