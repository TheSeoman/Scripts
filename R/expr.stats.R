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

load(PATHS$METH.DATA)
load(PATHS$METH.RANGES.DATA)

plot.cov <- function(data, type, title) {
  mean <- apply(data, 1, mean, na.rm = T)
  sd <- apply(data, 1, sd, na.rm = T)
  df <- data.frame(mean = mean, coef.var = sd/mean * 100)
  g <- ggplot(df, aes(mean, coef.var)) + geom_point() + theme(text = element_text(size=10)) +
    labs(title=title, x = paste0('Mean ', type), y = 'Coefficient of variation')
  return(g)
}

plot.var <- function(data, type, title) {
  mean <- apply(data, 1, mean, na.rm = T)
  var <- apply(data, 1, var, na.rm = T)
  df <- data.frame(mean = mean, var = var)
  g <- ggplot(df, aes(mean,var)) + geom_point() + theme(text = element_text(size=10)) +
    labs(title=title, x = paste0('Mean ', type), y = 'Variance')
  return(g)
}


plot.hist <- function(data, xlab, title, breaks) {
  data.flat <- melt(data)
  hist <- ggplot(data.flat, aes(value)) + geom_histogram(breaks = breaks) + 
    theme(text = element_text(size=10)) + labs(title=title, x = xlab, y = 'Count')
  return(hist)
}

expr.data <- f4.norm[names(expr.ranges), ]

expr.hist <- plot.hist(expr.data, 'Expression', 'A', seq(5, 15.5, by = 0.1))
expr.cov.scatter <- plot.cov(expr.data, 'expression', 'B')

pdf(file = paste0(PATHS$PLOT.DIR, 'expr_raw_hist_cov.pdf'), width = 6.3, height = 3)
grid.arrange(expr.hist, expr.cov.scatter, ncol = 2)
dev.off()

expr.res.hist <- plot.hist(expr.residuals, )

hervS2.expr.data <- expr.data[names(hervS2.expr.overlap$expr.ranges),]

pdf(file = paste0(PATHS$PLOT.DIR, 'expr_raw_res.pdf'), width = 6.3, height = 3)
grid.arrange(expr.hist, expr.res.hist, ncol = 2)
dev.off()

meth.data <- beta[names(meth.ranges),]
meth.hist <- plot.hist(meth.data, 'Methylation', 'A', seq(0, 1, by = 0.01))
meth.var.scatter <- plot.var(meth.data, 'methylation', 'B')

pdf(file = paste0(PATHS$PLOT.DIR, 'meth_raw_hist_var.pdf'), width = 6.3, height = 3)
grid.arrange(meth.hist, meth.var.scatter, nocl = 2)
dev.off()
