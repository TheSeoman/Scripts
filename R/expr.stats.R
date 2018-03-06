source('Scripts/R/paths.R')

require(GenomicRanges)
require(illuminaHumanv3.db)
require(reshape2)
require(grid)
require(gridExtra)
require(ggplot2)

load(PATHS$EXPR.DATA)
load(PATHS$EXPR.RANGES)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)

load(PATHS$METH.DATA)
load(PATHS$METH.RANGES.DATA)
load(PATHS$METH.RESIDUALS.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)

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

plot.bin.hist <- function(bin.df, bar.width, xlab, title){
  hist <- ggplot(bin.df, aes(x=bin, y=count)) + geom_bar(stat="identity", width = bar.width) +
    theme(text = element_text(size=10)) + labs(title=title, x = xlab, y = 'Count') 
  return(hist)
}

count.bins <- function(data, breaks) {
  i <- 1
  col.bins <- apply(data, 2, function(x) {
    cat(paste0('Processing sample: ', i), fill = T)
    i <<- i + 1
    bins <- table(cut(x, breaks))
    return(bins)
  })
  bins <- rowSums(col.bins)
  return(bins)
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

meth.bins <- count.bins(meth.dat, c(-Inf, seq(0.01, 1, by = 0.01)))
save(meth.bins, file = paste0(PATHS$DATA.DIR, 'Methylation/raw.meth.bins.RData'))

meth.bin.df <- data.frame(bin = seq(0.005, 1, by = 0.01), count = meth.bins)

meth.hist <- plot.bin.hist(meth.bin.df, 0.01, 'Methylation beta', 'A')
meth.var.scatter <- plot.var(meth.data, 'methylation beta', 'B')

na.row.count <- apply(meth.data, 1, function(row) {return(sum(is.na(row)))})

pdf(file = paste0(PATHS$PLOT.DIR, 'meth_raw_hist_var.pdf'), width = 6.3, height = 3)
png(file = paste0(PATHS$PLOT.DIR, 'meth_raw_hist_var.png'), width = 6.3/0.0138, height = 3/0.0138)
grid.arrange(meth.hist, meth.var.scatter, ncol = 2)
dev.off()

