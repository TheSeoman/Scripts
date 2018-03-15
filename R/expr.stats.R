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
load(PATHS$HERV.METH.OVERLAP.DATA)

load(PATHS$HERV.SNP.OVERLAP.DATA)


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

count.bins <- function(data, breaks, margin = 1) {
  i <- 1
  col.bins <- apply(data, margin, function(x) {
    cat(paste0('Processing sample: ', i), fill = T)
    i <<- i + 1
    bins <- table(cut(x, breaks))
    return(bins)
  })
  bins <- rowSums(col.bins)
  bin.df <- data.frame(bin = breaks[-1], count = bins)
  return(bin.df)
}

expr.data <- f4.norm[names(expr.ranges), ]

expr.hist <- plot.hist(expr.data, 'Expression', 'A', seq(5, 15.5, by = 0.1))
expr.cov.scatter <- plot.cov(expr.data, 'expression', 'B')

pdf(file = paste0(PATHS$PLOT.DIR, 'expr_raw_hist_cov.pdf'), width = 6.3, height = 3)
grid.arrange(expr.hist, expr.cov.scatter, ncol = 2)
dev.off()

#min(expr.residuals) == -6.736901, max(expr.residuals) == 7.072154
expr.res.hist <- plot.hist(expr.residuals, 'Expression residual', 'A', seq(-7.2, 7.2, by = 0.1))

hervS2.expr.data <- expr.data[names(hervS2.expr.overlap$expr.ranges),]
hervS2.expr.hist <- plot.hist(hervS2.expr.data, 'Expression', 'A', seq(5, 15.5, by = 0.1))
hervS2.expr.cov.scatter <- plot.cov(hervS2.expr.data, 'expression', 'B')
pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_expr_raw_hist_cov.pdf'), width = 6.3, height = 3)
grid.arrange(hervS2.expr.hist, hervS2.expr.cov.scatter, ncol = 2)
dev.off()


pdf(file = paste0(PATHS$PLOT.DIR, 'expr_raw_res.pdf'), width = 6.3, height = 3)
grid.arrange(expr.hist, expr.res.hist, ncol = 2)
dev.off()

meth.bins.df <- count.bins(meth.data, c(-Inf, seq(0.01, 1, by = 0.01)))
save(meth.bins.df, file = paste0(PATHS$DATA.DIR, 'Methylation/raw.meth.bins.RData'))

meth.hist <- plot.bin.hist(meth.bin.df, 0.01, 'Methylation beta', 'A')
meth.var.scatter <- plot.var(meth.data, 'methylation beta', 'B')

meth.res.bin.df <- count.bins(meth.residuals, seq(-1.25, 1.25, by = 0.025))
save(meth.res.bin.df, file = paste0(PATHS$DATA.DIR, 'Methylation/meth.res.bins.RData'))
meth.res.hist <- plot.bin.hist(meth.res.bin.df, 0.025, 'Methylation residual', 'B')

pdf(file = paste0(PATHS$PLOT.DIR, 'expr_meth_res_hist.pdf'), width = 6.3, height = 3)
grid.arrange(expr.res.hist, meth.res.hist, ncol = 2)
dev.off()

pdf(file = paste0(PATHS$PLOT.DIR, 'meth_raw_hist_var.pdf'), width = 6.3, height = 3)
png(file = paste0(PATHS$PLOT.DIR, 'meth_raw_hist_var.png'), width = 6.3/0.0138, height = 3/0.0138)
grid.arrange(meth.hist, meth.var.scatter, ncol = 2)
dev.off()

hervS2.meth.data <- meth.data[names(hervS2.meth.overlap$meth.ranges),]
hervS2.meth.bin.df <- count.bins(hervS2.meth.data, c(-Inf, seq(0, 1, by = 0.01)), 2)
hervS2.meth.hist <- plot.bin.hist(hervS2.meth.bin.df, 0.01, 'Methylation beta','A')
hervS2.meth.var.scatter <- plot.var(hervS2.meth.data, 'methylation beta', 'B')
pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_meth_raw_hist_var.pdf'), width = 6.3, height = 3) 
grid.arrange(hervS2.meth.hist, hervS2.meth.var.scatter, ncol = 2)
dev.off()

na.row.count <-
meth.probe.na.prop <- data.frame(value =  apply(meth.data, 1, function(row) {return(sum(is.na(row)))})/dim(meth.data)[2])
meth.sample.na.prop <-  data.frame(value =  apply(meth.data, 2, function(col) {return(sum(is.na(col)))})/dim(meth.data)[1])
meth.probe.na.hist <-  ggplot(meth.probe.na.prop, aes(value)) + geom_histogram(breaks = seq(0, 1, by = 0.01)) + 
  theme(text = element_text(size=10)) + labs(title='A', x = 'Proportion of NAs', y = 'Count') + scale_y_log10()
meth.sample.na.hist <-  ggplot(meth.sample.na.prop, aes(value)) + geom_histogram(breaks = seq(0, 0.06, by = 0.0006)) + 
  theme(text = element_text(size=10)) + labs(title='B', x = 'Proportion of NAs', y = 'Count')

pdf(file = paste0(PATHS$PLOT.DIR, 'meth_na_hist.pdf'), width = 6.3, height = 3)
grid.arrange(meth.probe.na.hist, meth.sample.na.hist, ncol = 2)
dev.off()

hervS2.snp.count <- data.frame(table(hervS2.snp.overlap$pairs$herv.id))
hervS2.snp.hist <-  ggplot(hervS2.snp.count, aes(Freq)) + geom_histogram(breaks = seq(0, 20, by = 1)) + 
  theme(text = element_text(size=10)) + labs(title='A', x = 'SNPs per HERV', y = 'Count')
#2423 HERVs with >40 SNPs not in plot
hervS2.snp.count$width <- width(hervS2.ranges[hervS2.snp.count$Var1])
hervS2.snp.count$density <- hervS2.snp.count$Freq/hervS2.snp.count$width

hervS2.snp.density.hist <- ggplot(hervS2.snp.count, aes(density)) + geom_histogram(breaks = seq(0, 0.05, by = 0.001)) + 
  theme(text = element_text(size=10)) + labs(title='B', x = 'SNPs per bp per HERV', y = 'Count')

pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_snp_hist.pdf'), width = 6.3, height = 3)
grid.arrange(hervS2.snp.hist, hervS2.snp.density.hist, ncol = 2)
dev.off()
