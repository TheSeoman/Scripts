source('Scripts/R/paths.R')

require(ggplot2)

load(PATHS$HERV.RANGES.DATA)

herv.stats <- data.frame(matrix(nrow = 3, ncol = 3))
rownames(herv.stats) <- c('S1', 'S2', 'S3')
colnames(herv.stats) <- c('count', 'mean.width', 'total.width')

for(set in c('S1', 'S2', 'S3')) {
  herv.ranges <- get(paste0('herv', set, '.ranges'))
  widths <- width(herv.ranges)
  herv.stats[set,] <- c(length(herv.ranges), mean(widths), sum(widths))
}

plot.herv.width <- function(herv.ranges) {
  widths <- width(herv.ranges)
  breaks <- seq(0, 10000, by = 100)
  hist <- ggplot(as.data.frame(widths), aes(widths)) + geom_histogram(breaks = breaks) + labs(x = "length (bp)") + scale_y_continuous(labels = fancy_scientific)
  return(hist)
}

hervS2.width.hist <- plot.herv.width(hervS2.ranges)

g <- ggplot(as.data.frame(hervS1.lengths), aes(hervS1.lengths)) + geom_histogram(breaks = seq(0, 9000, by = 100)) + xlab("length in bp")# + ggtitle("HERV S1 element length distribution")
g <- g + theme(text = element_text(size=14))

pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_width_hist.pdf'), width = 6.3, height = 3)
hervS2.width.hist
dev.off()
