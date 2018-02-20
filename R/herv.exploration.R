source('Scripts/R/paths.R')

load(PATHS$HERV.RANGES.DATA)

get.herv.lengths <- function(herv.ranges) {
  return(end(herv.ranges)-start(herv.ranges))
}

hervS1.lengths <- get.herv.lengths(hervS1.red.ranges)

hervS2.lengths <- get.herv.lengths(hervS2.ranges)

hervS1.lengths.df <- as.data.frame(hervS1.lengths)

g <- ggplot(as.data.frame(hervS1.lengths), aes(hervS1.lengths)) + geom_histogram(breaks = seq(0, 9000, by = 100)) + xlab("length in bp")# + ggtitle("HERV S1 element length distribution")
g <- g + theme(text = element_text(size=14))

pdf(file = paste0(PATHS$PLOT.DIR, 'hervS1_lengths_hist.pdf'), width = 6, height = 4)
g
dev.off()
