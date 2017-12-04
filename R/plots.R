source('Scripts/R/paths.R')

require(reshape2)
require(ggplot2)
expression.overlaps <- read.table(file = paste0(PATHS$DATA.DIR, 'overlaps/hervS1.expression.tsv'), sep = '\t', header = TRUE, row.names = 1)
plot.data <- melt(expression.overlaps)
plot.data <- cbind(plot.data, type = factor(rownames(expression.overlaps), levels =c("direct", "1kb", "2kb")))

png(paste0(PATHS$PLOT.DIR, 'hervS1.expr.overlap.hist.png'), width = 500, height = 500)
p <- ggplot(plot.data, aes(variable, value))
p <- p + geom_bar(stat = "identity", aes(fill = type), position = "dodge") + xlab("") + ylab("Count") + ggtitle("Expression probe - HERV Set 1 overlap")
p + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16), legend.text=element_text(size=13), legend.title=element_text(size=14))
dev.off()

methylation.overlaps <- read.table(file = paste0(PATHS$DATA.DIR, 'overlaps/hervS1.methylation.tsv'), sep = '\t', header = TRUE, row.names = 1)
plot.data <- melt(methylation.overlaps)
plot.data <- cbind(plot.data, type = factor(rownames(methylation.overlaps), levels = c("direct", "1kb", "2kb")))

png(paste0(PATHS$PLOT.DIR, 'hervS1.meth.overlap.hist.png'), width = 500, height = 500)
p <- ggplot(plot.data, aes(variable, value))
p <- p + geom_bar(stat = "identity", aes(fill = type), position = "dodge") + xlab("") + ylab("Count") + ggtitle("Methylation probe - HERV Set 1 overlaps")
p + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16), legend.text=element_text(size=13), legend.title=element_text(size=14))
dev.off()

load(PATHS$EXPR.OVERLAP.DATA)
expr <- data.frame(matrix(nrow = length(expr.S1.overlap$essay.ranges)))
rownames(expr) <- rownames(expr.S1.overlap$essay.data)
expr$mean <- apply(expr.S1.overlap$essay.data, 1, mean)
expr$sd <- apply(expr.S1.overlap$essay.data, 1, sd)
expr$cvar <- expr$sd/expr$mean * 100

png(paste0(PATHS$PLOT.DIR, 'hervS1.expr.cvar.hist.png'), width = 500, height = 300)
gp <- ggplot(data=expr, aes(expr$cvar)) + geom_histogram(breaks = seq(floor(min(expr$cvar)), ceiling(max(expr$cvar)), by = 0.25)) 
gp <- gp + xlab("Coefficent of variance") + ylab("Frequency") + ggtitle("Expression probes overlapping with HERV set 1")
gp + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16))
dev.off()

png(paste0(PATHS$PLOT.DIR, 'hervS1.expr.cvar.point.png'), width = 500, height = 260)
g <- ggplot(data=expr, aes(x = mean, cvar)) 
g <- g + geom_point() + xlab("Mean expression") + ylab("Coefficient of variance")
g + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16))
dev.off()

load(PATHS$METH.OVERLAP.DATA)
meth <- data.frame(matrix(nrow = length(meth.S1.overlap$essay.ranges)))
meth$mean <- apply(meth.S1.overlap$essay.data, 1, mean, na.rm = T)
meth$sd <- apply(meth.S1.overlap$essay.data, 1, sd, na.rm = T)
meth$cvar <- meth$sd/meth$mean * 100

png(paste0(PATHS$PLOT.DIR, 'hervS1.meth.cvar.hist.png'), width = 500, height = 300)
gp <- ggplot(data=meth, aes(meth$cvar)) + geom_histogram(breaks = seq(floor(min(meth$cvar)), ceiling(max(meth$cvar)), by = 4)) 
gp <- gp + xlab("Coefficent of variance") + ylab("Frequency") + ggtitle("Methylation probes overlapping with HERV set 1")
gp + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16))
dev.off()

png(paste0(PATHS$PLOT.DIR, 'hervS1.meth.cvar.point.png'), width = 500, height = 260)
g <- ggplot(data=meth, aes(x = mean, cvar)) 
g <- g + geom_point() + xlab("Mean methylation") + ylab("Coefficient of variance")
g + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16))
dev.off()

load(PATHS$HERV.SNP.RANGES.DATA)
load(PATHS$HERV.DATA)

overlap <- findOverlaps(hervS1.ranges, hervS1.snp.ranges, type = 'any')
plot.data <- data.frame(table(queryHits(overlap)))
colnames(plot.data) = c('index', 'count')
length <- end(hervS1.ranges[as.numeric(as.character(plot.data$index))]) - start(hervS1.ranges[as.numeric(as.character(plot.data$index))])
plot.data <- cbind(plot.data, length = length)
plot.data$freq <- plot.data$count/plot.data$length
plot.data[plot.data$count > 40, "count"] <- 41
plot.data[plot.data$freq > 0.05, "freq"] <- 0.051 

png(paste0(PATHS$PLOT.DIR, 'hervS1.snp.length.hist.png'), width = 500, height = 350)
p <- ggplot(data = plot.data, aes(freq)) 
p <- p + geom_histogram(breaks = seq(0, 0.051, by = 0.001)) + ylab("Frequency") + xlab("Number of SNPs per bp") + ggtitle("Number of SNPs per base in HERV element") + geom_vline(xintercept = 0.05, colour = "grey")
p + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16))
dev.off()

png(paste0(PATHS$PLOT.DIR, 'hervS1.snp.hist.png'), width = 500, height = 350)
p <- ggplot(data = plot.data, aes(count)) 
p <- p + geom_histogram(breaks = c(0:41)) + ylab("Frequency") + xlab("Number of SNPs") + ggtitle("Number of SNPs per HERV element") + geom_vline(xintercept = 40, colour = "grey")
p + theme(axis.text=element_text(size=13), axis.title=element_text(size=14), plot.title=element_text(size=16))
dev.off()