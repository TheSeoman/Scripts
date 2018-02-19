source('Scripts/R/paths.R')

load(PATHS$MAF001.RES.ME.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

load(PATHS$EQTM.ME.DATA)
load(PATHS$METH.RANGES.DATA)

load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$HERV.EQTM.OVERLAP.DATA)

library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)

get.chr.bin.ranges <- function(bin.size) {
  chrlen <- seqlengths(Hsapiens)
  
  chrs = paste0("chr", 1:22)
  chrlen <- chrlen[chrs]
  
  bin.start <- sapply(chrlen, function(c) {
    return(seq(1, c, by = bin.size))
  })
  
  seqnames <- unlist(lapply(chrs, function(chr) {
    return(rep(chr, length(bin.start[[chr]])))
  }))
  
  bin.ranges <- GRanges(seqnames = seqnames, ranges = IRanges(start = unlist(bin.start), width = bin.size))
  names(bin.ranges) <- 1:length(bin.ranges)
  
  return(bin.ranges)
}

get.pair.bins <- function(bin.ranges, pairs, e1.ranges, e2.ranges) {
  e1.bin.overlaps <- findOverlaps(e1.ranges[pairs[, 1]], bin.ranges)
  e2.bin.overlaps <- findOverlaps(e2.ranges[pairs[, 2]], bin.ranges)
  
  mappable.pairs <- intersect(queryHits(e1.bin.overlaps), queryHits(e2.bin.overlaps))
  
  e1.bin <- subjectHits(e1.bin.overlaps)[queryHits(e1.bin.overlaps) %in% mappable.pairs]
  e2.bin <- subjectHits(e2.bin.overlaps)[queryHits(e2.bin.overlaps) %in% mappable.pairs]
  
  pair.bins <- cbind.data.frame(e1.bin, e2.bin)
  
  bin.table <- data.frame(table(pair.bins), stringsAsFactors = F)
  bin.table[, 1] <- as.character(bin.table[, 1])
  bin.table[, 2] <- as.character(bin.table[, 2])
  
  for( i in as.character(1:length(bin.ranges))) {
    if (!i %in% bin.table[, 1]) {
      bin.table[dim(bin.table)[1] + 1, ] <- c(i, '1', 0)
    }
    if (!i %in% bin.table[, 2]) {
      bin.table[dim(bin.table)[1] + 1, ] <- c('1', i, 0)
    }
  }
  bin.table$Freq <- log10(as.numeric(bin.table$Freq))
  return(bin.table)
}

get.single.bins <- function(bin.ranges, e.ranges) {
  e.bin.overlaps <- findOverlaps(e.ranges, bin.ranges)
  e.bin <- table(subjectHits(e.bin.overlaps))
  e.bin <- log10(e.bin)
  return(e.bin)
}

plot.pair.bins <- function(pair.bins, e1.name, e2.name, title) {
  g <- ggplot(pair.bins, aes(e1.bin, e2.bin)) + geom_tile(aes(fill = Freq)) + theme(axis.text=element_blank(), text = element_text(size=20), legend.text = element_text(size=10), legend.title = element_text(size=12)) + xlab(e1.name) + ylab(e2.name)
  g <- g + scale_fill_gradient(low = '#DDDDDD', high = 'black', na.value = 'white')
  g <- g + ggtitle(title) + labs(fill = expression(log[10](Freq)))
  return(g)
}

bin.ranges <- get.chr.bin.ranges(1e7)

eqtl.pairs <- rbind.data.frame(eqtl.me$cis$eqtls[, 2:1], eqtl.me$trans$eqtls[, 2:1], stringsAsFactors = F)
eqtl.bins <- get.pair.bins(bin.ranges, eqtl.pairs, expr.ranges, snp.ranges)
expr.bins <- get.single.bins(bin.ranges, expr.ranges)
svg(filename=paste0(PATHS$PLOT.DIR, 'all.eqtl.heatmap'), width = 5, height = 4.5)
png(file=paste0(PATHS$PLOT.DIR, 'all.eqtl.heatmap.png'), width = 700, height = 650)
plot.pair.bins(eqtl.bins, 'Expression', 'Genotype', 'Global eQTL pair positions')
dev.off()

eqtm.expr.ids <- c(as.character(eqtm.me$cis$eqtls$gene), as.character(eqtm.me$trans$eqtls$gene))
eqtm.meth.ids <- c(as.character(eqtm.me$cis$eqtls$snps), as.character(eqtm.me$trans$eqtls$snps))
eqtm.pairs <- cbind.data.frame(eqtm.expr.ids, eqtm.meth.ids, stringsAsFactors = F)
eqtm.bins <- get.pair.bins(bin.ranges, eqtm.pairs, expr.ranges, meth.ranges)
png(file=paste0(PATHS$PLOT.DIR, 'all.eqtm.heatmap.png'), width = 700, height = 650)
plot.pair.bins(eqtm.bins, 'Expression', 'Methylation', 'Global eQTM pair positions')
dev.off()

hervS1.eqtl.expr.pairs <- rbind.data.frame(hervS1.eqtl.overlap$cis.expr[, 2:1], hervS1.eqtl.overlap$trans.expr[, 2:1], stringsAsFactors = F)
hervS1.eqtl.expr.bins <- get.pair.bins(bin.ranges, hervS1.eqtl.expr.pairs, expr.ranges, snp.ranges)

hervS1.eqtl.either.pairs <- rbind.data.frame(hervS1.eqtl.overlap$cis.either[, 2:1], hervS1.eqtl.overlap$trans.either[, 2:1], stringsAsFactors = F)
hervS1.eqtl.either.bins <- get.pair.bins(bin.ranges, hervS1.eqtl.either.pairs, expr.ranges, snp.ranges)

png(file=paste0(PATHS$PLOT.DIR, 'hervS1.either.eqtl.heatmap.png'), width = 700, height = 650)
plot.pair.bins(bin.table, 'Expression', 'Genotype', 'hervS1 eQTL pair positions')
dev.off()

hervS1.eqtm.either.pairs <- rbind.data.frame(hervS1.eqtm.overlap$cis.either[, 2:1], hervS1.eqtm.overlap$trans.either[, 2:1], stringsAsFactors = F)
hervS1.eqtm.either.bins <- get.pair.bins(bin.ranges, hervS1.eqtm.either.pairs, expr.ranges, meth.ranges)

png(file=paste0(PATHS$PLOT.DIR, 'hervS1.either.eqtm.heatmap.png'), width = 700, height = 650)
plot.pair.bins(bin.table, 'Expression', 'Genotype', 'hervS1 eQTM pair positions')
dev.off()


g + scale_x_discrete(limits = c('1', '300'))
