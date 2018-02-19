source('Scripts/R/paths.R')

load(PATHS$MAF001.RES.ME.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)

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
  bin.table$Freq <- log10(bin.table$Freq)
  
  return(bin.table)
}

get.single.bins <- function(bin.ranges, e.ranges) {
  e.bin.overlaps <- findOverlaps(e.ranges, bin.ranges)
  e.bin <- table(subjectHits(e.bin.overlaps))
  e.bin <- log10(e.bin)
  return(e.bin)
}

plot.pair.bins <- function(pair.bins, e1.name, e2.name) {
  g <- ggplot(eqtl.bins, aes(e1.bin, e2.bin)) + geom_tile(aes(fill = Freq)) + theme(axis.text=element_blank()) + xlab(e1.name) + ylab(e2.name)
  return(g)
}

bin.ranges <- get.chr.bin.ranges(1e7)


eqtl.pairs <- rbind.data.frame(eqtl.me$cis$eqtls[, 2:1], eqtl.me$trans$eqtls[, 2:1])
eqtl.bins <- get.pair.bins(bin.ranges, eqtl.pairs, expr.ranges, snp.ranges)

expr.bins <- get.single.bins(bin.ranges, expr.ranges)

plot.pair.bins(eqtl.bins, 'Expression', 'Genotype')
