source('Scripts/R/paths.R')

load(PATHS$MAF001.RES.ME.DATA)
load(PATHS$EXPR.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

load(PATHS$EQTM.ME.DATA)
load(PATHS$METH.RANGES.DATA)

load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$HERV.EQTM.OVERLAP.DATA)

load(PATHS$HERV.EXPR.OVERLAP.DATA)
load(PATHS$HERV.SNP.OVERLAP.DATA)

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
  
  bin.names <- as.character(1:length(bin.ranges))
  missing.e1.bins <- bin.names[!bin.names %in% bin.table$e1.bin]
  missing.e2.bins <- bin.names[!bin.names %in% bin.table$e2.bin]
  
  e1.add <- expand.grid(missing.e1.bins, bin.names)
  e1.add$Freq <- 0
  colnames(e1.add) <- c('e1.bin', 'e2.bin', 'Freq')
  
  e2.add <- expand.grid(bin.names[!bin.names %in% missing.e1.bins], missing.e2.bins)
  e2.add$Freq <- 0
  colnames(e2.add) <- c('e1.bin', 'e2.bin', 'Freq')
  
  full.bin.table <- rbind.data.frame(bin.table, e1.add)
  full.bin.table <- rbind.data.frame(full.bin.table, e2.add)
  full.bin.table[, 1] <- as.numeric(full.bin.table[, 1])
  full.bin.table[, 2] <- as.numeric(full.bin.table[, 2])
  full.bin.table <- full.bin.table[order(full.bin.table$e1.bin, full.bin.table$e2.bin),]
  full.bin.table$Freq <- as.numeric(full.bin.table$Freq)
  return(full.bin.table)
}

get.single.bins <- function(bin.ranges, e.ranges) {
  e.bin.overlaps <- findOverlaps(e.ranges, bin.ranges)
  e.bin <- table(subjectHits(e.bin.overlaps))
  bin.names <- as.character(1:length(bin.ranges))
  missing.bins <- bin.names[!bin.names %in% names(e.bin)]
  e.bin[missing.bins] <- 0
  e.bin <- e.bin[order(as.numeric(names(e.bin)))]
  # e.bin <- log10(e.bin)
  return(e.bin)
}

get.potential.pair.bins <- function(all.e1.bins, all.e2.bins, herv.e1.bins = NULL, herv.e2.bins = NULL) {
  potential.pair.bins <- expand.grid(1:300, 1:300)
  potential.pair.bins[, 3] <- 0
  colnames(potential.pair.bins) <- c('e1.bin', 'e2.bin', 'Freq')
  if (!is.null(herv.e1.bins) & !is.null(herv.e2.bins)) {
    for( i in 1:dim(potential.pair.bins)[1]) {
      e1.bin <- potential.pair.bins[i, 1]
      e2.bin <- potential.pair.bins[i, 2]
      potential.pair.bins[i, 3] <- ((all.e1.bins[e1.bin] - herv.e1.bins[e1.bin]) * all.e2.bins[e2.bin]) + (all.e1.bins[e1.bin] * (all.e2.bins[e2.bin] - herv.e2.bins[e2.bin]))
    }
  } else {
    for( i in 1:dim(potential.pair.bins)[1]) {
      e1.bin <- potential.pair.bins[i, 1]
      e2.bin <- potential.pair.bins[i, 2]
      potential.pair.bins[i, 3] <- all.e1.bins[e1.bin] * all.e2.bins[e2.bin]
    }  
  }
  return(potential.pair.bins)
}

normalize.pair.bins <- function(pair.bins, potential.pair.bins) {
  for(i in 1:dim(pair.bins)[1]) {
    if(pair.bins[i, 3] != 0){
      pair.bins[i, 3] <- pair.bins[i, 3] / potential.pair.bins[i, 3]
    }
  }
  pair.bins$Freq <- log10(pair.bins$Freq)
  return(pair.bins)
}

plot.pair.bins <- function(pair.bins, e1.name, e2.name, title) {
  g <- ggplot(pair.bins, aes(e1.bin, e2.bin)) + geom_tile(aes(fill = Freq)) + theme(axis.text=element_blank(), text = element_text(size=14), legend.text = element_text(size=10), legend.title = element_text(size=12)) + xlab(e1.name) + ylab(e2.name)
  g <- g + scale_fill_gradient(low = 'white', high = 'red', na.value = 'lightgrey')
  g <- g + ggtitle(title) + labs(fill = expression(log[10](Frac)))
  return(g)
}

bin.ranges <- get.chr.bin.ranges(1e7)

eqtl.pairs <- rbind.data.frame(eqtl.me$cis$eqtls[, 2:1], eqtl.me$trans$eqtls[, 2:1], stringsAsFactors = F)
eqtl.bins <- get.pair.bins(bin.ranges, eqtl.pairs, expr.ranges, snp.ranges)
eqtl.potential.bins <- get.potential.pair.bins(expr.bins, snp.bins)

norm.eqtl.bins <- eqtl.bins
norm.eqtl.bins$Freq <- eqtl.bins$Freq/eqtl.potential.bins$Freq

expr.bins <- get.single.bins(bin.ranges, expr.ranges)
snp.bins <- get.single.bins(bin.ranges, snp.ranges)

norm.eqtl.bins <- normalize.pair.bins(eqtl.bins, expr.bins, snp.bins)

pdf(file=paste0(PATHS$PLOT.DIR, 'all.norm.eqtl.heatmap.pdf'), width = 5, height = 4.5)
png(file=paste0(PATHS$PLOT.DIR, 'all.norm.eqtl.heatmap.png'), width = 500, height = 450)
plot.pair.bins(eqtl.bins, 'Expression', 'Genotype', 'Fractions of significant SNP-expression probe pairs')
dev.off()

eqtm.expr.ids <- c(as.character(eqtm.me$cis$eqtls$gene), as.character(eqtm.me$trans$eqtls$gene))
eqtm.meth.ids <- c(as.character(eqtm.me$cis$eqtls$snps), as.character(eqtm.me$trans$eqtls$snps))
eqtm.pairs <- cbind.data.frame(eqtm.expr.ids, eqtm.meth.ids, stringsAsFactors = F)
eqtm.bins <- get.pair.bins(bin.ranges, eqtm.pairs, expr.ranges, meth.ranges)
png(file=paste0(PATHS$PLOT.DIR, 'all.eqtm.heatmap.png'), width = 700, height = 650)
plot.pair.bins(eqtm.bins, 'Expression', 'Methylation', 'Global eQTM pair positions')
dev.off()

hervS2.eqtl.either.pairs <- rbind.data.frame(hervS2.eqtl.overlap$cis.either[, 2:1], hervS2.eqtl.overlap$trans.either[, 2:1], stringsAsFactors = F)
hervS2.eqtl.either.bins <- get.pair.bins(bin.ranges, hervS2.eqtl.either.pairs, expr.ranges, snp.ranges)

hervS2.expr.bins <- get.single.bins(bin.ranges, hervS2.expr.overlap$expr.ranges)
hervS2.snp.bins <- get.single.bins(bin.ranges, hervS2.snp.overlap$snp.ranges)

hervS2.eqtl.potential.bins <- get.potential.pair.bins(expr.bins, snp.bins, hervS2.expr.bins, hervS2.snp.bins)

norm.hervS2.eqtl.either.bins <- hervS2.eqtl.either.bins
norm.hervS2.eqtl.either.bins$Freq <- hervS2.eqtl.either.bins$Freq/hervS2.eqtl.potential.bins$Freq
norm.hervS2.eqtl.either.bins[is.na(norm.hervS2.eqtl.either.bins$Freq), 'Freq'] <- NA
norm.hervS2.eqtl.either.bins$Freq <- log10(norm.hervS2.eqtl.either.bins$Freq)

png(file=paste0(PATHS$PLOT.DIR, 'hervS1.either.eqtl.heatmap.png'), width = 700, height = 650)
plot.pair.bins(norm.hervS2.eqtl.either.bins, 'Expression', 'Genotype', 'B')
dev.off()


png(file=paste0(PATHS$PLOT.DIR, 'hervS1.either.eqtm.heatmap.png'), width = 700, height = 650)
plot.pair.bins(bin.table, 'Expression', 'Genotype', 'hervS1 eQTM pair positions')
dev.off()


g + scale_x_discrete(limits = c('1', '300'))
