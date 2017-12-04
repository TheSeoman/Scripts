source('Scripts/R/paths.R')

require(GenomicRanges)
require(illuminaHumanv3.db)
require(FDb.InfiniumMethylation.hg19)

get.expression.ranges <- function () {
  require(illuminaHumanv3.db)
  allLocs <- unlist(as.list(illuminaHumanv3GENOMICLOCATION))
  start <- as.numeric(unlist(lapply(allLocs , function(x)
    strsplit(as.character(x),":")[[1]][2])));
  validPos <- !is.na(start);
  start <- start[validPos];

  chrs <- unlist(lapply(allLocs, function(x)
    strsplit(as.character(x),":")[[1]][1]))[validPos];
  end <- as.numeric(unlist(lapply(allLocs , function(x)
    strsplit(as.character(x),":")[[1]][3])))[validPos];
  strand <- substr(unlist(lapply(allLocs , function(x)
    strsplit(as.character(x),":")[[1]][4])), 1, 1)[validPos];
  ids <- names(allLocs)[validPos];
  gr <- GRanges(chrs, ranges=IRanges(start,end), strand=strand, ids);
  return(gr)
}

calc.overlap.data <- function (herv.ranges, essay.ranges, essay.data) {
  overlap.hits <- findOverlaps(herv.ranges, essay.ranges, type = 'any')
  herv.overlap.ranges <- herv.ranges[unique(queryHits(overlap.hits))]
  essay.overlap.ranges <- essay.ranges[unique(subjectHits(overlap.hits))]
  if (!is.null(essay.ranges$ids)) {
    essay.overlap.ids <- unique(intersect(essay.overlap.ranges$ids, rownames(essay.data)))
  } else {
    essay.overlap.ids <- unique(intersect(names(essay.overlap.ranges), rownames(essay.data)))
  }
  essay.overlap.data <- essay.data[essay.overlap.ids,]
  out <- list()
  out$hits <- overlap.hits
  out$herv.ranges <- herv.overlap.ranges
  out$essay.ranges <- essay.overlap.ranges
  out$essay.data <- essay.overlap.data
  return(out)
}

combine.overlaps <- function (herv.ranges, essay1.ranges, essay2.ranges ,overlap1, overlap2) {
  herv.indices <- intersect(queryHits(overlap1$hits), queryHits(overlap2$hits))
  essay1.indices <- unique(subjectHits(overlap1$hits)[queryHits(overlap1$hits) %in% herv.indices])
  essay2.indices <- unique(subjectHits(overlap2$hits)[queryHits(overlap2$hits) %in% herv.indices])
  out <- list()
  out$herv.ranges <- herv.ranges[herv.indices]
  out$essay1.ranges <- essay1.ranges[essay1.indices]
  out$essay2.ranges <- essay2.ranges[essay2.indices]
  out$essay1.data <- overlap1$essay.data[unique(out$essay1.ranges$ids),]
  out$essay2.data <- overlap2$essay.data[unique(names(out$essay2.ranges)),]
  return (out)
}

print.overlap.info <- function(overlap) {
  message(paste0("Overlap info:\n# Overlaps: ", length(overlap$hits), 
                 "\n# hERVs: ", length(overlap$herv.ranges),
                 "\n# probes: ", dim(overlap$essay.data)[1]))
}

create.granges <- function (annotation, flanking = 0) {
  return (
    GRanges(
      seqnames = annotation$Chr,
      ranges = IRanges(annotation$start - flanking, end = annotation$end + flanking),
      strand = annotation$strand
    )
  )
}

plot.coefficient.of.variance <- function(data, data.type, main, file) {
  mean <- apply(data, 1, mean, na.rm = T)
  sd <- apply(data, 1, sd, na.rm = T)
  cvar <- sd/mean * 100
  hist(cvar, xlab = 'coefficient of variance', main = main, breaks = 50)
  plot(mean, cvar, xlab = paste0('mean ', data.type), ylab = 'coefficent of variance', pch = 3) 
}

hervS1.annot <- read.table(PATHS$HERVS1.ANNOT, sep = "\t", col.names = c('Chr', 'start', 'end', 'strand', 'family', 'class', 'superfamily'))
hervS2.annot <- read.table(PATHS$HERVS2.ANNOT, sep = "\t", col.names = c('Chr', 'start', 'end', 'strand', 'family', 'class', 'superfamily'))
hervS3.annot <- read.table(PATHS$HERVS3.ANNOT, sep = "\t", col.names = c('Chr', 'start', 'end', 'strand', 'family', 'class', 'superfamily'))

hervS1.ranges <- create.granges(hervS1.annot)
hervS2.ranges <- create.granges(hervS2.annot)
hervS3.ranges <- create.granges(hervS3.annot)

hervS1.1kb.ranges <- create.granges(hervS1.annot, 1000)
hervS2.1kb.ranges <- create.granges(hervS2.annot, 1000)
hervS3.1kb.ranges <- create.granges(hervS3.annot, 1000)

hervS1.2kb.ranges <- create.granges(hervS1.annot, 2000)
hervS2.2kb.ranges <- create.granges(hervS2.annot, 2000)
hervS3.2kb.ranges <- create.granges(hervS3.annot, 2000)

save(hervS1.ranges, hervS2.ranges, hervS3.ranges, hervS1.1kb.ranges, hervS2.1kb.ranges, hervS3.1kb.ranges, hervS1.2kb.ranges, hervS2.2kb.ranges, hervS3.2kb.ranges, file = PATHS$HERV.DATA)

expr.ranges <- get.expression.ranges()
load(PATHS$EXPR.DATA)
expr.data <- f4.norm

#meth.ranges <- features(FDb.InfiniumMethylation.hg19)
meth.ranges <- getPlatform(platform='HM450', genome='hg19')
meth.ranges <- meth.ranges[grep('cg|ch', meth.ranges$probeType)]
load(PATHS$METH.DATA)
meth.data <- data.frame(transform(beta))
rm(beta)

expr.S1.overlap <- calc.overlap.data(hervS1.ranges, expr.ranges, expr.data)
expr.S2.overlap <- calc.overlap.data(hervS2.ranges, expr.ranges, expr.data)
expr.S3.overlap <- calc.overlap.data(hervS3.ranges, expr.ranges, expr.data)

print.overlap.info(expr.S1.overlap)
print.overlap.info(expr.S2.overlap)
print.overlap.info(expr.S3.overlap)

pdf('/media/data/Masterarbeit/Plots/herv.expr.var.pdf', width = 8, height = 8)
layout(matrix(c(1:4), 2, 2, byrow = F))
plot.coefficient.of.variance(expr.S1.overlap$essay.data, 'expression', 'Probes overlapping with hERV set 1', '')
plot.coefficient.of.variance(expr.S2.overlap$essay.data, 'expression', 'Probes overlapping with hERV set 2', '')
dev.off()

expr.S1.1kb.overlap <- calc.overlap.data(hervS1.1kb.ranges, expr.ranges, expr.data)
expr.S2.1kb.overlap <- calc.overlap.data(hervS2.1kb.ranges, expr.ranges, expr.data)
expr.S3.1kb.overlap <- calc.overlap.data(hervS3.1kb.ranges, expr.ranges, expr.data)

print.overlap.info(expr.S1.1kb.overlap)
print.overlap.info(expr.S2.1kb.overlap)
print.overlap.info(expr.S3.1kb.overlap)

annotation[annotation$ProbeID %in% expr.S1.overlap$essay.ranges$ids,]
annotation[annotation$ProbeID %in% expr.S1.1kb.overlap$essay.ranges$ids,]
annotation[annotation$ProbeID %in% expr.S1.2kb.overlap$essay.ranges$ids,]

pdf('/media/data/Masterarbeit/Plots/herv.expr.1kb.var.pdf', width = 8, height = 8)
layout(matrix(c(1:4), 2, 2, byrow = F))
plot.coefficient.of.variance(expr.S1.1kb.overlap$essay.data, 'expression', 'Probes overlapping with hERV set 1 + flanking 1kb', '')
plot.coefficient.of.variance(expr.S2.1kb.overlap$essay.data, 'expression', 'Probes overlapping with hERV set 2 + flanking 1kb', '')
dev.off()

expr.S1.2kb.overlap <- calc.overlap.data(hervS1.2kb.ranges, expr.ranges, expr.data)
expr.S2.2kb.overlap <- calc.overlap.data(hervS2.2kb.ranges, expr.ranges, expr.data)
expr.S3.2kb.overlap <- calc.overlap.data(hervS3.2kb.ranges, expr.ranges, expr.data)

print.overlap.info(expr.S1.2kb.overlap)
print.overlap.info(expr.S2.2kb.overlap)
print.overlap.info(expr.S3.2kb.overlap)

save(expr.S1.overlap, expr.S2.overlap, expr.S3.overlap, expr.S1.1kb.overlap, expr.S2.1kb.overlap, expr.S3.1kb.overlap, expr.S1.2kb.overlap, expr.S2.2kb.overlap, expr.S3.2kb.overlap, file = PATHS$EXPR.OVERLAP.DATA))


meth.S1.overlap <- calc.overlap.data(hervS1.ranges, meth.ranges, meth.data)
meth.S2.overlap <- calc.overlap.data(hervS2.ranges, meth.ranges, meth.data)
meth.S3.overlap <- calc.overlap.data(hervS3.ranges, meth.ranges, meth.data)

print.overlap.info(meth.S1.overlap)
print.overlap.info(meth.S2.overlap)
print.overlap.info(meth.S3.overlap)

pdf('/media/data/Masterarbeit/Plots/herv.meth.var.pdf', width = 8, height = 8)
layout(matrix(c(1:4), 2, 2, byrow = F))
plot.coefficient.of.variance(meth.S1.overlap$essay.data, 'methylation', 'Probes overlapping with hERV set 1', '')
plot.coefficient.of.variance(meth.S2.overlap$essay.data, 'methylation', 'Probes overlapping with hERV set 2', '')
dev.off()

meth.S1.1kb.overlap <- calc.overlap.data(hervS1.1kb.ranges, meth.ranges, meth.data)
meth.S2.1kb.overlap <- calc.overlap.data(hervS2.1kb.ranges, meth.ranges, meth.data)
meth.S3.1kb.overlap <- calc.overlap.data(hervS3.1kb.ranges, meth.ranges, meth.data)

print.overlap.info(meth.S1.1kb.overlap)
print.overlap.info(meth.S2.1kb.overlap)
print.overlap.info(meth.S3.1kb.overlap)

meth.S1.2kb.overlap <- calc.overlap.data(hervS1.2kb.ranges, meth.ranges, meth.data)
meth.S2.2kb.overlap <- calc.overlap.data(hervS2.2kb.ranges, meth.ranges, meth.data)
meth.S3.2kb.overlap <- calc.overlap.data(hervS3.2kb.ranges, meth.ranges, meth.data)

print.overlap.info(meth.S1.2kb.overlap)
print.overlap.info(meth.S2.2kb.overlap)
print.overlap.info(meth.S3.2kb.overlap)

save(meth.S1.overlap, meth.S2.overlap, meth.S3.overlap, meth.S1.1kb.overlap, meth.S2.1kb.overlap, meth.S3.1kb.overlap, meth.S1.2kb.overlap, meth.S2.2kb.overlap, meth.S3.2kb.overlap, file = PATHS$METH.OVERLAP.DATA)

both.S1.overlap <- combine.overlaps(hervS1.ranges, expr.ranges, meth.ranges, expr.S1.overlap, meth.S1.overlap)
both.S2.overlap <- combine.overlaps(hervS2.ranges, expr.ranges, meth.ranges, expr.S2.overlap, meth.S2.overlap)
both.S3.overlap <- combine.overlaps(hervS3.ranges, expr.ranges, meth.ranges, expr.S3.overlap, meth.S3.overlap)

pdf('/media/data/Masterarbeit/Plots/herv.meth.var.pdf', width = 8, height = 8)
layout(matrix(c(1:4), 2, 2, byrow = F))
plot.coefficient.of.variance(meth.S1.overlap$essay.data, 'methylation', 'Probes overlapping with hERV set 1', '')
plot.coefficient.of.variance(meth.S2.overlap$essay.data, 'methylation', 'Probes overlapping with hERV set 2', '')
dev.off()

both.S1.1kb.overlap <- combine.overlaps(hervS1.ranges, expr.ranges, meth.ranges, expr.S1.1kb.overlap, meth.S1.1kb.overlap)
both.S2.1kb.overlap <- combine.overlaps(hervS2.ranges, expr.ranges, meth.ranges, expr.S2.1kb.overlap, meth.S2.1kb.overlap)
both.S3.1kb.overlap <- combine.overlaps(hervS3.ranges, expr.ranges, meth.ranges, expr.S3.1kb.overlap, meth.S3.1kb.overlap)

both.S1.2kb.overlap <- combine.overlaps(hervS1.ranges, expr.ranges, meth.ranges, expr.S1.2kb.overlap, meth.S1.2kb.overlap)
both.S2.2kb.overlap <- combine.overlaps(hervS2.ranges, expr.ranges, meth.ranges, expr.S2.2kb.overlap, meth.S2.2kb.overlap)
both.S3.2kb.overlap <- combine.overlaps(hervS3.ranges, expr.ranges, meth.ranges, expr.S3.2kb.overlap, meth.S3.2kb.overlap)

plot.coefficient.of.variance(meth.S1.overlap$essay.data, 'methylation', 'Probes overlapping with hERV set 1', '')
