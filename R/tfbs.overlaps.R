source('Scripts/R/paths.R')

require('GenomicRanges')

load(PATHS$EXPR.RANGES.DATA)
load(PATHS$METH.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)
load(PATHS$TFBS.RANGES.DATA)

enlarge.ranges <- function(ranges, flanking) {
  enlarged.ranges <- GRanges(
    seqnames = seqnames(ranges),
    ranges = IRanges(start = start(ranges) - flanking, end = end(ranges) + flanking),
    strand = strand(ranges),
    name = ranges$name,
    score = ranges$score
  )
  names(enlarged.ranges) <- names(ranges)
  return (enlarged.ranges)
}

get.tfbs.overlaps <- function (tfbs.ranges, essay.ranges, flanking, data.type) {
  out <- list()
  if( flanking > 0) {
    essay.ranges <- enlarge.ranges(essay.ranges, flanking)
  }
  overlap.hits <- findOverlaps(tfbs.ranges, essay.ranges, type = 'any')
  out$pairs <- cbind.data.frame(names(tfbs.ranges[queryHits(overlap.hits)]), names(essay.ranges[subjectHits(overlap.hits)]), stringsAsFactors = FALSE)
  colnames(out$pairs) <- c('tfbs.id', paste0(data.type, '.id'))
  out$tfbs.ranges <- tfbs.ranges[unique(names(tfbs.ranges[queryHits(overlap.hits)]))]
  out[[paste0(data.type, '.ranges')]] <- essay.ranges[unique(names(essay.ranges[subjectHits(overlap.hits)]))]
  return(out)
}

meth.tfbs.overlap <- get.tfbs.overlaps(blood.tfbs.ranges, meth.ranges, 100, 'meth')
save(meth.tfbs.overlap, file = PATHS$METH.TFBS.OVERLAP.DATA)

expr.tfbs.overlap <- get.tfbs.overlaps(blood.tfbs.ranges, expr.ranges, 0, 'expr')
save(expr.tfbs.overlap, file = PATHS$EXPR.TFBS.OVERLAP.DATA)

snp.tfbs.overlap <- get.tfbs.overlaps(blood.tfbs.ranges, snp.ranges, 0, 'snp')
save(snp.tfbs.overlap, file = PATHS$SNP.TFBS.OVERLAP.DATA)
  