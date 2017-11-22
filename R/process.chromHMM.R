source('Scripts/R/paths.R')

require(GenomicRanges)

load(PATHS$HERVS2.CHROMMHMM.DATA)
hervS2.elements = ls(hervS2.annotation)

load(PATHS$HERV.DATA)
hervS1.elements = sapply(hervS1.ranges, function(range) {
  irange = ranges(range)
  element = paste0(seqnames(range), ':', start(irange), '-', end(irange))
  return(element)
})

hervS3.elements = sapply(hervS3.ranges, function(range) {
  irange = ranges(range)
  element = paste0(seqnames(range), ':', start(irange), '-', end(irange))
  return(element)
})