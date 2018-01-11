source('Scripts/R/paths.R')

## annotate with the roadmap chromHMM states
load(PATHS$CHROMHMM.SAMPLE.DATA)

chromHMM.range.annotation <- function (ranges, id, dir = paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/"), suffix = "_15_coreMarks_mnemonics.bed.bgz") {
  library(Rsamtools)
  library(GenomicRanges)
  message(paste0('Processing ', id))
  file = paste0(dir, id, suffix)
  avail = as.logical(seqnames(ranges) %in% seqnamesTabix(file))
  avail.ann = scanTabix(file, param=ranges[avail])
  ann.split = sapply(avail.ann, function(x) strsplit(x, "\t") )  
  return (ann.split)
}

load(PATHS$SNP.RANGES.DATA)
args <- commandArgs(TRUE)
id <- ids[as.integer(args[1])]
message(paste0('Calculate chromHMM states for snp.ranges in ', id))
annotation = chromHMM.range.annotation(snp.ranges, id)

save(annotation, file = paste0(PATHS$CHROMHMM.OUT.DIR, 'snp.ranges/', id, '.chromHMM.RData'))


