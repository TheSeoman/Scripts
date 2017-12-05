source('Scripts/R/paths.R')

load(PATHS$HERVS2.CHROMHMM.DATA)
load(PATHS$HERV.SNP.RANGES.DATA)

require(GenomicRanges)

get.ranges.from.annotation <- function (annotation) {
  annotation.ranges <- lapply(annotation, function(sample){
    ann.list <- unlist(sample)
    ann.table = data.frame(matrix(ann.list[!is.na(ann.list)], ncol=4, byrow=TRUE), stringsAsFactors=FALSE) 
    ann.table = ann.table[!duplicated(ann.table),]
    ann.table[,2] = as.numeric(ann.table[,2])
    ann.table[,3] = as.numeric(ann.table[,3]) - 1 
    ann.ranges = GRanges(seqnames = ann.table[,1], ranges = IRanges(start = as.numeric(ann.table[,2]), end = as.numeric(ann.table[,3])), state = ann.table[,4])
    
    return(ann.ranges)
  })  
}

hervS2.annotation.ranges <- get.ranges.from.annotation(hervS2.annotation)

overlap <- findOverlaps(hervS2.snp.ranges, hervS2.annotation.ranges[[1]])

hervS2.snp.ranges$state = hervS2.annotation.ranges[[1]][subjectHits(overlap)]$state
