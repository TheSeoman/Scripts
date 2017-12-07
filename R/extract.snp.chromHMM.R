source('Scripts/R/paths.R')

load(PATHS$HERVS2.2KB.CHROMHMM.DATA)
load(PATHS$HERV.SNP.RANGES.DATA)
load(PATHS$CHROMHMM.SAMPLE.DATA)

require(GenomicRanges)

get.ranges.from.annotation <- function (annotation) {
  annotation.ranges <- lapply(annotation, function(sample){
    ann.list <- unlist(sample)
    ann.table <- data.frame(matrix(ann.list[!is.na(ann.list)], ncol<-4, byrow<-TRUE), stringsAsFactors<-FALSE) 
    ann.table <- ann.table[!duplicated(ann.table),]
    ann.table[,2] <- as.numeric(ann.table[,2]) + 1
    ann.table[,3] <- as.numeric(ann.table[,3]) 
    ann.ranges <- GRanges(seqnames <- ann.table[,1], ranges <- IRanges(start <- as.numeric(ann.table[,2]), end <- as.numeric(ann.table[,3])), state <- ann.table[,4])
    
    return(ann.ranges)
  })
  return(annotation.ranges)
}

get.snp.annotation <- function (snp.ranges, annotation.ranges, sample.ids) {
  snp.annotation.list <- lapply(annotation.ranges, function (sample.annotation.ranges) {
    overlap <- findOverlaps(snp.ranges, sample.annotation.ranges)
    return (sample.annotation.ranges[subjectHits(overlap)]$state)
  })
  snp.annotation <- cbind.data.frame(snp.annotation.list)
  rownames(snp.annotation) <- names(snp.ranges)
  colnames(snp.annotation) <- ids
  return(snp.annotation)
} 

hervS2.2kb.annotation.ranges <- get.ranges.from.annotation(hervS2.2kb.annotation)
hervS2.2kb.snp.annotation <- get.snp.annotation(hervS2.2kb.snp.ranges, hervS2.2kb.annotation.ranges, ids)

hervS3.snp.annotation <- hervS2.2kb.snp.annotation[names(hervS3.snp.ranges),]

