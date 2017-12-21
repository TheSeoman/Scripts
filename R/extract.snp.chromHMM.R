source('Scripts/R/paths.R')

load(PATHS$HERVS2.2KB.CHROMHMM.DATA)
load(PATHS$HERV.SNP.RANGES.DATA)
load(PATHS$CHROMHMM.SAMPLE.DATA)

require(GenomicRanges)

get.ranges.from.annotation <- function (annotation) {
  annotation.ranges <- lapply(annotation, function(sample){
    ann.list <- unlist(sample)
    ann.table <- data.frame(matrix(ann.list[!is.na(ann.list)], ncol=4, byrow = TRUE), stringsAsFactors = FALSE) 
    ann.table <- ann.table[!duplicated(ann.table),]
    ann.table[,2] <- as.numeric(ann.table[,2]) + 1
    ann.table[,3] <- as.numeric(ann.table[,3]) 
    ann.ranges <- GRanges(seqnames = ann.table[,1], ranges = IRanges(start = as.numeric(ann.table[,2]), end = as.numeric(ann.table[,3])), state = ann.table[,4])
    
    return(ann.ranges)
  })
  return(annotation.ranges)
}

get.snp.annotation <- function (snp.ranges, annotation.ranges, sample.ids) {
  snp.annotation.list <- lapply(annotation.ranges, function (sample.annotation.ranges) {
    overlap <- findOverlaps(snp.ranges, sample.annotation.ranges)
    return (sample.annotation.ranges[subjectHits(overlap)]$state)
  })
  snp.annotation <- data.frame(snp.annotation.list)
  rownames(snp.annotation) <- names(snp.ranges)
  colnames(snp.annotation) <- ids
  return(snp.annotation)
} 

# all snps in hervs
hervS2.2kb.annotation.ranges <- get.ranges.from.annotation(hervS2.2kb.annotation)
hervS2.2kb.snp.annotation <- get.snp.annotation(hervS2.2kb.snp.ranges, hervS2.2kb.annotation.ranges, ids)

save(hervS2.2kb.snp.annotation, file = PATHS$HERVS2.2KB.SNP.CHROMHMM.DATA)

hervS3.snp.annotation <- hervS2.2kb.snp.annotation[names(hervS3.snp.ranges),]

# snps in meQTLs related to hervs
load(PATHS$HERV.MEQTL.OVERLAP.DATA)
hervS1.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS1.meqtl.overlap$snp$snp),]
hervS1.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS1.meqtl.overlap$meth$cpg),]
hervS1.1kb.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS1.1kb.meqtl.overlap$snp$snp),]
hervS1.1kb.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS1.1kb.meqtl.overlap$meth$cpg),]
hervS1.2kb.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS1.2kb.meqtl.overlap$snp$snp),]
hervS1.2kb.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS1.2kb.meqtl.overlap$meth$cpg),]

hervS2.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS2.meqtl.overlap$snp$snp),]
hervS2.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS2.meqtl.overlap$meth$cpg),]
hervS2.1kb.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS2.1kb.meqtl.overlap$snp$snp),]
hervS2.1kb.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS2.1kb.meqtl.overlap$meth$cpg),]
hervS2.2kb.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS2.2kb.meqtl.overlap$snp$snp),]
hervS2.2kb.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS2.2kb.meqtl.overlap$meth$cpg),]

hervS3.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS3.meqtl.overlap$snp$snp),]
hervS3.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS3.meqtl.overlap$meth$cpg),]
hervS3.1kb.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS3.1kb.meqtl.overlap$snp$snp),]
hervS3.1kb.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS3.1kb.meqtl.overlap$meth$cpg),]
hervS3.2kb.meqtl.snp.annotation <- hervS2.2kb.snp.annotation[unique(hervS3.2kb.meqtl.overlap$snp$snp),]
hervS3.2kb.meqtl.meth.annotation <- hervS2.2kb.snp.annotation[unique(hervS3.2kb.meqtl.overlap$meth$cpg),]



save(hervS1.meqtl.snp.annotation, hervS1.meqtl.meth.annotation, hervS1.1kb.meqtl.snp.annotation, hervS1.1kb.meqtl.meth.annotation , hervS1.2kb.meqtl.snp.annotation, hervS1.2kb.meqtl.meth.annotation, 
     hervS2.meqtl.snp.annotation, hervS2.meqtl.meth.annotation, hervS2.1kb.meqtl.snp.annotation, hervS2.1kb.meqtl.meth.annotation , hervS2.2kb.meqtl.snp.annotation, hervS2.2kb.meqtl.meth.annotation,
     hervS3.meqtl.snp.annotation, hervS3.meqtl.meth.annotation, hervS3.1kb.meqtl.snp.annotation, hervS3.1kb.meqtl.meth.annotation , hervS3.2kb.meqtl.snp.annotation, hervS3.2kb.meqtl.meth.annotation,
     file = PATHS$HERV.MEQTL.CHROMHMM.ANNOTATION.DATA)
