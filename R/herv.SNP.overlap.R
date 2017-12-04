source('Scripts/R/paths.R')


create.granges.from.snpinfo <- function(snpinfo) {
  return (
    GRanges(seqnames = snpinfo$chr,
            ranges = IRanges(
              names = rownames(snpinfo),
              start = as.numeric(snpinfo$pos), 
              width = 1
            ), 
            strand = c('*'),
            orig = Rle(snpinfo$orig),
            alt = Rle(snpinfo$alt)
          )
  )
}

create.snpinfo.from.granges <- function(snp.ranges) {
  snp.info <- cbind(as.character(seqnames(snp.ranges)), pos = start(ranges(snp.ranges)))
  snp.info <- cbind(snp.info, orig = as.character(snp.ranges$orig))
  snp.info <- cbind(snp.info, alt = as.character(snp.ranges$alt))
  rownames(snp.info) = names(snp.ranges)
  return(data.frame(snp.info))
}

get.overlap.snp.ranges <- function (herv.ranges, snp.ranges) {
  overlap.hits <- findOverlaps(herv.ranges, snp.ranges, type = 'any')
  herv.snp.ranges <- snp.ranges[unique(subjectHits(overlap.hits))]
  return (herv.snp.ranges)
}

load(PATHS$HERV.DATA)

load(PATHS$HERVS2.2KB.SNP.INFO.DATA)
hervS2.2kb.snp.ranges <- create.granges.from.snpinfo(hervS2.2kb.snp.info)

hervS1.snp.ranges <- get.overlap.snp.ranges(hervS1.ranges, hervS2.2kb.snp.ranges)
hervS1.snp.info <- create.snpinfo.from.granges(hervS1.snp.ranges)
hervS1.1kb.snp.ranges <- get.overlap.snp.ranges(hervS1.1kb.ranges, hervS2.2kb.snp.ranges)
hervS1.1kb.snp.info <- create.snpinfo.from.granges(hervS1.1kb.snp.ranges)
hervS1.2kb.snp.ranges <- get.overlap.snp.ranges(hervS1.2kb.ranges, hervS2.2kb.snp.ranges)
hervS1.2kb.snp.info <- create.snpinfo.from.granges(hervS1.2kb.snp.ranges)

hervS2.snp.ranges <- get.overlap.snp.ranges(hervS2.ranges, hervS2.2kb.snp.ranges)
hervS2.snp.info <- create.snpinfo.from.granges(hervS2.snp.ranges)
hervS2.1kb.snp.ranges <- get.overlap.snp.ranges(hervS2.1kb.ranges, hervS2.2kb.snp.ranges)  
hervS2.1kb.snp.info <- create.snpinfo.from.granges(hervS2.1kb.snp.ranges)

hervS3.snp.ranges <- get.overlap.snp.ranges(hervS3.ranges, hervS2.2kb.snp.ranges)
hervS3.snp.info <- create.snpinfo.from.granges(hervS3.snp.ranges)
hervS3.1kb.snp.ranges <- get.overlap.snp.ranges(hervS3.1kb.ranges, hervS2.2kb.snp.ranges)
hervS3.1kb.snp.info <- create.snpinfo.from.granges(hervS3.1kb.snp.ranges)
hervS3.2kb.snp.ranges <- get.overlap.snp.ranges(hervS3.2kb.ranges, hervS2.2kb.snp.ranges)
hervS3.2kb.snp.info <- create.snpinfo.from.granges(hervS3.2kb.snp.ranges)

save(hervS1.snp.info, hervS1.1kb.snp.info, hervS1.2kb.snp.info, hervS2.snp.info, hervS2.1kb.snp.info, 
     hervS2.2kb.snp.info, hervS3.snp.info, hervS3.1kb.snp.info, hervS3.2kb.snp.info, file = PATHS$HERV.SNP.INFO.DATA)

save(hervS1.snp.ranges, hervS1.1kb.snp.ranges, hervS1.2kb.snp.ranges, hervS2.snp.ranges, hervS2.1kb.snp.ranges, 
     hervS2.2kb.snp.ranges, hervS3.snp.ranges, hervS3.1kb.snp.ranges, hervS3.2kb.snp.ranges, file = PATHS$HERV.SNP.RANGES.DATA)
