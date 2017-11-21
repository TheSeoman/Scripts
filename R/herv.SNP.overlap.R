source('Scripts/R/paths.R')

split.ranges <- function (ranges, split.size) {
  out <- list()
  for (i in 0:floor(length(ranges)/split.size)) {
    out[i+1] <- ranges[(i*split.size+1):min((i+1)*split.size, length(ranges))]
  }
  return(out)
}

split.scan.combine <- function (ranges, split.size) {
  out <- list()
  split <- split.ranges(ranges, split.size)
  snps <- list()
  unique.snp.ids <- c()
  for (i in 1:length(split)) {
    snps[i] <- scanSNPs(split[i])
    unique.snp.ids <- unique(c(unique.snp.ids, rownames(snps[i]$snpInfo)))
  }
  
}

#' Scans genotype files for SNPs within the provided genomic ranges
#'
#' @param ranges GRanges object containing the ranges which to scan for SNPs
#' @param dir Base directory in which genotype information is stored
#' @param genotype.file Path to and including file which contains the genotypes,
#' relative to the base directory
#' @param id.file Path to and including file which contains the individual ids,
#' relative to the base directory
#' 
#' @author Johann Hawe
#'
scan.snps <- function(ranges) {
  library(data.table)
  
  # create system command using tabix...
  ranges.str <- paste(ranges, collapse=" ")
  
  cmd <- paste0("tabix ", F.SNP, " ", ranges.str)
  data <- try(fread(cmd, sep="\t", data.table=F), silent=T)
  if(inherits(data, "try-error")){
    cat("No SNPs found in specified regions.\n")
    return(list())
  } else {
    data <- data[!duplicated(data[,2]),,drop=F]
    message(paste("Processed", nrow(data), "SNPs." ))
    
    # process the genotype information to get integers/factors
    for(i in 6:ncol(data)){
      data[,i] <- factor(round(as.numeric(data[,i])))
    }
    
    ## create colnames using individual codes
    ids <- read.table(F.ID, stringsAsFactors=F, colClasses="character")
    colnames(data)<- c("chr", "name", "pos", "orig", "alt", ids[,1])
    rownames(data) <- data$name
    
    return(list(snpInfo=data[,c(1,3,4,5)], snps=data[,6:ncol(data)]))
  }
}

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

plot.hist.overlap.snps <- function (herv.ranges, snp.ranges) {
  overlap.hits <- findOverlaps(herv.ranges, snp.ranges, type = 'any')
  hist(table(as.factor(queryHits(overlap.hits))), breaks = c(1:300), xlim = c(0, 30))
}

load(PATHS$HERV.DATA)

hervS1.filtered.ranges <- hervS1.ranges[grepl('^chr\\d+$', seqnames(hervS1.ranges))]
hervS1.1kb.filtered.ranges <- hervS1.1kb.ranges[grepl('^chr\\d+$', seqnames(hervS1.1kb.ranges))]
hervS2.filtered.ranges <- hervS2.ranges[grepl('^chr\\d+$', seqnames(hervS2.ranges))]

hervS3.filtered.ranges <- hervS3.ranges[grepl('^chr\\d+$', seqnames(hervS3.ranges))]

hervS1.SNPs <- scanSNPs(hervS1.filtered.ranges)
hervS1.SNPs$ranges <- create.granges.from.snpinfo(hervS1.SNPs$snpInfo)
save(hervS1.SNPs, file = PATHS$HERVS1.SNP.DATA)
plot.hist.overlap.snps(hervS1.ranges, hervS1.SNPs$ranges)

hervS1.1kb.SNPs <- scanSNPs(hervS1.1kb.filtered.ranges)
hervS1.1kb.SNPs$ranges <- create.granges.from.snpinfo(hervS1.1kb.SNPs$snpInfo)

hervS2.SNPs <- scanSNPs(hervS2.filtered.ranges)

hervS3.SNPs <- scan.snps(hervS3.filtered.ranges)

hervS3.SNP.ranges <- create.granges.from.snpinfo(hervS3.SNPs$snpInfo)

save(hervS3.SNPs, file = PATHS$HERVS3.SNP.DATA)


hervS1.SNPs.present <- S1.SNPs$snps[apply(S1.SNPs$snps, 1, function(row) any(row != 0)),]


