source('Scripts/R/paths.R')

library(data.table)

scan.snps <- function(ranges) {

  # create system command using tabix...
  ranges.str <- paste(sapply(ranges, function(range) {
    irange = ranges(range)
    element = paste0(seqnames(range), ':', start(irange), '-', end(irange))
    return(element)
  }), collapse=" ")
  
  cmd <- paste0("tabix ", PATHS$F.SNP, " ", ranges.str)
  data <- try(fread(cmd, sep="\t", data.table=F), silent=F)
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
    ids <- read.table(PATHS$F.SNP.SAMPLES, stringsAsFactors=F, colClasses="character")
    colnames(data)<- c("chr", "name", "pos", "orig", "alt", ids[,1])
    rownames(data) <- data$name
    
    return(list(snpInfo=data[,c(1,3,4,5)], snps=data[,6:ncol(data)]))
  }
}

load(PATHS$HERV.S2.2KB.DATA)

args <- commandArgs(TRUE)
batch.number <- as.integer(args[1])
batch.count <- as.integer(args[2])
batch.size <- ceiling(length(hervS2.2kb.ranges) / batch.count)

ranges <- hervS2.2kb.ranges[((batch.number-1)*batch.size+1):(min(batch.number*batch.size, length(hervS2.2kb.ranges)))]

snps <- scan.snps(ranges)

if(length(snps) > 0) {
  saveRDS(snps, file = paste0(PATHS$DATA.DIR, 'SNPs/S2.2kb/batch', batch.number, '.rds'))
}


