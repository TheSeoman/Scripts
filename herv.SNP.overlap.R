DATA.DIR = '/media/data/Masterarbeit/data/'

F.SNPS = paste0(DATA.DIR, 'SNPs/full_sorted.bgz')
F.SAMPLES = paste0(DATA.DIR,"F4/individuals.txt")

HERV.DATA <- paste0(DATA.DIR, 'herv/ranges.RData')


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

scanSNPs <- function(ranges) {
  require(Rsamtools);
  require(data.table);
  
  # rather naive implementation for now	
  snpTabix <- TabixFile(file=F.SNPS);
  result <- scanTabix(snpTabix, param=ranges);
  message(paste0("scanTabix on ", F.SNPS, " for ", length(ranges), " ranges."))
  counts <- countTabix(snpTabix, param=ranges);
  message(paste0("countTabix on ", F.SNPS, " for ", length(ranges), " ranges."))
  result <- result[names(counts[counts>0])];
  result <- unlist(result);
  
  if(!is.null(result)) {
    # get data frame
    if(length(result) > 4000){
      cuts <- seq(0,length(result), by=floor(length(result)/50));
      data <- c();
      for(i in 2:length(cuts)){
        temp <- fread(paste(result[(cuts[i-1]+1):cuts[i]],collapse="\n"), 
                      sep="\t", data.table=F);
        data <- rbind(data,temp);
        message(paste0("read ", i-1, "/", length(cuts)-1, " chunks."));
      }
      temp <- read.table(text=result[(cuts[length(cuts)]+1):length(result)], 
                         sep="\t", colClasses="character");
      data <- rbind(data,temp);
      rm(temp);
    }
    else {
      data <- read.table(text=result, sep="\t", colClasses="character");
    }
    if(nrow(data) != length(result)){
      stop("sanity check failed.");
    }
    data <- unique(data);
    message(paste("Processed", nrow(data), "SNPs." ));
    
    rm(result);
    rm(counts);
    
    # process the genotype information to get integers/factors
    for(i in 6:ncol(data)){
      data[,i] <- factor(round(as.numeric(data[,i])));
    }
    message("Rounded genotype information.")
    
    #create colnames using individual codes
    ids <- read.table(F.SAMPLES, stringsAsFactors=F, colClasses="character");
    colnames(data)<- c("chr", "name", "pos", "orig", "alt", ids[,1])
    rownames(data) <- data$name;
    message("Created colnames.")
    
    return(list(snpInfo=data[,c(1,3,4,5)], snps=data[,6:ncol(data)]));
  } else {
    cat("No SNPs found in specified regions.\n")
    return(list());
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

load(HERV.DATA)

hervS1.filtered.ranges <- hervS1.ranges[grepl('^chr\\d+$', seqnames(hervS1.ranges))]
hervS1.1kb.filtered.ranges <- hervS1.1kb.ranges[grepl('^chr\\d+$', seqnames(hervS1.1kb.ranges))]
hervS2.filtered.ranges <- hervS2.ranges[grepl('^chr\\d+$', seqnames(hervS2.ranges))]

hervS3.filtered.ranges <- hervS3.ranges[grepl('^chr\\d+$', seqnames(hervS3.ranges))]

hervS1.SNPs <- scanSNPs(hervS1.filtered.ranges)
hervS1.SNPs$ranges <- create.granges.from.snpinfo(hervS1.SNPs$snpInfo)
save(hervS1.SNPs, file = paste0(DATA.DIR, 'SNPs/hervS1.SNP.RData'))
plot.hist.overlap.snps(hervS1.ranges, hervS1.SNPs$ranges)

hervS1.1kb.SNPs <- scanSNPs(hervS1.1kb.filtered.ranges)
hervS1.1kb.SNPs$ranges <- create.granges.from.snpinfo(hervS1.1kb.SNPs$snpInfo)

hervS2.SNPs <- scanSNPs(hervS2.filtered.ranges)

hervS3.SNPs <- scanSNPs(hervS3.filtered.ranges)

hervS3.SNP.ranges <- create.granges.from.snpinfo(hervS3.SNPs$snpInfo)



save(hervS3.SNPs, file = paste0(DATA.DIR, 'SNPs/hervS3.SNP.RData'))


hervS1.SNPs.present <- S1.SNPs$snps[apply(S1.SNPs$snps, 1, function(row) any(row != 0)),]


