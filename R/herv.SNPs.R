# file paths not yet extracted

DATA.DIR = '~/data/'
KORA.DIR = '/storage/groups/groups_epigenereg/analyses/PV_K14115g_Heinig/'

F.SNPS = paste0(KORA.DIR, 'results/20160204/genoF4/dosage_combined/MAF001/full_sorted.bgz')
F.SAMPLES = paste0(KORA.DIR, 'results/20160204/genoF4/individuals.txt')

HERV.DATA <- paste0(DATA.DIR, 'herv/ranges.RData')

F.S3.SNPs <- paste0(DATA.DIR, 'SNPs/hervS3.SNP.RData')
F.S3.2kb.SNPs <- paste0(DATA.DIR, 'SNPs/hervS3.1kb.SNP.RData')

require(GenomicRanges)

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

load(HERV.DATA)
print('Filtering ranges.')
hervS3.filtered.ranges <- hervS3.ranges[grepl('^chr\\d+$', seqnames(hervS3.ranges))]
hervS3.2kb.filtered.ranges <- hervS3.2kb.ranges[grepl('^chr\\d+$', seqnames(hervS3.2kb.ranges))]

print('Start scanning for SNPs in hervS3.filtered.ranges')
hervS3.SNPs <- scanSNPs(hervS3.filtered.ranges)
message(paste0('Saving hervS3.SNPs to ', F.S3.SNPs))
save(hervS3.SNPs, file = F.S3.SNPs)

#print('Start scanning for SNPs in hervS3.2kb.filtered.ranges')
#hervS3.2kb.SNPs <- scanSNPs(hervS3.filtered.ranges)
#message(paste0('Saving hervS3.SNPs to ', F.S3.2kb.SNPs))
#save(hervS3.2kb.SNPs, file = F.S3.2kb.SNPs)
