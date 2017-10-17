DATA.DIR = '/media/data/Masterarbeit/data/'

F.SNPS = paste0(DATA.DIR, 'SNPs/full_sorted.bgz')

HERV.DATA <- paste0(DATA.DIR, 'herv/ranges.RData')


scanSNPs <- function(ranges) {
  require(Rsamtools);
  require(data.table);
  
  # rather naive implementation for now	
  snpTabix <- TabixFile(file=F.SNPS);
  result <- scanTabix(snpTabix, param=ranges);
  counts <- countTabix(snpTabix, param=ranges);
  result <- result[names(counts[counts>0])];
  result <- unlist(result);
  
  ranges <- ranges[grepl('chr\\d+$|chrX|chrY', seqnames(ranges))]
  
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
    
    #create colnames using individual codes
    ids <- read.table(paste0(DATA.DIR,"F4/individuals.txt"), stringsAsFactors=F, colClasses="character");
    colnames(data)<- c("chr", "name", "pos", "orig", "alt", ids[,1])
    rownames(data) <- data$name;
    
    return(list(snpInfo=data[,c(1,3,4,5)], snps=data[,6:ncol(data)]));
  } else {
    cat("No SNPs found in specified regions.\n")
    return(list());
  }
}
load(HERV.DATA)

hervS1.filtered.ranges <- hervS1.ranges[grepl('^chr\\d+$', seqnames(hervS1.ranges))]

S1.SNPs <- scanSNPs(hervS1.filtered.ranges)



