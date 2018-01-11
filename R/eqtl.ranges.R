source('Scripts/R/paths.R')

load(PATHS$HERV.MAF001.ME)
get.expression.ranges <- function () {
  require(illuminaHumanv3.db)
  allLocs <- unlist(as.list(illuminaHumanv3GENOMICLOCATION))
  start <- as.numeric(unlist(lapply(allLocs , function(x)
    strsplit(as.character(x),":")[[1]][2])));
  validPos <- !is.na(start);
  start <- start[validPos];
  
  chrs <- unlist(lapply(allLocs, function(x)
    strsplit(as.character(x),":")[[1]][1]))[validPos];
  end <- as.numeric(unlist(lapply(allLocs , function(x)
    strsplit(as.character(x),":")[[1]][3])))[validPos];
  strand <- substr(unlist(lapply(allLocs , function(x)
    strsplit(as.character(x),":")[[1]][4])), 1, 1)[validPos];
  ids <- names(allLocs)[validPos];
  gr <- GRanges(chrs, ranges=IRanges(start,end), strand=strand);
  names(gr) <- ids
  return(gr)
}

get.snp.ranges <- function(){
  require(GenomicRanges)
  snp.pos <- read.table(PATHS$F.SNP.POS, sep = '\t', header = TRUE)
  snp.ranges <- GRanges(seqnames = snp.pos$chr, ranges <- IRanges(start = as.numeric(snp.pos$pos), width = 1))  
  names(snp.ranges) <- snp.pos$snp
  return(snp.ranges)
}

expr.ranges <- get.expression.ranges()
save(expr.ranges, file = PATHS$EXPR.RANGES.DATA)

snp.ranges <- get.snp.ranges()
save(snp.ranges, file = PATHS$SNP.RANGES.DATA)

all.snp.ids <- unique(c(as.character(me$cis$eqtls$snps), as.character(me$trans$eqtls$snps)))
eqtl.expr.ids <- unique(c(as.character(me$cis$eqtls$gene), as.character(me$trans$eqtls$gene)))
eqtl.expr.ranges <- expr.ranges[expr.ranges$ids %in% eqtl.expr.ids]
