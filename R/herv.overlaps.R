source('Scripts/R/paths.R')

require(GenomicRanges)
require(illuminaHumanv3.db)
require(FDb.InfiniumMethylation.hg19)
require(rtracklayer)

if(!file.exists(PATHS$EXPR.RANGES.DATA)) {
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
    gr <- GRanges(chrs, ranges=IRanges(start,end), strand=strand)
    return(gr)
  }
  expr.ranges <- get.expression.ranges()
  save(expr.ranges, file = PATHS$EXPR.RANGES.DATA)
} else {
  load(PATHS$EXPR.RANGES)
}

if (!file.exists(PATHS$METH.RANGES.DATA)) {
  meth.ranges <- getPlatform(platform = 'HM450', genome = 'hg19')
  meth.ranges <- meth.ranges[grep('cg|ch', meth.ranges$probeType)]
  save(meth.ranges, file = PATHS$METH.RANGES.DATA)
} else {
  load(PATHS$METH.RANGES.DATA)
}

load(PATHS$SNP.RANGES.DATA)
load(PATHS$TFBS.RANGES.DATA)

calc.overlap.data <- function (herv.ranges, essay.ranges, essay.data, data.type) {
  overlap.hits <- findOverlaps(herv.ranges, essay.ranges, type = 'any')
  herv.overlap.ranges <- herv.ranges[unique(queryHits(overlap.hits))]
  essay.overlap.ranges <- essay.ranges[unique(subjectHits(overlap.hits))]
  out <- list()
  out$pairs <- cbind.data.frame(names(herv.ranges[queryHits(overlap.hits)]), names(essay.ranges[subjectHits(overlap.hits)]), stringsAsFactors = FALSE)
  colnames(out$pairs) <- c('herv.id', paste0(data.type, '.id'))
  out$herv.ranges <- herv.overlap.ranges
  out[[paste0(data.type, '.ranges')]] <- essay.overlap.ranges
  if(!is.null(essay.data)) {
    essay.overlap.ids <- unique(intersect(names(essay.overlap.ranges), rownames(essay.data)))
    essay.overlap.data <- essay.data[essay.overlap.ids,]
    out[[paste0(data.type, '.data')]] <- essay.overlap.data
  }
  return(out)
}

combine.overlaps <- function (overlap1, overlap2, overlap1.type, overlap2.type) {
  herv.ids <- intersect(overlap1$pairs$herv.id, overlap2$pairs$herv.id)
  out <- list()
  out$herv.ranges <- overlap1$herv.ranges[herv.ids]
  overlap1.pairs <- overlap1$pairs[overlap1$pairs$herv.id %in% herv.ids, ]
  overlap2.pairs <- overlap2$pairs[overlap2$pairs$herv.id %in% herv.ids, ]
  out$pairs <- data.frame(matrix(nrow = 0, ncol = 3))
  for(herv.id in herv.ids) {
    for(essay1.id in overlap1.pairs[overlap1.pairs$herv.id == herv.id, 2]) {
      for(essay2.id in overlap2.pairs[overlap2.pairs$herv.id == herv.id, 2]) {
        out$pairs <- rbind.data.frame(out$pairs, c(herv.id, essay1.id, essay2.id), stringsAsFactors = FALSE)
      }
    }
  }
  colnames(out$pairs) <- c('herv.id', paste0(overlap1.type, '.id'), paste0(overlap2.type, '.id'))
  
  out[[paste0(overlap1.type, '.ranges')]] <- overlap1[[paste0(overlap1.type, '.ranges')]][unique(overlap1.pairs[, 2])]
  out[[paste0(overlap2.type, '.ranges')]] <- overlap1[[paste0(overlap1.type, '.ranges')]][unique(overlap1.pairs[, 2])]
  out[[paste0(overlap1.type, '.data')]] <- overlap1[[paste0(overlap1.type, '.data')]][unique(overlap1.pairs[, 2]),]
  out[[paste0(overlap2.type, '.data')]] <- overlap1[[paste0(overlap1.type, '.data')]][unique(overlap1.pairs[, 2]),]
  return (out)
}

print.overlap.info <- function(overlap) {
  message(paste0("Overlap info:\n# Overlaps: ", length(overlap$hits), 
                 "\n# hERVs: ", length(overlap$herv.ranges),
                 "\n# probes: ", dim(overlap$essay.data)[1]))
}


enlarge.ranges <- function(ranges, flanking) {
    enlarged.ranges <- GRanges(
        seqnames = seqnames(ranges),
        ranges = IRanges(start = start(ranges) - flanking, end = end(ranges) + flanking),
        strand = strand(ranges),
        name = ranges$name,
        score = ranges$score
    )
    names(enlarged.ranges) <- names(ranges)
    return (enlarged.ranges)
}

name.ranges.by.coordinates <- function(ranges) {
  names <- paste0(as.character(seqnames(ranges)), ':', start(ranges), '-', end(ranges))
  names(ranges) <- names
  return (ranges)
}

if(!file.exists(PATHS$HERV.RANGES.DATA)) {
  hervS1.ranges <- import(PATHS$HERVS1.ANNOT, format = 'BED')
  hervS1.ranges <- name.ranges.by.coordinates(hervS1.ranges)
  hervS2.ranges <- import(PATHS$HERVS2.ANNOT, format = 'BED')
  hervS2.ranges <- name.ranges.by.coordinates(hervS2.ranges)
  hervS3.ranges <- import(PATHS$HERVS3.ANNOT, format = 'BED')
  hervS3.ranges <- name.ranges.by.coordinates(hervS3.ranges)
  
  hervS1.1kb.ranges <- enlarge.ranges(hervS1.ranges, 1000)
  hervS2.1kb.ranges <- enlarge.ranges(hervS2.ranges, 1000)
  hervS3.1kb.ranges <- enlarge.ranges(hervS3.ranges, 1000)
  
  hervS1.2kb.ranges <- enlarge.ranges(hervS1.ranges, 2000)
  hervS2.2kb.ranges <- enlarge.ranges(hervS2.ranges, 2000)
  hervS3.2kb.ranges <- enlarge.ranges(hervS3.ranges, 2000)
  
  save(
    hervS1.ranges,
    hervS2.ranges,
    hervS3.ranges,
    hervS1.1kb.ranges,
    hervS2.1kb.ranges,
    hervS3.1kb.ranges,
    hervS1.2kb.ranges,
    hervS2.2kb.ranges,
    hervS3.2kb.ranges,
    file = PATHS$HERV.RANGES.DATA
  )
  save(hervS2.ranges, file = PATHS$HERVS2.RANGES.DATA)
  save(hervS2.2kb.ranges, file = PATHS$HERVS2.2KB.RANGES.DATA)
} else {
  load(PATHS$HERV.RANGES.DATA)
}

load(PATHS$EXPR.DATA)
expr.data <- f4.norm

load(PATHS$METH.DATA)
meth.data <- data.frame(beta)
rm(beta)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    # expr.overlap.name <- paste0('herv', set, flanking, '.expr.overlap')
    # assign(expr.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), expr.ranges, expr.data, 'expr')) 
    # meth.overlap.name <- paste0('herv', set, flanking, '.meth.overlap')
    # assign(meth.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), meth.ranges, meth.data, 'meth'))
    # snp.overlap.name <- paste0('herv', set, flanking, '.snp.overlap')
    # assign(snp.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), snp.ranges, NULL, 'snp'))
    tfbs.overlap.name <- paste0('herv', set, flanking, '.tfbs.overlap')
    assign(tfbs.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), blood.tfbs.ranges, NULL, 'tfbs'))
  }
}

save(
  hervS1.expr.overlap,
  hervS2.expr.overlap,
  hervS3.expr.overlap,
  hervS1.1kb.expr.overlap,
  hervS2.1kb.expr.overlap,
  hervS3.1kb.expr.overlap,
  hervS1.2kb.expr.overlap,
  hervS2.2kb.expr.overlap,
  hervS3.2kb.expr.overlap,
  file = PATHS$HERV.EXPR.OVERLAP.DATA
)

save(
  hervS1.meth.overlap,
  hervS2.meth.overlap,
  hervS3.meth.overlap,
  hervS1.1kb.meth.overlap,
  hervS2.1kb.meth.overlap,
  hervS3.1kb.meth.overlap,
  hervS1.2kb.meth.overlap,
  hervS2.2kb.meth.overlap,
  hervS3.2kb.meth.overlap,
  file = PATHS$HERV.METH.OVERLAP.DATA
)

save(
  hervS1.snp.overlap,
  hervS2.snp.overlap,
  hervS3.snp.overlap,
  hervS1.1kb.snp.overlap,
  hervS2.1kb.snp.overlap,
  hervS3.1kb.snp.overlap,
  hervS1.2kb.snp.overlap,
  hervS2.2kb.snp.overlap,
  hervS3.2kb.snp.overlap,
  file = PATHS$HERV.SNP.OVERLAP.DATA
)

save(
  hervS1.tfbs.overlap,
  hervS2.tfbs.overlap,
  hervS3.tfbs.overlap,
  hervS1.1kb.tfbs.overlap,
  hervS2.1kb.tfbs.overlap,
  hervS3.1kb.tfbs.overlap,
  hervS1.2kb.tfbs.overlap,
  hervS2.2kb.tfbs.overlap,
  hervS3.2kb.tfbs.overlap,
  file = PATHS$HERV.TFBS.OVERLAP.DATA
)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    overlap.name <- paste0('herv', set, flanking, '.both.overlap') 
    cat(overlap.name, fill = TRUE)
    assign(overlap.name, combine.overlaps(get(paste0('herv', set, flanking, '.expr.overlap')),
                                          get(paste0('herv', set, flanking, '.meth.overlap')),
                                          'expr', 'meth'))
  }
}
