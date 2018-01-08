source('Scripts/R/paths.R')

library(GenomicRanges)

get.ranges.from.annotation <- function (annotation) {
  annotation.ranges <- lapply(annotation, function(sample) {
    ann.list <- unlist(sample)
    ann.table <-
      data.frame(matrix(ann.list[!is.na(ann.list)], ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
    ann.table <- ann.table[!duplicated(ann.table), ]
    ann.table[, 2] <- as.numeric(ann.table[, 2]) + 1
    ann.table[, 3] <- as.numeric(ann.table[, 3])
    ann.ranges <-
      GRanges(
        seqnames = ann.table[, 1],
        ranges = IRanges(
          start = as.numeric(ann.table[, 2]),
          end = as.numeric(ann.table[, 3])
        ),
        state = ann.table[, 4]
      )
    
    return(ann.ranges)
  })
  return(annotation.ranges)
}

combine.single.1nt.chromHMM <-
  function (sample.dir, samples, ranges, start = 1, end = 27) {
    annotation.list <- lapply(list.files(sample.dir)[start:end], function(file) {
      message(paste0('loading: ', file))
      load(paste0(sample.dir, file))
      temp <- unlist(lapply(annotation, function (ann) ann[[1]][4]))
      return(temp)
    })
    ranges.map <- names(ranges)  
    names(ranges.map) <- paste0(seqnames(ranges), ':', start(ranges), '-', end(ranges))
    
    samples <- samples[order(samples)]
    
    i <- 1
    ann.list <- lapply(annotation.list, function(annotation) {
      message(paste0('processing: ', samples[i]))
      annotation.loc <- names(annotation)
      missing.loc <- ranges.map[!(names(ranges.map) %in% annotation.loc)]
      names(annotation) <- ranges.map[annotation.loc]
      annotation[missing.loc] <- NA
      annotation <- annotation[order(names(annotation))]
      i <- i + 1
      return(annotation)
    })
    
    
    message('Creating matrix...')
    annotation.combined <-
      data.frame(matrix(nrow = length(ranges), ncol = length(annotation.list)))
    
    rownames(annotation.combined) <- names(ranges)[order(names(ranges))]
    colnames(annotation.combined) <- samples[start:end]
    
    for (i in 1:length(ann.list)) {
      annotation.combined[, i] <- ann.list[[i]]
    }
    
    return(annotation.combined)
  }

combine.single.broad.chromHMM <- function(sample.dir, samples, start = 1, end = 27) {
  annotation.list <- lapply(list.files(sample.dir)[start:end], function(file) {
    message(paste0('loading: ', file))
    load(paste0(sample.dir, file))
    ann.loc <- ls(annotation)
    for (loc in ann.loc) {
      split <- unlist(strsplit(loc, ':|-'))
      annotation[loc][[1]][[1]][2] <- split[2]
      annotation[loc][[length(annotation[loc])]][[1]][3] <- split[3]
    }
    return(annotation)
  })
}

load(PATHS$CHROMHMM.SAMPLE.DATA)
# meth
# load(PATHS$METH.RANGES.DATA)
# sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/meth.ranges/')
# meth.chromhmm.states <- combine.single.1nt.chromHMM(sample.dir, ids, meth.ranges)
# save(meth.chromhmm.states, file = PATHS$METH.CHROMHMM.DATA)

# snp
load(PATHS$SNP.RANGES.DATA)
sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/snp.ranges/')
snp.chromhmm.states <- combine.single.1nt.chromHMM(sample.dir, ids, snp.ranges, 1, 1)
save(snp.chromhmm.states, file = PATHS$SNP.CHROMHMM.DATA)

# expr
# sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/expr.ranges/')
# expr.chromhmm.annotation <- combine.single.broad.chromHMM(sample.dir, ids)
# save(expr.chromhmm.annotation, file = PATHS$EXPR.CHROMHMM.DATA)
