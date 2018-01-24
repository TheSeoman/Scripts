source('Scripts/R/paths.R')

require(GenomicRanges)

combine.single.1nt.chromHMM <-
  function (sample.dir, samples, ranges, start = 1, end = 27) {
    annotation.list <- lapply(list.files(sample.dir)[start:end], function(file) {
      message(paste0('loading: ', file))
      load(paste0(sample.dir, file))
      temp <- unlist(lapply(annotation, function (ann) ann[4]))
      return(temp)
    })
    ranges.map <- names(ranges)  
    names(ranges.map) <- paste0(seqnames(ranges), ':', start(ranges), '-', end(ranges))
    
    samples <- samples[order(samples)]
    
    ann.list <- lapply(annotation.list, function(annotation) {
      annotation.loc <- names(annotation)
      message('...')
      missing.loc <- ranges.map[!(names(ranges.map) %in% annotation.loc)]
      names(annotation) <- ranges.map[annotation.loc]
      annotation[missing.loc] <- NA
      annotation <- annotation[order(names(annotation))]
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

combine.single.broad.chromHMM <- function(sample.dir, samples, ranges, start = 1, end = 27) {
  # ranges.map <- names(ranges)
  # names(ranges.map) <- paste0(seqnames(ranges), ':', start(ranges), '-', end(ranges))
  
  annotation.list <- lapply(list.files(sample.dir)[start:end], function(file) {
    message(paste0('loading: ', sample.dir, file))
    load(paste0(sample.dir, file))
    # annotation <- annotation[!duplicated(names(annotation))]
    # ann.loc <- names(annotation)
    # for (loc in ann.loc) {
    #   split <- unlist(strsplit(loc, ':|-'))
    #   annotation[loc][[1]][[1]][2] <- split[2]
    #   annotation[loc][[length(annotation[loc])]][[1]][3] <- split[3]
    # }
    # names(annotation) <- ranges.map[ann.loc]
    return(annotation)
  })
}

get.majority.chromHMM <- function (annotation.list, samples, ranges) {
  majority.annotation.list <- lapply(annotation.list, function(sample) {
    temp <- lapply(sample, function(ann) {
      if (length(ann) == 1) {
        return(ann[[1]][4])
      } else {
        lengths <- unlist(lapply(ann, function(x) as.numeric(x[3])-as.numeric(x[2])))
        return(ann[[which(lengths==max(lengths)[1])]][4])
      }
    })
    return(temp)
  })
  
  complete.list <- lapply(majority.annotation.list, function(annotation) {
    message('...')
    annotation.ids <- names(annotation)
    missing.ids <- names(ranges)[!names(ranges) %in% annotation.ids]
    annotation[missing.ids] <- NA
    annotation <- annotation[order(names(annotation))]
    return(annotation)
  })
  
  message('Creating matrix...')
  annotation.combined <-
    data.frame(matrix(nrow = length(ranges), ncol = length(complete.list)))
  
  rownames(annotation.combined) <- names(ranges)
  colnames(annotation.combined) <- samples[start:end]
  
  for (i in 1:length(complete.list)) {
    annotation.combined[, i] <- unlist(complete.list[[i]])
  }
  
  return(annotation.combined)
}



load(PATHS$CHROMHMM.SAMPLE.DATA)
# meth
# load(PATHS$METH.RANGES.DATA)
# sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/meth.ranges/')
# meth.chromhmm.states <- combine.single.1nt.chromHMM(sample.dir, ids, meth.ranges)
# save(meth.chromhmm.states, file = PATHS$METH.CHROMHMM.DATA)

# snp
# load(PATHS$SNP.RANGES.DATA)
# sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/snp.ranges/')
# snp.chromhmm.states <- combine.single.1nt.chromHMM(sample.dir, ids, snp.ranges, 1, 27)
# save(snp.chromhmm.states, file = PATHS$SNP.CHROMHMM.DATA)

# expr
# load(PATHS$EXPR.RANGES.DATA)
# sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/expr.ranges/')
# expr.chromhmm.annotation.list <- combine.single.broad.chromHMM(sample.dir, ids, expr.ranges)
# save(expr.chromhmm.annotation.list, file = PATHS$EXPR.CHROMHMM.LIST.DATA)
# expr.chromhmm.annotation <- get.majority.chromHMM(expr.chromhmm.annotation.list, ids, expr.ranges)
# save(expr.chromhmm.annotation, file = PATHS$EXPR.CHROMHMM.DATA)

# hervS2
load(PATHS$HERVS2.RANGES.DATA)
sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/S2/')
hervS2.annotation.list <- combine.single.broad.chromHMM(sample.dir, ids, hervS2.ranges)
save(hervS2.annotation.list, file = PATHS$HERVS2.CHROMHMM.DATA)