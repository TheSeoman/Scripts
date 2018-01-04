source('Scripts/R/paths.R')

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
    annotation.list <- lapply(list.files(sample.dir)[c(start:end)], function(file) {
        message(paste0('loading: ', file))
        load(paste0(sample.dir, file))
        return(annotation)
      })
    annotation.ranges.list <-
      get.ranges.from.annotation(annotation.list)
    
    samples <- samples[order(samples)]
    
    annotation.combined <-
      data.frame(matrix(nrow = length(ranges), ncol = length(annotation.list)))
    rownames(annotation.combined) <- names(ranges)
    colnames(annotation.combined) <- samples[c(start:end)]
    prev.id <- ''
    for (i in c(1:length(annotation.list))) {
      message(paste0('Processing: ', samples[i]))
      overlap <-
        findOverlaps(ranges, annotation.ranges.list[[i]])
      for (j in c(1:length(overlap))) {
        if (j %% round(length(overlap) / 1000) == 0) {
          message(paste0(samples[i], ': ', round(j / length(overlap) * 1000) / 10, '%'))
        }
        colnames(annotation.combined) <- samples
        for (i in c(1:length(annotation.list))) {
          overlap <- findOverlaps(ranges, annotation.ranges.list[[i]])
          for (j in c(1:length(overlap))) {
            id <- names(ranges)[queryHits(overlap)[j]]
            if (id != prev.id) {
              annotation.combined[id, samples[i]] <-
                annotation.ranges.list[[i]][subjectHits(overlap)[j]]$state
              prev.id <- id
            }
          }
        }
      }
    }
    return(annotation.combined)
  }

# run
# args <- commandArgs()
load(PATHS$CHROMHMM.SAMPLE.DATA)
load(PATHS$METH.RANGES.DATA)
sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/meth.ranges/')

meth.chromhmm.states <- combine.single.1nt.chromHMM(sample.dir, ids, meth.ranges)

save(meth.chromhmm.states, file = PATHS$METH.CHROMHMM.DATA)
