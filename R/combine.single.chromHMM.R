source('Scripts/R/paths.R')

get.ranges.from.annotation <- function (annotation) {
  annotation.ranges <- lapply(annotation, function(sample){
    ann.list <- unlist(sample)
    ann.table <- data.frame(matrix(ann.list[!is.na(ann.list)], ncol=4, byrow = TRUE), stringsAsFactors = FALSE) 
    ann.table <- ann.table[!duplicated(ann.table),]
    ann.table[,2] <- as.numeric(ann.table[,2]) + 1
    ann.table[,3] <- as.numeric(ann.table[,3]) 
    ann.ranges <- GRanges(seqnames = ann.table[,1], ranges = IRanges(start = as.numeric(ann.table[,2]), end = as.numeric(ann.table[,3])), state = ann.table[,4])
    
    return(ann.ranges)
  })
  return(annotation.ranges)
}

sample.dir = paste0(PATHS$DATA.DIR, 'chromHMM/S2.2kb/')

hervS2.2kb.annotation = lapply(list.files(sample.dir), function(file) {
  load(paste0(sample.dir, file)) 
  return(annotation)
  })

save(hervS2.2kb.annotation, file = PATHS$HERVS2.2KB.CHROMHMM.DATA)

load(PATHS$METH.RANGES.DATA)
load(PATHS$CHROMHMM.SAMPLE.DATA)
sample.dir <- paste0(PATHS$DATA.DIR, 'chromHMM/meth.ranges/')
meth.ranges.annotation.list <- lapply(list.files(sample.dir), function(file) {
  load(paste0(sample.dir, file))
  return(annotation)
})

meth.ranges.annotation <- data.frame(lapply(meth.ranges.annotation.list, function(sample) {
  return(sapply(sample, function(probe) probe[[1]][4]))
}))

c1 <- names(meth.ranges.annotation[[1]])
c25 <- names(meth.ranges.annotation[[25]])
length(setdiff(c1, c25))

combine.single.1nt.chromHMM(sample.dir, samples, ranges) {
  annotation.list <- lapply(list.files(sample.dir), function(file) {
    load(paste0(sample.dir, file))
    return(annotation)
  })
  annotation.ranges.list <- get.ranges.from.annotation(annotation.list)
  
  samples <- samples[order(samples)]
  
  annotation.combined <- data.frame(matrix(nrow = length(ranges), ncol = length(annotation.list)))
  rownames(annotation.combined) <- names(ranges)
  colnames(annotation.combined) <- samples
  for(i in c(1:length(annotation.list))) {
    overlap <- findOverlaps(ranges, annotation.ranges.list[[i]])
    for(j in c(1:length(overlap))) {
      id <- names(ranges)[queryHits(overlap)[j]]
      if (id != prev.id) {
        annotation.combined[id, samples[i]] <- annotation.ranges.list[[i]][subjectHits(overlap)[j]]$state  
        prev.id <- id
      }
     }
  }
  
  
}