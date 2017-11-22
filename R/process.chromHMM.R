source('Scripts/R/paths.R')

require(GenomicRanges)

sum.up.chromHMM.states = function (ann.list, ranges) {
  
  message(paste0('# areas: ', length(ann.list), ' #elements: ', length(unlist(ann.list))))
  for (element in c(1:length(ann.list))) {
    ann.list[[i]][1][[1]][2] = start(ranges[i])
    ann.list[[i]][length(ann.list[[i]])][[1]][3] = end(ranges[i])
  }
  
  
  ann.table = data.frame(matrix(unlist(ann.list), ncol=4, byrow=TRUE),stringsAsFactors=FALSE) 
  colnames(ann.table) = c("chr", "start", "end", "state")
  ann.table$length = as.numeric(ann.table$end) - as.numeric(ann.table$start)
  
  ann.summed = data.frame(matrix(ncol = 15, nrow = 0))
  colnames(ann.summed) = levels
  for(level in levels) {
    ann.summed[1, level] <- sum(ann.table$length[ann.table$state == level]) 
  }
  
  return (ann.summed)
}

load(PATHS$HERVS2.CHROMMHMM.DATA)
hervS2.elements = ls(hervS2.annotation)

load(PATHS$HERV.DATA)
hervS1.elements = sapply(hervS1.ranges, function(range) {
  irange = ranges(range)
  element = paste0(seqnames(range), ':', start(irange), '-', end(irange))
  return(element)
})

hervS1.annotation <- lapply(hervS2.annotation, function (sample.states) {
  return(sample.states[hervS2.elements %in% hervS1.elements])
})
save(hervS1.annotation, file = PATHS$HERVS1.CHROMMHMM.DATA)

hervS3.elements = sapply(hervS3.ranges, function(range) {
  irange = ranges(range)
  element = paste0(seqnames(range), ':', start(irange), '-', end(irange))
  return(element)
})

hervS3.annotation <- lapply(hervS2.annotation, function (sample.states) {
  return(sample.states[hervS2.elements %in% hervS3.elements])
})
save(hervS3.annotation, file = PATHS$HERVS3.CHROMMHMM.DATA)


