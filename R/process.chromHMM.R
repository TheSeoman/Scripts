source('Scripts/R/paths.R')

require(GenomicRanges)

sum.up.chromHMM.states = function (ann.list, ranges) {
  
  message(paste0('# areas: ', length(ann.list), ' #elements: ', length(unlist(ann.list))))
  for (element in c(1:length(ann.list))) {
    element.name = ls(ann.list[i])
    element.split = strsplit(element.name, split = '[:-]')
    ann.list[[i]][[1]][2] = element.split[[1]][2]
    ann.list[[i]][[length(ann.list[[i]])]][3] = element.split[[1]][3]
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
hervS2.elements = ls(hervS2.annotation[[1]])

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
  return(sample.states[hervS3.elements])
})
save(hervS3.annotation, file = PATHS$HERVS3.CHROMMHMM.DATA)


hervS3.summed.annotation = sum.up.chromHMM.states(hervS3.annotation, hervS3.ranges)

hervS3.ranges.withann = hervS3.ranges[hervS3.elements %in% ls(hervS3.annotation[[1]])]
sp <- strsplit(ls(hervS3.annotation[[1]]), split = '[-:]')

