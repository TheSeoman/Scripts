source('Scripts/R/paths.R')

require('GenomicRanges')

cat('Loading herv-ranges and chromHMM annotations', fill = TRUE)
load(PATHS$HERV.RANGES.DATA)
load(PATHS$HERVS2.2KB.CHROMHMM.DATA)

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

get.annotation.from.ranges <- function(query.ranges, annotation.ranges.list) {
  subject.annotation <- lapply(annotation.ranges.list, function(annotation.ranges) {
    temp <- list()
    overlaps <- findOverlaps(query.ranges, annotation.ranges)
    last.name <- ''
    new.index <- 1
    for(i in 1:length(overlaps)) {
      if ((i %% round(length(overlaps) / 100)) == 0) {
        cat(paste0(round(i/length(overlaps)*100), '% done'), fill = TRUE)
      }
      query.range <- query.ranges[queryHits(overlaps)[i]]
      query.name <- names(query.range)
      ann.range <- annotation.ranges[subjectHits(overlaps)[i]]
      if(query.name != last.name) {
        new.index <- 1
        ann.row <- c(as.character(seqnames(ann.range)), start(query.range), min(end(ann.range), end(query.range)), ann.range$state)
        temp[[query.name]] <- list(ann.row)
        last.name <- query.name
      } else {
        new.index <- new.index + 1
        ann.row <- c(as.character(seqnames(ann.range)), start(ann.range), min(end(ann.range), end(query.range)), ann.range$state)
        temp[[query.name]][[new.index]] <- ann.row
      }
    }    
    
  })
  return(subject.annotation)
}

get.annotation.by.names <- function(query.ranges, annotation.list) {
  names <- paste0(as.character(seqnames(query.ranges)), ':', start(query.ranges), '-', end(query.ranges))
  return (lapply(annotation.list, function(ann) ann[names]))
}


cat('Generating ranges from chromHMM annotations', fill = TRUE)
hervS2.2kb.annotation.ranges.list <- get.ranges.from.annotation(hervS2.2kb.annotation)
cat('Extracting annotations for hervS2.1kb', fill = TRUE)
hervS2.1kb.annotation <- get.annotation.from.ranges(hervS2.1kb.ranges, hervS2.2kb.annotation.ranges.list)
save(hervS2.1kb.annotation, file = paste0(PATHS$DATA.DIR, 'chromHMM/extractions/hervS2.1kb.annotation.RData'))
cat('Extracting annotations for hervS2', fill = TRUE)
hervS2.annotation <- get.annotation.from.ranges(hervS2.ranges, hervS2.2kb.annotation.ranges.list)
save(hervS2.annotation, file = paste0(PATHS$DATA.DIR, 'chromHMM/extractions/hervS2.annotation.RData'))
save(hervS2.annotation, hervS2.1kb.annotation, hervS2.2kb.annotation, file = PATHS$HERVS2.CHROMHMM.DATA)

# cat('Extracting annotations for hervS1', fill = TRUE)
# hervS1.2kb.annotation <- get.annotation.by.names(hervS1.2kb.ranges, hervS2.2kb.annotation)
# hervS1.1kb.annotation <- get.annotation.by.names(hervS1.1kb.ranges, hervS2.1kb.annotation)
# hervS1.annotation <- get.annotation.by.names(hervS1.ranges, hervS2.annotation)
# save(hervS1.annotation, hervS1.1kb.annotation, hervS1.2kb.annotation, file = PATHS$HERVS1.CHROMHMM.DATA)
# 
# cat('Extracting annotations for hervS3', fill = TRUE)
# hervS3.2kb.annotation <- get.annotation.by.names(hervS3.2kb.ranges, hervS1.2kb.annotation)
# hervS3.1kb.annotation <- get.annotation.by.names(hervS3.1kb.ranges, hervS1.1kb.annotation)
# hervS3.annotation <- get.annotation.by.names(hervS3.ranges, hervS1.annotation)
# save(hervS3.annotation, hervS3.1kb.annotation, hervS3.2kb.annotation, file = PATHS$HERVS3.CHROMHMM.DATA)
