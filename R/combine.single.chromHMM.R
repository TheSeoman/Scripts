source('Scripts/R/paths.R')

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

meth.ranges.annotation <- data.frame(meth.ranges.annotation.list)
