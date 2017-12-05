source('Scripts/R/paths.R')

require(GenomicRanges)

## annotate with the roadmap chromHMM states
roadmap.samples = read.csv(paste0(PATHS$ROADMAP.DIR, "sample_info.txt"),
                           sep="\t",
                           stringsAsFactors=F)

mnemonics =
  read.csv(paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/mnemonics.txt"),
           sep="\t")
levels = with(mnemonics, paste(STATE.NO., MNEMONIC, sep="_"))
labels = mnemonics[,"DESCRIPTION"]

use = roadmap.samples$ANATOMY == "BLOOD"
ids = roadmap.samples[use,"Epigenome.ID..EID."]
cells = roadmap.samples[use,"Standardized.Epigenome.name"]

## match this to the houseman cell types
type = rep("other", sum(use))
type[grep("CD4", roadmap.samples[use,5])] = "CD4T"
type[grep("CD8", roadmap.samples[use,5])] = "CD8T"
type[grep("CD15", roadmap.samples[use,5])] = "Gran"
type[grep("CD56", roadmap.samples[use,5])] = "NK"
type[grep("CD19", roadmap.samples[use,5])] = "Bcell"
type[grep("CD14", roadmap.samples[use,5])] = "Mono"


#call for each sample
sum.up.chromHMM.states = function (ann.list) {
  
  res = lapply(ann.list, function(sample) {
  message(paste0('# areas: ', length(sample), ' #elements: ', length(unlist(sample))))
  for (i in c(1:length(sample))) {
    element.name = ls(sample[i])
    element.split = strsplit(element.name, split = '[:-]')
    sample[[i]][[1]][2] = element.split[[1]][2]
    sample[[i]][[length(sample[[i]])]][3] = element.split[[1]][3]
  }
  
  temp <- unlist(sample)
  
  ann.table = data.frame(matrix(temp[!is.na(temp)], ncol=4, byrow=TRUE), stringsAsFactors=FALSE) 
  colnames(ann.table) = c("chr", "start", "end", "state")
  ann.table$length = as.numeric(ann.table$end) - as.numeric(ann.table$start)
  
  ann.summed = data.frame(matrix(ncol = length(levels), nrow = 0))
  colnames(ann.summed) = levels
  for(level in levels) {
    ann.summed[1, level] <- sum(ann.table$length[ann.table$state == level]) 
  }
  
  return (ann.summed)
  })
  
  annotation <- data.frame(matrix(unlist(res), ncol = 15, byrow = TRUE))
  colnames(annotation) <- levels
  rownames(annotation) <- ids
  
  return(annotation)
}

load(PATHS$HERVS2.CHROMHMM.DATA)
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
save(hervS1.annotation, file = PATHS$HERVS1.CHROMHMM.DATA)

hervS3.elements = sapply(hervS3.ranges, function(range) {
  irange = ranges(range)
  element = paste0(seqnames(range), ':', start(irange), '-', end(irange))
  return(element)
})

hervS3.annotation <- lapply(hervS2.annotation, function (sample.states) {
  return(sample.states[hervS3.elements])
})
save(hervS3.annotation, file = PATHS$HERVS3.CHROMHMM.DATA)

load(PATHS$HERVS1.CHROMHMM.DATA)
hervS1.summed.annotation = sum.up.chromHMM.states(hervS1.annotation)
hervS3.summed.annotation = sum.up.chromHMM.states(hervS3.annotation)

#plots etc
library(ggplot2)
library(reshape2)

plot.summed.annotation <- function (annotation) {
annotation.prop <- annotation/rowSums(annotation)
plot.data = melt(annotation.prop)
colnames(plot.data) = c("State", "Proportion")
plot.data = cbind(plot.data, Cell=ids)
plot.data = plot.data[order(plot.data$State), ]
plot.data = cbind(plot.data, Type=type)

print(qplot(State, Cell, fill=Proportion, geom="tile", data=plot.data)
      +  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5))
      + facet_grid(Type ~ ., scale="free_y", space="free_y"))
}

plot.weighted.summed.annotation <- function (annotation) {
  houseman = read.csv(PATHS$F.HOUSEMAN, sep = ";")
  colnames(annotation) <- labels
  by.type = apply(annotation / rowSums(annotation), 2, tapply, type, mean)
  pop.mean = colMeans(houseman[, -1])
  overall.state = t(pop.mean %*% by.type[names(pop.mean), ])
  colnames(overall.state) = "Weighted.proportion"
  overall.state = data.frame(State = rownames(overall.state), overall.state)
  ggplot(data=overall.state) + geom_bar(aes(State, y=Weighted.proportion), stat="identity")  + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=13), axis.title=element_text(size=14), plot.title=element_text(size=16)) + ggtitle('chromHMM state proportions over all HERV set 1 elements')
}

png(paste0(PATHS$PLOT.DIR, 'hervS1.chromHMM.weighted.png'), width = 700, height = 700)
plot.weighted.summed.annotation(hervS1.summed.annotation)
dev.off()

