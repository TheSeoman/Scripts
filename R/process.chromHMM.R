source('Scripts/R/paths.R')

require(GenomicRanges)

if (!file.exists(PATHS$CHROMHMM.META.DATA)) {
  roadmap.samples = read.csv(
    paste0(PATHS$ROADMAP.DIR, "sample_info.txt"),
    sep = "\t",
    stringsAsFactors = F
  )
  
  mnemonics =
    read.csv(paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/mnemonics.txt"),
             sep = "\t")
  levels = with(mnemonics, paste(STATE.NO., MNEMONIC, sep = "_"))
  labels = mnemonics[, "DESCRIPTION"]
  
  use = roadmap.samples$ANATOMY == "BLOOD"
  ids = roadmap.samples[use, "Epigenome.ID..EID."]
  cells = roadmap.samples[use, "Standardized.Epigenome.name"]
  
  ## match this to the houseman cell types
  type = rep("other", sum(use))
  type[grep("CD4", roadmap.samples[use, 5])] = "CD4T"
  type[grep("CD8", roadmap.samples[use, 5])] = "CD8T"
  type[grep("CD15", roadmap.samples[use, 5])] = "Gran"
  type[grep("CD56", roadmap.samples[use, 5])] = "NK"
  type[grep("CD19", roadmap.samples[use, 5])] = "Bcell"
  type[grep("CD14", roadmap.samples[use, 5])] = "Mono"
  type <- type[order(ids)]
  ids <- ids[order(ids)]
  save(ids, type, levels, file = PATHS$CHROMHMM.META.DATA)
} else {
  load(PATHS$CHROMHMM.META.DATA)
}

#call for each sample
sum.up.chromHMM.states = function (ann.list) {
  res = lapply(ann.list, function(sample) {
  message(paste0('# areas: ', length(sample), ' #elements: ', length(unlist(sample))))
  ann.table = data.frame(matrix(unlist(sample), ncol=4, byrow=TRUE), stringsAsFactors=FALSE) 
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

weight.chromhmm.states <- function(annotation.table) {
  
}

load(PATHS$SNP.CHROMHMM.DATA)

load(PATHS$HERVS2.CHROMHMM.DATA)
hervS2.elements = names(hervS2.annotation[[1]])

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
  annotation.prop <- annotation / rowSums(annotation)
  plot.data = melt(annotation.prop)
  colnames(plot.data) = c("State", "Proportion")
  plot.data = cbind(plot.data, Cell = ids)
  plot.data = plot.data[order(plot.data$State),]
  plot.data = cbind(plot.data, Type = type)
  
  print(qplot(State, Cell, fill=Proportion, geom="tile", data=plot.data)
      +  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5))
      + facet_grid(Type ~ ., scale="free_y", space="free_y"))
}


plot.weighted.summed.annotation <- function (annotation, title) {
  houseman = read.csv(PATHS$F.HOUSEMAN, sep = ";")
  colnames(annotation) <- labels
  by.type = apply(annotation, 2, tapply, type, mean)
  pop.mean = colMeans(houseman[, -1])
  overall.state = t(pop.mean %*% by.type[names(pop.mean), ])
  colnames(overall.state) = "Weighted.proportion"
  overall.state = data.frame(State = rownames(overall.state), overall.state)
  ggplot(data=overall.state) + geom_bar(aes(State, y=Weighted.proportion), stat="identity")  + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=13), axis.title=element_text(size=14), plot.title=element_text(size=16)) + ggtitle(title)
}

plot.weighted.normalized.annotation <- function (annotation, background) {
  backgroud.prop = background / rowSums(background)
  annotation = annotation / background.prop 
  houseman = read.csv(PATHS$F.HOUSEMAN, sep = ";")
  colnames(annotation) <- labels
  by.type = apply(annotation, 2, tapply, type, mean)
  pop.mean = colMeans(houseman[, -1])
  overall.state = t(pop.mean %*% by.type[names(pop.mean), ])
  colnames(overall.state) = "Weighted.proportion"
  overall.state = data.frame(State = rownames(overall.state), overall.state)
  ggplot(data=overall.state) + geom_bar(aes(State, y=Weighted.proportion), stat="identity")  + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5, size=13), axis.title=element_text(size=14), plot.title=element_text(size=16)) + ggtitle(title)
}


png(paste0(PATHS$PLOT.DIR, 'hervS1.chromHMM.weighted.png'), width = 700, height = 700)
plot.weighted.summed.annotation(hervS1.summed.annotation, 'chromHMM state proportions over all HERV set 1 elements')
dev.off()

load(PATHS$HERV.MEQTL.CHROMHMM.ANNOTATION.DATA)

count.sample.annotation <- function (annotation) {
  counted.annotation <-
    apply(annotation, 2, function(sample) {
      counted = data.frame(matrix(ncol = length(levels), nrow = 0))
      colnames(counted) = levels
      for (level in levels) {
        counted[1, level] <- length(sample[sample == level])
      }
      return (counted)
    })
  counted.annotation <- data.frame(matrix(unlist(counted.annotation), ncol = 15, byrow = TRUE))
  colnames(counted.annotation) <- labels
  rownames(counted.annotation) <- ids
  
  return (counted.annotation)
}

hervS1.meth.annotation <- meth.chromhmm.states[unique(names(meth.S1.overlap$essay.ranges)),]

counted.annotation <- get.counted.annotation(test)
png(paste0(PATHS$PLOT.DIR, 'hervS1.meqtl.snp.chromHMM.weighted.png'), width = 500, height = 500)
plot.weighted.summed.annotation(counted.annotation, 'chromHMM state proportions over meQTL-snps in HERV set 1')
dev.off()

count.snp.annotation <- function (snp.annotation) {
  weighted <- apply(snp.annotation, 1, function(states) {
    row <- vector(mode = 'numeric', length = 15)
    temp <- matrix(nrow = 15, ncol = 27)
    temp[] <- 0
    
    for( i in 1:27 ) {
      temp[which(levels == states[i]) , i] <- 1
    }
    by.type = apply(t(temp), 2, tapply, type, mean)
    pop.mean = colMeans(houseman[, -1])
    pop.mean = pop.mean/sum(pop.mean)
    weighted.state = pop.mean %*% by.type[names(pop.mean), ]
    return(weighted.state)
  })
  rownames(weithed) <- levels
  return(t(weighted))  
  
}

load(PATHS$SNP.CHROMHMM.DATA)

snp.chromhmm.counted <- get.counted.annotation(snp.chromhmm.states)

