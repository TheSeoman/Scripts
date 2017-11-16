
DATA.DIR = "/media/data/Masterarbeit/data/"

HERV.DATA <- paste0(DATA.DIR, 'herv/ranges.RData')

METH.OVERLAP.DATA <- paste0(DATA.DIR, 'overlaps/methylation.RData')

roadmap.dir = paste0(DATA.DIR, "roadmap/") #"/storage/groups/groups_epigenreg/roadmap/"
out.dir = paste0(DATA.DIR, "chromHMM/")

library(ggplot2)
library(reshape2)

## annotate with the roadmap chromHMM states
roadmap.samples = read.csv(paste0(roadmap.dir, "sample_info.txt"),
                           sep="\t",
                           stringsAsFactors=F)

mnemonics =
  read.csv(paste0(roadmap.dir, "chromHMM/15state/mnemonics.txt"),
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


load(HERV.DATA)
load(METH.OVERLAP.DATA)



herv.annotation.proportions <- chromHMM.range.annotation.proportions(hervS3.ranges, ids)
plot.data = melt(test3)
colnames(plot.data) = c("State", "Cell", "Proportion")
plot.data = plot.data[order(plot.data$State), ]
plot.data = cbind(plot.data, type=type)

print(qplot(State, Cell, fill=Proportion, geom="tile", data=plot.data)
      +  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5))
      + facet_grid(type ~ ., scale="free_y", space="free_y"))


herv.cpg.annotation.proportions <- chromHMM.annotation(meth.S3.overlap$herv.ranges, ids)
cpg.herv.annotation <- chromHMM.annotation(meth.S3.overlap$essay.ranges, ids)

write.table(annotation,
            file=file.path(out.dir, "cpgs_with_epigenetic_states.txt"),
            sep="\t",
            quote=F)

herv.annotation.counts = apply(herv.annotation, 2, function(ann)
  table(factor(ann, levels=levels, labels=labels)))

cpg.herv.annoatation.counts = apply(cpg.herv.annotation, 2, function(ann)
  table(factor(ann, levels=levels, labels=labels)))

plot.data = t(annotation.counts / length(cpg.ranges))

plot.data = melt(plot.data)
colnames(plot.data) = c("Cell", "State", "Proportion")
plot.data = cbind(plot.data, type=type)


pdf(file=file.path(out.dir, "epigenetic_state_annotation.pdf"))
print(qplot(State, Cell, fill=Proportion, geom="tile", data=plot.data)
      +  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5))
      + facet_grid(type ~ ., scale="free_y", space="free_y"))
dev.off()

gran = plot.data[plot.data$type == "Gran",c("State", "Proportion")]
gran = gran[order(gran$State),]

write.table(gran,
            file=file.path(out.dir,
                           
                           "epigenetic_state_annotation_granulocytes.txt"),
            sep="\t",
            quote=F,
            row.names=F)

pdf(file=file.path(out.dir,
                   "epigenetic_state_annotation_granulocytes.pdf"))
ggplot(data=gran) +
  geom_bar(aes(State, Proportion), stat="identity") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5)) +
  labs(title="CpG annotations in Neutrophils (Granulozytes)")
dev.off()

## check the average cell composition in our samples from the houseman
estimates
houseman =
  read.csv(paste0(DATA.DIR, "Houseman/KF4_QN_estimated_cell_distribution_meanimpute473_lessThanOneTRUE.csv"),
           sep=";")

pdf(file=file.path(out.dir, "houseman_boxplot_kora.pdf"))
boxplot.matrix(as.matrix(houseman[,-1]),
               xlab="Cell type",
               ylab="Estimated fraction")
dev.off()

## weight the percentages by the cell type estimates
by.type = apply(annotation.counts / length(cpg.ranges),
                1,
                tapply,
                type,
                mean)
pop.mean = colMeans(houseman[,-1])

overall.state = t(pop.mean %*% by.type[names(pop.mean),])
colnames(overall.state) = "Weighted.proportion"
overall.state = data.frame(State=rownames(overall.state), overall.state)

write.table(overall.state,
            file=file.path(out.dir,
                           "epigenetic_state_annotation_weighted.txt"),
            sep="\t",
            quote=F,
            row.names=F)

pdf(file=file.path(out.dir, "epigenetic_state_annotation_weighted.pdf"))
ggplot(data=overall.state) +
  geom_bar(aes(State, y=Weighted.proportion), stat="identity") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust=0.5))
dev.off()

## also make a probabilistic annotation (weighted by cell types)
weighted = t(apply(annotation, 1, function(ann) {
  ## make a cell x state indicator matrix
  tab = sapply(levels, function(l) as.numeric(ann == l))
  ## average by cell type
  by.type = apply(tab, 2, tapply, type, mean)
  ## weight by population means
  weighted = pop.mean %*% by.type[names(pop.mean),]
  return(weighted)
}))
colnames(weighted) = labels

write.table(weighted,
            file=file.path(out.dir,
                           "cpgs_with_epigenetic_states_weighted.txt"),
            sep="\t",
            quote=F)



-----------------------
  #' Annotate positions with their epigenetic states
  #'
  #' @param cpg.ranges GRanges object with the positions to annotate
#' @param ids character vector with the roadmap epigenome ids to use
#' @param dir the directory where chromHMM files are stored
#'
#' @return character matrix with length(cpg.ranges) rows and length(ids) columns
#'         with the state annotation of the range in each epigenome
#'
#' @export
chromHMM.annotation <- function(cpg.ranges, 
                                ids, 
                                
                                dir=paste0(roadmap.dir, "chromHMM/15state/"), 
                                suffix="_15_coreMarks_mnemonics.bed.bgz")
{
  library(Rsamtools)
  
  annotation = sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(cpg.ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=cpg.ranges[avail])
    avail.ann = sapply(avail.ann, function(x) strsplit(x, "\t")[[1]][4])
    ann = rep(NA, length(cpg.ranges))
    ann[avail] = avail.ann
    return(ann)
  })
  colnames(annotation) = ids
  rownames(annotation) = names(cpg.ranges)
  
  return(annotation)
}

chromHMM.range.annotation.proportions <- function (ranges, ids, dir = paste0(roadmap.dir, "chromHMM/15state/"), suffix = "_15_coreMarks_mnemonics.bed.bgz") {
  library(Rsamtools)
  library(GenomicRanges)
  library(sqldf)
  annotation <- sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=ranges[avail])
    ann.split = sapply(avail.ann, function(x) strsplit(x, "\t") )  

    message(paste0('# areas: ', length(ann.split), ' #elements: ', length(unlist(ann.split))))
    for (i in c(1:length(ann.split))) {
      ann.split[[i]][1][[1]][2] = start(ranges[i])
      ann.split[[i]][length(ann.split[[1]])][[1]][3] = end(ranges[i])
    }
  
  
    ann.table = data.frame(matrix(unlist(ann.split), ncol=4, byrow=TRUE),stringsAsFactors=FALSE) 
    colnames(ann.table) = c("chr", "start", "end", "state")
    ann.table$length = as.numeric(ann.table$end) - as.numeric(ann.table$start)
    
    ann.summed = data.frame(matrix(ncol = 15, nrow = 0))
    colnames(ann.summed) = levels
    for(level in levels) {
      ann.summed[1, level] <- sum(ann.table$length[ann.table$state == level]) 
    }
    
    return (ann.summed)
  }) 

  annotation <- data.frame(matrix(unlist(annotation, nrow = 15)))
  rownames(annotation) <- levels
  colnames(annotation) <- ids

  return(annotation)
}






