DATA.DIR = "/home/icb/julian.schmidt/data/"

HERV.DATA <- paste0(DATA.DIR, 'herv/ranges.RData')

METH.OVERLAP.DATA <- paste0(DATA.DIR, 'overlaps/methylation.RData')
roadmap.dir = "/storage/groups/groups_epigenereg/datasets/roadmap/"
out.dir = paste0(DATA.DIR, "chromHMM/")

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

chromHMM.range.annotation <- function (ranges, ids, dir = paste0(roadmap.dir, "chromHMM/15state/"), suffix = "_15_coreMarks_mnemonics.bed.bgz") {
  library(Rsamtools)
  library(GenomicRanges)
  annotation <- sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=ranges[avail])
    ann.split = sapply(avail.ann, function(x) strsplit(x, "\t") )  
    
    message(paste0('Processing ', id))
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

load(HERV.DATA)
message('Calculate chromHMM states for hervS1')
hervS1.annotation.proportions <- chromHMM.range.annotation(hervS1.ranges, ids)
message('Calculate chromHMM states for hervS2')
hervS2.annotation.proportions <- chromHMM.range.annotation(hervS2.ranges, ids)
message('Calculate chromHMM states for hervS3')
hervS3.annotation.proportions <- chromHMM.range.annotation(hervS3.ranges, ids)

save(hervS1.annotation.proportions, hervS2.annotation.proportions, hervS3.annotation.proportions, file = paste0(out.dir, "herv.proportions.RData"))


