source('Scripts/R/paths.R')

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
                                dir=paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/"), 
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

chromHMM.range.annotation <- function (ranges, ids, dir = paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/"), suffix = "_15_coreMarks_mnemonics.bed.bgz") {
  library(Rsamtools)
  library(GenomicRanges)
  annotation <- lapply(ids, function(id) {
    message(paste0('Processing ', id))
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=ranges[avail])
    ann.split = sapply(avail.ann, function(x) strsplit(x, "\t") )  
    
    return (ann.split)
  }) 
  
  return(annotation)
}

load(HERV.DATA)
message('Calculate chromHMM states for hervS2')
hervS2.annotation <- chromHMM.range.annotation(hervS2.ranges, ids)

save(hervS2.annotation, file = paste0(PATHS$CHROMHMM.OUT.DIR, "hervS2.chromHMM.states.RData"))


