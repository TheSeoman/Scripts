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

chromHMM.range.annotation <- function (ranges, id, dir = paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/"), suffix = "_15_coreMarks_mnemonics.bed.bgz") {
  library(Rsamtools)
  library(GenomicRanges)
  message(paste0('Processing ', id))
  file = paste0(dir, id, suffix)
  avail = as.logical(seqnames(ranges) %in% seqnamesTabix(file))
  avail.ann = scanTabix(file, param=ranges[avail])
  ann.split = sapply(avail.ann, function(x) strsplit(x, "\t") )  
  return (ann.split)
}

load(HERV.S2.DATA)
# message('Calculate chromHMM states for hervS1')
# hervS1.annotation <- chromHMM.range.annotation(hervS1.ranges, ids)
args <- commandArgs(TRUE)
id <- ids[as.integer(args[1])]
message(paste0('Calculate chromHMM states for hervS2 in ', id))
annotation.name = paste('hervS2', id, 'chromHMM', sep = '.') 
assign(annotation.name, chromHMM.range.annotation(hervS2.ranges, id))
# message('Calculate chromHMM states for hervS3')
# hervS3.annotation <- chromHMM.range.annotation(hervS3.ranges, ids)

annotation <- get(annotation.name)
save(annotation, file = paste0(PATHS$CHROMHMM.OUT.DIR, id, '.chromHMM.RData'))


