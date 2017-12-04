source('Scripts/R/paths.R')

## annotate with the roadmap chromHMM states
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

## match this to the houseman cell
type = rep("other", sum(use))
type[grep("CD4", roadmap.samples[use, 5])] = "CD4T"
type[grep("CD8", roadmap.samples[use, 5])] = "CD8T"
type[grep("CD15", roadmap.samples[use, 5])] = "Gran"
type[grep("CD56", roadmap.samples[use, 5])] = "NK"
type[grep("CD19", roadmap.samples[use, 5])] = "Bcell"
type[grep("CD14", roadmap.samples[use, 5])] = "Mono"

chromHMM.snp.annotation <-
  function (snp.ranges,
            id,
            dir = paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/"),
            suffix = "_15_coreMarks_mnemonics.bed.bgz") {
    library(Rsamtools)
    library(GenomicRanges)
    message(paste0('Processing ', id))
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(snp.ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param = snp.ranges[avail])
    avail.ann = sapply(avail.ann, function(x)
      strsplit(x, "\t")[[1]][4])
    ann = rep(NA, length(snp.ranges))
    ann[avail] = avail.ann
    return(ann)
  }

load(PATHS$HERV.SNP.RANGES.DATA)
args <- commandArgs(TRUE)
start <- as.integer(args[1])
end <- as.integer(args[2])
for (i in c(start:end)) {
  id <- ids[i]
  message(paste0('Calculate chromHMM states for snps in hervS3+-2kb in ', id))
  annotation = chromHMM.snp.annotation(hervS2.2kb.snp.ranges, id)
  saveRDS(annotation,
          file = paste0(PATHS$CHROMHMM.OUT.DIR, 'S2.2kb/', id, '.snp.chromHMM.RData'))
}