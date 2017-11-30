source('Scripts/R/paths.R')

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

snp.chromHMM.annotation <- function(snp.ranges, 
                                ids, 
                                dir=paste0(PATHS$ROADMAP.DIR, "chromHMM/15state/"), 
                                suffix="_15_coreMarks_mnemonics.bed.bgz")
{
  library(Rsamtools)
  
  annotation = sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(snp.ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=snp.ranges[avail])
    avail.ann = sapply(avail.ann, function(x) strsplit(x, "\t")[[1]][4])
    ann = rep(NA, length(snp.ranges))
    ann[avail] = avail.ann
    return(ann)
  })
  colnames(annotation) = ids
  rownames(annotation) = names(snp.ranges)
  
  return(annotation)
}

load(PATHS$HERV.SNP.RANGES.DATA)
hervS1.snp.annotation = snp.chromHMM.annotation(hervS1.snp.ranges, ids)

save(hervS1.snp.annotation, file = PATHS$HERVS1.SNP.CHROMHMM.DATA)


