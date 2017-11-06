load("/media/data/Masterarbeit/data/KF3_mvalue_qn_bmiq.RData")
load("/media/data/Masterarbeit/data/f4.methylation.rdat")

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

load("/media/data/Masterarbeit/data/F4/Expression/f4_annotation.Rdata")

f4.norm.trans = t(f4.norm)
covData <- merge(f4.norm.trans, covars.f4, by.x = 'row.names', by.y = 'ZZ.NR')

exprAnnot <- read.table('/media/data/Masterarbeit/data/HumanHT-12_V3_0_R3_11283641_A.txt', skip = 8, sep = "\t", header = TRUE)

library(sqldf)
x <- sqldf("
  SELECT h.Chr, h.hervID, h.startQuery, h.endQuery, h.strand, h.class, a.ProbeId, a.startPosition, a.endPosition 
  FROM herv h, annotation a
  WHERE a.Chr = h.Chr
  AND ((a.startPosition >= h.startQuery AND a.startPosition <= h.endQuery)
  OR (a.endPosition >= h.startQuery AND a.endPosition <= h.endQuery))
           ")




# Methylation
library(FDb.InfiniumMethylation.hg19)
hm450.hg19 <- getPlatform(platform='HM450', genome='hg19')

hm_overlap <- findOverlaps(hm450.hg19, herv_ranges, type = 'within')
hm_overlap_c <- countOverlaps(herv_ranges, hm450.hg19)


load('/media/data/Masterarbeit/data/F4/KORAF4_illuminamethylation450k_qn_bmiq_n1727/KF4_beta_qn_bmiq.RData')
meth.herv.ids <- names(hm450.hg19)[queryHits(hm_overlap)]
f4.meth.herv <- beta[meth.herv.ids,]
