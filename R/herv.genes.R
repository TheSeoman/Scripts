source('Scripts/R/paths.R')

library(GenomicFeatures)
require(org.Hs.eg.db)

load(PATHS$HERV.DATA)
hg19.knownGene.db <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene", goldenPath_url="http://hgdownload.soe.ucsc.edu/goldenPath")

knownGene.ranges <- genes(hg19.knownGene.db)

hervS1.gene.overlaps <- findOverlaps(hervS1.ranges, knownGene.ranges, type = 'any')

knownGene.ranges[unique(subjectHits(hervS1.gene.overlaps))]

