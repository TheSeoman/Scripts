DATA.DIR = '/media/data/Masterarbeit/data/'
MEQTL.DATA = paste0(DATA.DIR, 'meQTL/cosmopairs_combined_151216.RData')
METH.OVERLAP.DATA <- paste0(DATA.DIR, 'overlaps/methylation.RData')

find.meQTL.overlap <- function(meth.ranges, meQTL.cosmo) {
  cpg <- unique(names(meth.ranges))
  cpg.enriched <- cpg[cpg %in% meQTL.cosmo$cpg]  
  pairs.enriched <- meQTL.cosmo$pair[meQTL.cosmo$cpg %in% cpg.enriched]
  distinct.snps <- unique(meQTL.cosmo$snp[cosmo$cpg %in% cpg.enriched])
  cosmo_herv <- meQTL.cosmo[meQTL.cosmo$cpg %in% cpg.enriched]
}

load(MEQTL.DATA)
load(METH.OVERLAP.DATA)

hervS1.meQTL.overlap <- find.meQTL.overlap(meth.S1.overlap$essay.ranges, cosmo)









