DATA.DIR = '/media/data/Masterarbeit/data/'
MEQTL.DATA = paste0(DATA.DIR, 'meQTL/cosmopairs_combined_151216.RData')
METH.OVERLAP.DATA <- paste0(DATA.DIR, 'overlaps/methylation.RData')

MEQTL.HERV.DATA <- paste0(DATA.DIR, 'meQTL/meqtl.herv.RData')
MEQTL.HERV.1KB.DATA <- paste0(DATA.DIR, 'meQTL/meqtl.herv.1kb.RData')
MEQTL.HERV.2KB.DATA <- paste0(DATA.DIR, 'meQTL/meqtl.herv.2kb.RData')

find.meQTL.overlap <- function(meth.ranges, meQTL.cosmo) {
  cpg <- unique(names(meth.ranges))
  cpg.enriched <- cpg[cpg %in% meQTL.cosmo$cpg]  
  pairs.enriched <- meQTL.cosmo$pair[meQTL.cosmo$cpg %in% cpg.enriched]
  distinct.snps <- unique(meQTL.cosmo$snp[cosmo$cpg %in% cpg.enriched])
  cosmo.herv <- meQTL.cosmo[meQTL.cosmo$cpg %in% cpg.enriched,]
  out <- list()
  out$pairs <- pairs.enriched
  out$cpg <- cpg.enriched
  out$snps <- distinct.snps
  out$cosmo <- cosmo.herv
  return(out)
}

print.herv.meqtl.overlap.info <- function(herv.meqtl.overlap){
  message(paste0('Unique relevant cpgs: ', length(herv.meqtl.overlap$cpg)))
  message(paste0('Unique relevant snps: ', length(herv.meqtl.overlap$snps)))
  message(paste0('Relevant pairs: ', length(herv.meqtl.overlap$pairs)))
}

load(MEQTL.DATA)
load(METH.OVERLAP.DATA)

hervS1.meQTL.overlap <- find.meQTL.overlap(meth.S1.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS1.meQTL.overlap)
hervS2.meQTL.overlap <- find.meQTL.overlap(meth.S2.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS2.meQTL.overlap)
hervS3.meQTL.overlap <- find.meQTL.overlap(meth.S3.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS3.meQTL.overlap)
save(hervS1.meQTL.overlap, hervS2.meQTL.overlap, hervS3.meQTL.overlap, file = MEQTL.HERV.DATA)

hervS1.1kb.meQTL.overlap <- find.meQTL.overlap(meth.S1.1kb.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS1.1kb.meQTL.overlap)
hervS2.1kb.meQTL.overlap <- find.meQTL.overlap(meth.S2.1kb.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS2.1kb.meQTL.overlap)
hervS3.1kb.meQTL.overlap <- find.meQTL.overlap(meth.S3.1kb.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS3.1kb.meQTL.overlap)
save(hervS1.1kb.meQTL.overlap, hervS2.1kb.meQTL.overlap, hervS3.1kb.meQTL.overlap, file = MEQTL.HERV.1KB.DATA)

hervS1.2kb.meQTL.overlap <- find.meQTL.overlap(meth.S1.2kb.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS1.2kb.meQTL.overlap)
hervS2.2kb.meQTL.overlap <- find.meQTL.overlap(meth.S2.2kb.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS2.2kb.meQTL.overlap)
hervS3.2kb.meQTL.overlap <- find.meQTL.overlap(meth.S3.2kb.overlap$essay.ranges, cosmo)
print.herv.meqtl.overlap.info(hervS3.2kb.meQTL.overlap)
save(hervS1.2kb.meQTL.overlap, hervS2.2kb.meQTL.overlap, hervS3.2kb.meQTL.overlap, file = MEQTL.HERV.2KB.DATA)

