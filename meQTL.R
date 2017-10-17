DATA.DIR = '/media/data/Masterarbeit/data/'
MEQTL.DATA = paste0(DATA.DIR, 'meQTL/cosmopairs_combined_151216.RData')

load(MEQTL.DATA)

cpg <- unique(names(meth.S1.overlap$essay.ranges))

cpg_enriched <- cpg[cpg %in% cosmo$cpg]

pairs_enriched <- cosmo$pair[cosmo$cpg %in% cpg_enriched]

distinct_snps <- unique(cosmo$snp[cosmo$cpg %in% cpg_enriched])

cosmo_herv <- cosmo[cosmo$cpg %in% cpg_enriched]
