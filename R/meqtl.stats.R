source('Scripts/R/paths.R')

load(PATHS$MEQTL.COSMO.DATA)
load(PATHS$MEQTL.TRANS.DATA)

cis.cosmo <- cosmo[cosmo$snp.chr == cosmo$cpg.chr, ]

total.snps <- unique(as.character(cosmo$snp))
total.cpgs <- unique(as.character(cosmo$cpg))

trans.snps <- unique(as.character(trans.cosmo$snp))
trans.cpgs <- unique(as.character(trans.cosmo$cpg))
