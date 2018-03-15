source('Scripts/R/paths.R')

load(PATHS$MEQTL.COSMO.DATA)
load(PATHS$MEQTL.TRANS.DATA)
load(PATHS$HERV.MEQTL.OVERLAP.DATA)
load(PATHS$HERV.MEQTL.TRANS.OVERLAP.DATA)
load(PATHS$HERV.MEQTL.OVERLAP.DATA)


cis.cosmo <- cosmo[cosmo$snp.chr == cosmo$cpg.chr, ]

total.snps <- unique(as.character(cosmo$snp))
total.cpgs <- unique(as.character(cosmo$cpg))

trans.snps <- unique(as.character(trans.cosmo$snp))
trans.cpgs <- unique(as.character(trans.cosmo$cpg))

meqtl.snp.counts <- table(cosmo$snp)
meqtl.snp.counts <- meqtl.snp.counts[meqtl.snp.counts > 0]
meqtl.snp.counts <- data.frame(meqtl.snp.counts)

meqtl.meth.counts <- table(cosmo$cpg)
meqtl.meth.counts <- meqtl.meth.counts[meqtl.meth.counts > 0]
meqtl.meth.counts <- data.frame(meqtl.meth.counts)


meqtl.snp.hist <-  ggplot(meqtl.snp.counts, aes(Freq)) + geom_histogram(breaks = seq(0, 350, by = 3)) + 
  theme(text = element_text(size=10)) + labs(title='B', x = 'Proportion of NAs', y = 'Count') + scale_y_log10()

hervS2.meqtl.either.pairs <- hervS2.meqtl.overlap$either
hervS2.meqtl.either.cpgs <- as.character(unique(hervS2.meqtl.either.pairs$cpg))
hervS2.meqtl.either.snps <- as.character(unique(hervS2.meqtl.either.pairs$snp))

hervS2.meqtl.meth.pairs <- hervS2.meqtl.overlap$meth
hervS2.meqtl.meth.cpgs <- as.character(unique(hervS2.meqtl.meth.pairs$cpg))
hervS2.meqtl.meth.snps <- as.character(unique(hervS2.meqtl.meth.pairs$snp))

hervS2.meqtl.snp.pairs <- hervS2.meqtl.overlap$snp
hervS2.meqtl.snp.cpgs <- as.character(unique(hervS2.meqtl.snp.pairs$cpg))
hervS2.meqtl.snp.snps <- as.character(unique(hervS2.meqtl.snp.pairs$snp))

top.meth <- which.max(table(as.character(hervS2.meqtl.overlap$meth$cpg)))

hervS2.meqtl.cpgs <- as.character(unique(hervS2.meqtl.overlap$meth$cpg))
hervS2.meqtl.snps <- as.character(unique(hervS2.meqtl.overlap$snp$snp))

hervS2.meqtl.trans.either.pairs <- hervS2.meqtl.trans.overlap$either
hervS2.meqtl.trans.either.cpgs <- as.character(unique(hervS2.meqtl.trans.either.pairs$cpg))
hervS2.meqtl.trans.either.snps <- as.character(unique(hervS2.meqtl.trans.either.pairs$snp))

hervS2.meqtl.trans.meth.pairs <- hervS2.meqtl.trans.overlap$meth
hervS2.meqtl.trans.meth.cpgs <- as.character(unique(hervS2.meqtl.trans.meth.pairs$cpg))
hervS2.meqtl.trans.meth.snps <- as.character(unique(hervS2.meqtl.trans.meth.pairs$snp))

hervS2.meqtl.trans.snp.pairs <- hervS2.meqtl.trans.overlap$snp
hervS2.meqtl.trans.snp.cpgs <- as.character(unique(hervS2.meqtl.trans.snp.pairs$cpg))
hervS2.meqtl.trans.snp.snps <- as.character(unique(hervS2.meqtl.trans.snp.pairs$snp))
