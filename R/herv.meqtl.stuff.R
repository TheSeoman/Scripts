source('Scripts/R/paths.R')

load(PATHS$HERV.SNP.OVERLAP.DATA)
load(PATHS$HERV.METH.OVERLAP.DATA)

load(PATHS$HERV.MEQTL.TRANS.OVERLAP.DATA)
load(PATHS$MEQTL.TRANS.PAIRS.DATA)


head(hervS1.meth.overlap$pairs)

meth.hervs <- hervS1.meth.overlap$pairs$herv.id
names(meth.hervs) <- hervS1.meth.overlap$pairs$meth.id

snp.hervs <- hervS1.snp.overlap$pairs$herv.id
names(snp.hervs) <- hervS1.snp.overlap$pairs$snp.id

meth.meqtl.pairs <- hervS1.meqtl.trans.overlap$meth[, 1:2]
meth.meqtl.pairs$herv <- meth.hervs[as.character(meth.meqtl.pairs$cpg)]

snp.meqtl.pairs <- hervS1.meqtl.trans.overlap$snp[, 1:2]
snp.meqtl.pairs$herv <- snp.hervs[as.character(snp.meqtl.pairs$snp)]

hervS1.meth.overlap$pairs[hervS1.meth.overlap$pairs$meth.id %in% hervS1.meth.overlap$pairs$meth.id[duplicated(hervS1.meth.overlap$pairs$meth.id)],]

