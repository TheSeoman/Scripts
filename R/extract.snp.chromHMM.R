source('Scripts/R/paths.R')

load(PATHS$CHROMHMM.SAMPLE.DATA)
load(PATHS$METH.CHROMHMM.DATA)
load(PATHS$SNP.CHROMHMM.DATA)
load(PATHS$HERV.MEQTL.OVERLAP.DATA)

require(GenomicRanges)

extract.meqtl.annotation <- function (meqtl.overlap, meth.annotation, snp.annotation) {
  #annotation of meqtl-snps, where the snps themselves lie in herv elements
  out$snp.snp <- snp.annotation[unique(meqtl.overlap$snp$snp),]
  #annotation of meqtl-cpgs, where the associated snps lie in herv elements
  out$snp.cpg <- meth.annotation[unique(meqtl.overlap$snp$cpg),]
  out$cpg.snp <- snp.annotation[unique(meqtl.overlap$meth$snp),]
  out$cpg.cpg <- meth.annotation[unique(meqtl.overlap$meth$cpg),]
  out$both.snp <- snp.annotation[unique(meqtl.overlap$both$snp),]
  out$both.cpg <- meth.annotation[unique(meqtl.overlap$both$cpg),]
  out$either.snp <- snp.annotation[unique(meqtl.overlap$either$snp),]
  out$either.cpg <- meth.annotation[unique(meqtl.overlap$either$cpg),]
  return(out)
}

# snps in meQTLs related to hervs

