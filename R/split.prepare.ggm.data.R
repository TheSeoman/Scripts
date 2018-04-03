source('Scripts/R/paths.R')
source('Scripts/R/prepare.ggm.data.R')
# source('Scripts/R/enrich.ggm.data.string.R')

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 50

prepare.ggm.data.snp.herv(set = 'hervS2', filter = 'snp', seed = 'meqtl', string = T, snp.count.threshold = 5, flanking = 2.5e5, batch = batch, batch.size = batch.size)
# enrich.ggm.data.with.string.cpg.herv(set = 'hervS2', filter = 'meth', seed = 'meqtl', flanking = 2.5e5, batch, batch.size)

