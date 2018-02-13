source('Scripts/R/paths.R')
source('Scripts/R/enrich.ggm.data.string.R')

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 25

enrich.ggm.data.with.string(set = 'hervS1', filter = 'snp', seed = 'meqtl', flanking = 2.5e5, batch, batch.size)
