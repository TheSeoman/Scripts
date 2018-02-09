source('Scripts/R/paths.R')
source('Scripts/R/prepare.ggm.data.R')

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 25

prepare.ggm.data(set = 'hervS1', filter = 'snp', seed = 'meqtl', string = F, snp.count.threshold = 5, flanking = 5e5, batch, batch.size)
