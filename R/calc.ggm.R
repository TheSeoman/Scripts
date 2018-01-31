source('Scripts/R/paths.R')

require(BDgraph)



set <- 'hervS1'
filter <- 'snp'
seed <- 'meqtl'

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '/')

args <- commandArgs(TRUE)
snp.index <- as.integer(args[1])

load(paste0(GGM.DIR, 'snps.RData'))
snp <- snps[snp.index]


load(paste0(GGM.DIR, 'data/', snp, '.RData'))

iter <- 50000
burnin <- iter/2
cores <- 4

cat(paste0('Calculating ggm for ', snp), fill = T)

ggm <- bdgraph(ggm.data, n = dim(ggm.data)[1], method = "gcgm", algorithm = "bdmcmc", iter = iter, 
                burnin = burnin, g.start = "empty", prior.df = 3, g.prior = 0.5,
                multi.update = NULL, save.all = T, cores = cores )

save(ggm, file = paste0(GGM.DIR, 'ggm/', snp, '.RData'))
