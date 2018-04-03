source('Scripts/R/paths.R')

require(BDgraph)

set <- 'hervS2'
filter <- 'snp'
seed <- 'meqtl'
flanking <- 2.5e5
string <- T

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
dir.create(paste0(GGM.DIR, 'ggm/'), showWarnings = F, recursive = T)

args <- commandArgs(TRUE)
index <- as.integer(args[1])

if (filter == 'snp') {
  load(paste0(GGM.DIR, 'extra.snps.RData'))
  id <- snps[index]
} else {
  load(paste0(GGM.DIR, 'cpgs.RData'))
  id <- paste(cpg.sets[[index]], collapse = '|')
}

load(paste0(GGM.DIR, 'data/', id, '.RData'))

iter <- 50000
burnin <- iter/2
cores <- 4

cat(paste0('Calculating ggm for ', id), fill = T)

ggm <- bdgraph(ggm.data, n = dim(ggm.data)[1], method = "gcgm", algorithm = "bdmcmc", iter = iter, 
                burnin = burnin, g.start = "empty", prior.df = 3, g.prior = 0.5,
                multi.update = NULL, save.all = T, cores = cores )

save(ggm, file = paste0(GGM.DIR, 'ggm/', id, '.RData'))
