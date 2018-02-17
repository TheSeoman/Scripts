source('Scripts/R/paths.R')

num.batches <- 20
batch.size <- 25

set <- 'hervS1'
filter <- 'snp'
seed <- 'meqtl'
flanking <- 2.5e5
string <- T

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
load(paste0(GGM.DIR, 'snps.RData'))
load(paste0(GGM.DIR, 'data.meta1.RData'))
load(paste0(GGM.DIR, 'data.overview1.RData'))

total.data.meta <- data.meta
total.data.overview <- data.overview

for(batch in 2:num.batches) {
  load(paste0(GGM.DIR, 'data.meta', batch,'.RData'))
  load(paste0(GGM.DIR, 'data.overview', batch, '.RData'))
  
  range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(snps)))
  for (i in range) {
    total.data.meta[[i]] <- data.meta[[i]]
  }
  total.data.overview[range,] <- data.overview[range, ]  
}

data.overview <- total.data.overview
data.meta <- total.data.meta

snps <- rownames(data.overview[data.overview$path.genes > 0,])

save(data.overview, file = paste0(GGM.DIR, 'data.overview.RData'))
save(data.meta, file = paste0(GGM.DIR, 'data.meta.RData'))
save(snps, file = paste0(GGM.DIR, 'path.snps.RData'))
