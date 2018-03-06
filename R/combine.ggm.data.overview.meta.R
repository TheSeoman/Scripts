source('Scripts/R/paths.R')

num.batches <- 15
batch.size <- 10

set <- 'hervS2'
filter <- 'meth'
seed <- 'meqtl'
flanking <- 2.5e5
string <- T

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')

if(filter == 'snp') {
  load(paste0(GGM.DIR, 'snps.RData'))
  ids <- snps
} else {
  load(paste0(GGM.DIR, 'cpgs.RData'))
  ids <- cpg.sets
}

load(paste0(GGM.DIR, 'data.meta1.RData'))
load(paste0(GGM.DIR, 'data.overview1.RData'))

total.data.meta <- data.meta
total.data.overview <- data.overview

for(batch in 2:num.batches) {
  load(paste0(GGM.DIR, 'data.meta', batch,'.RData'))
  load(paste0(GGM.DIR, 'data.overview', batch, '.RData'))
  
  range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(ids)))
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

missing.cpg.sets <- rownames(data.overview[data.overview$path.genes == 1, ])
save(missing.cpg.sets, file = paste0(GGM.DIR, 'missing.cpgs.RData'))
for (id in missing.cpg.sets){
  load(paste0(GGM.DIR, 'data/',id, '.RData'))
  colnames(ggm.data)[dim(ggm.data)[2]] <- data.meta[[id]]$path.genes[1]
  save(ggm.data, file = paste0(GGM.DIR, 'data/',id, '.RData'))
}



