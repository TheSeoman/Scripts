source('Scripts/R/paths.R')

num.batches <- 38
batch.size <- 50

set <- 'hervS2'
filter <- 'snp'
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
  
  cat(paste0('Processing batch: ', batch, ' of length ', length(data.meta)), fill = T)
  
  range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(ids)))
  for (name in names(data.meta)) {
    total.data.meta[[name]] <- data.meta[[name]]
  }
  total.data.overview[range,] <- data.overview[range, ]  
}

data.overview <- total.data.overview
data.meta <- total.data.meta

for(id in names(data.meta)) {
  data.overview[id, 'total.entities'] <- length(unique(unlist(data.meta[[id]]))) + 1
}

load(paste0(GGM.DIR, 'data.meta.no.comb.RData'))
load(paste0(GGM.DIR, 'data.overview.no.comb.RData'))

best.snps <- lapply(data.meta, function(x) x$best.snps)
duplicated.best.snps <- unique(best.snps[duplicated(best.snps)])
for(duplicated.snp in duplicated.best.snps) {
  sets <- names(best.snps[best.snps == duplicated.snp])
  cat(sets, fill = T)
  new.name <- paste(sets, collapse = '|')
  new.meta <- list()
  for(set in sets) {
    for(name in names(data.meta[[sets[1]]])) {
      new.meta[[name]] <- unique(c(new.meta[[name]], data.meta[[set]][[name]]))
    }
  }
  new.meta$extra.meth.ids <- new.meta$extra.meth.ids[!new.meta$extra.meth.ids %in% new.meta$seed.meth.ids]
  new.meta$extra.meth.genes <- new.meta$extra.meth.genes[!new.meta$extra.meth.genes %in% new.meta$seed.meth.genes]
  new.meta$extra.meth.no.gene.probes <- new.meta$extra.meth.no.gene.probes[!new.meta$extra.meth.no.gene.probes %in% new.meta$seed.meth.no.gene.probes]
  new.meta$extra.tfbs.genes <- new.meta$extra.tfbs.genes[!new.meta$extra.tfbs.genes %in% new.meta$seed.tfbs.genes]
  
  data.meta[sets] <- NULL
  data.meta[[new.name]] <- new.meta
  
  values <- c(length(new.meta$seed.meth.ids), length(new.meta$extra.meth.ids), length(new.meta$seed.tfbs.genes), length(new.meta$extra.meth.ids),
              length(new.meta$seed.meth.genes), length(new.meta$seed.meth.no.gene.probes), length(new.meta$extra.meth.genes), length(new.meta$extra.meth.no.gene.probes))
  names <- c('seed.cpgs', 'extra.cpgs', 'seed.meth.tfs', 'extra.meth.tfs', 'seed.meth.genes', 'seed.meth.no.gene.probes', 'extra.meth.genes', 'extra.meth.no.gene.probes')
  
  for(i in 1:8) {
    data.overview[sets[1], names[i]] <- values[i]
  }
  data.overview <- data.overview[!rownames(data.overview) %in% sets[-1], ]
  rownames(data.overview)[rownames(data.overview) == sets[1]] <- new.name
}

c <- sapply(rownames(data.overview), function(x) {strsplit(x, '\\|')})
cpg.sets <- as.list(c)
names(cpg.sets) <- 1:length(c)
save(cpg.sets, file = paste0(GGM.DIR, 'cpgs.RData')) 

save(data.overview, file = paste0(GGM.DIR, 'data.overview.RData'))
data.meta <- data.meta[rownames(data.overview)]
for(name in names(data.meta)) {
  #data.meta[[name]]$extra.meth.ids <- unique(data.meta[[name]]$extra.meth.ids)
  data.overview[name, 'extra.cpgs'] <- length(data.meta[[name]]$extra.meth.ids)
  data.overview[name, 'total.entities'] <- length(unique(unlist(data.meta[[name]])))
}

save(data.meta, file = paste0(GGM.DIR, 'data.meta.RData'))

missing.cpg.sets <- rownames(data.overview[data.overview$path.genes == 1, ])
save(missing.cpg.sets, file = paste0(GGM.DIR, 'missing.cpgs.RData'))
for (id in missing.cpg.sets){
  load(paste0(GGM.DIR, 'data/',id, '.RData'))
  colnames(ggm.data)[dim(ggm.data)[2]] <- data.meta[[id]]$path.genes[1]
  save(ggm.data, file = paste0(GGM.DIR, 'data/',id, '.RData'))
}



