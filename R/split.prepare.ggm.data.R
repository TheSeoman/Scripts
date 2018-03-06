source('Scripts/R/paths.R')
source('Scripts/R/enrich.ggm.data.string.R')

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 10

enrich.ggm.data.with.string.cpg.herv(set = 'hervS2', filter = 'meth', seed = 'meqtl', flanking = 2.5e5, batch, batch.size)

missing.sets <- list()
for(seed.meth.ids in cpg.sets) {
  set.name <- paste(seed.meth.ids, collapse = '|')
  i <- 1
  if(!file.exists(paste0(GGM.STRING.DIR, 'data/', set.name, '.RData'))) {
    missing.sets[[i]] <- seed.meth.ids
    i <- i + 1
    cat(set.name, fill = T)
  }
}

missing.sets <- m.sets

missing.sets[[1]] <- m.sets[[1]][[1]]
missing.sets[[2]] <- m.sets[[2]][[1]]
missing.sets[[3]] <- m.sets[[3]][[1]]
missing.sets[[4]] <- m.sets[[4]][[1]]
missing.sets[[5]] <- m.sets[[5]][[1]]

save(missing.sets, file = paste0(GGM.STRING.DIR, 'missing.sets.RData'))
