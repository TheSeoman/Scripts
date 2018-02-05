source('Scripts/R/paths.R')

require('illuminaHumanv3.db')

set <- 'hervS1'
filter <- 'snp'
seed <- 'meqtl'
flanking <- 2.5e5
string  <- F

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
load(paste0(GGM.DIR, 'snps.RData'))
load(paste0(GGM.DIR, 'data.meta.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))

load(PATHS$HERV.EQTL.OVERLAP.DATA)

load(PATHS$EQTM.ME.DATA)
cis.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('cis.', filter)]]
trans.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('trans.', filter)]]

trans.eqtm.pairs <- eqtm.me$trans$eqtls[, c('snps', 'gene')]
cis.eqtm.pairs <- eqtm.me$cis$eqtls[, c('snps', 'gene')]
eqtm.pairs <- rbind.data.frame(trans.eqtm.pairs, cis.eqtm.pairs, stringsAsFactors = F)

snp <- snps[1]

genes <- unlist(as.list(illuminaHumanv3SYMBOL))
genes <- genes[!is.na(genes)]


ggms <- list()
ggm.overview <- matrix(nrow = length(snps), ncol = 6)
rownames(ggm.overview) <- snps
colnames(ggm.overview) <- c('entities', 'edges', 'meQTLs', 'meQTL edges', 'eQTLs', 'eQTL edges')

for (snp in snps) {
  cat(paste0('Processing snp: ', snp), fill = T)
  load(paste0(GGM.DIR, 'ggm/', snp, '.RData'))
  ggms[[snp]] <- ggm
  entities <- data.overview[snp, 'total.entities']
  edges <- sum(ggm$last_graph)/2
  meqtls <- data.overview[snp, 'cpgs']
  meqtl.edges <- sum(ggm$last_graph[1, 2:(meqtls+1)])
  
  eqtl.expr.ids <- unique(c(as.character(cis.eqtl.pairs[cis.eqtl.pairs$snps == snp, 'gene']), as.character(trans.eqtl.pairs[trans.eqtl.pairs$snps == snp, 'gene'])))
  eqtl.expr.genes <- genes[eqtl.expr.ids]
  
  eqtl.indices <- which(colnames(ggm$last_graph) %in% c(eqtl.expr.ids, eqtl.expr.genes))
  eqtls <- length(eqtl.indices)
  eqtl.edges <- sum(ggm$last_graph[snp, eqtl.indices])
  
  ggm.overview[snp,] <- c(entities, edges, meqtls, meqtl.edges, eqtls, eqtl.edges)
}

save(ggm.overview, file = paste0(GGM.DIR, 'ggm.overview.RData'))
