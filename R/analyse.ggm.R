source('Scripts/R/paths.R')

require('illuminaHumanv3.db')

set <- 'hervS1'
filter <- 'snp'
seed <- 'meqtl'

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '/')
load(paste0(GGM.DIR, 'snps.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))

load(PATHS$HERV.EQTL.OVERLAP.DATA)
cis.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('cis.', filter)]]
trans.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('trans.', filter)]]

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
