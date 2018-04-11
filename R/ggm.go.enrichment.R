source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')



load(PATHS$EXPR.GENE.ANNOT.DATA)
gene.universe <- unique(probe2gene)

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/hervS2.meqtl.snp.250kb.string/')

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 50


cutoff <- 0.9
load(paste0(GGM.DIR, 'data.meta.RData'))

load(paste0(GGM.DIR, 'filtered.ggm.overview.', cutoff, '.RData'))
load(paste0(GGM.DIR, 'ggm.meta', cutoff, '.RData'))

snp.ggm.go.enrichment <- function(snps) {
  global.enrichments <- list()
  
    for(snp in snps[((batch-1)*batch.size+1):min(length(snps), batch.size*batch)]) {
      cat(paste0(which(snps == snp), ': Processing snp ', snp), fill = T)
      g.meta <- ggm.meta[[snp]]
      d.meta <- data.meta[[snp]]
      ggm.genes <- unique(c(rownames(g.meta$cpg.genes)[g.meta$cpg.genes$con],
                            rownames(g.meta$snp.genes)[g.meta$snp.genes$con],
                            rownames(g.meta$tfs)[g.meta$tfs$con],
                            rownames(g.meta$path.genes)[g.meta$path.genes$con]
      ))
      ggm.genes <- ggm.genes[ggm.genes %in% gene.universe]
      data.genes <- unique(c(d.meta$snp.genes, d.meta$meth.genes, d.meta$tfbs.genes, d.meta$path.genes))
      global.enrichments[[snp]] <- go.enrichment(ggm.genes, gene.universe, gsc, c('BP'))
    }
  save(global.enrichments, file = paste0(GGM.DIR, 'ggm.bp.enrichment.', batch, '.RData'))
}

snp.ggm.go.enrichment(rownames(filtered.snp.ggm.overview))

combine.go.enrichments <- function(batches, batch.size) {
  load(paste0(GGM.DIR, 'ggm.bp.enrichment.1.RData'))
  total.global.enrichments <- global.enrichments
  for(batch in 2:batches) {
    load(paste0(GGM.DIR, 'ggm.bp.enrichment.', batch, '.RData'))
    total.global.enrichments <- append(total.global.enrichments, global.enrichments)
  }  
  
  global.enrichments <- total.global.enrichments
  save(global.enrichments, file = paste0(GGM.DIR, 'ggm.bp.enrichment.RData'))
  
}