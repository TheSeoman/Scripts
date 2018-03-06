source('Scripts/R/paths.R')
source('Scripts/R/lib.R')
source('Scripts/R/preprocess-ggm-data.methods.R')

load(PATHS$EXPR.RESIDUALS.DATA)
probe2gene <- unlist(as.list(illuminaHumanv3SYMBOL))
probe2gene <- probe2gene[!is.na(probe2gene)]

enrich.ggm.data.with.string.cpg.herv <- function(set = 'hervS2', filter = 'meth', seed = 'meqtl', flanking = 2.5e5, batch = NULL, batch.size = NULL) {
  GGM.NO.STRING.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb/')
  GGM.STRING.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb.string/')
  
  load(paste0(GGM.NO.STRING.DIR, 'cpgs.RData'))
  load(paste0(GGM.NO.STRING.DIR, 'data.meta.RData'))
  load(paste0(GGM.NO.STRING.DIR, 'data.overview.RData'))
  
  if(is.null(batch) | is.null(batch.size)) {
    range <- c(1:length(cpg.sets))
  } else {
    range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(cpg.sets)))    
  }
  
  load.string.db()
  
  for(seed.meth.ids in cpg.sets[range]) {
    set.name <- paste(seed.meth.ids, collapse = '|')
    
    meta <- data.meta[[set.name]]
    
    meth.ids <- c(seed.meth.ids, meta$extra.meth.ids)
    
  
    
    cat(paste0('Processing cpg-set: ', set.name), fill = T)
    load(paste0(paste0(GGM.NO.STRING.DIR, 'data/', set.name, '.RData')))
    snps <- meta$best.snps
    snp.genes <- meta$snp.genes
    
    tfbs.ann <- get.chipseq.context(meth.ids)
    cpgs.with.tfbs <- meth.ids[meth.ids %in% rownames(tfbs.ann[rowSums(tfbs.ann)>0, , drop = F])]
    
    snp.genes.in.string <- snp.genes[snp.genes %in% nodes(STRING.DB)]
    
    if(length(cpgs.with.tfbs) > 0 & length(snp.genes.in.string) > 0) {
      string.db <- STRING.DB
      string.db <- add.to.graphs(list(string.db), snps, snp.genes, cpgs.with.tfbs, tfbs.ann)[[1]]
      
      tfs = unique(unlist(adj(string.db, cpgs.with.tfbs)))
      
      nodeset = c(nodes(STRING.DB), setdiff(tfs, "KAP1"), snp.genes.in.string, cpgs.with.tfbs)
      string.db = subGraph(intersect(nodes(string.db), nodeset), string.db)
      
      path.genes <- get.string.shortest.paths(cis = cpgs.with.tfbs, 
                                              trans=unique(c(snp.genes.in.string, tfs)), 
                                              snp.genes=snp.genes.in.string,
                                              string.db)
      cat(paste(path.genes, collapse = ', '), fill = T)
    } else {
      cat('No cpgs.with.tfbs', fill = T)
      path.genes <- character(0)
    }
    
    
    extra.genes <- path.genes[!path.genes %in% c(snp.genes, meta$seed.meth.genes, meta$seed.tfbs.genes, meta$extra.meth.genes, meta$extra.tfbs.genes)]
    
    if(length(extra.genes) > 0 ){
      path.gene.data.list <- lapply(extra.genes, function(gene) {
        probe.ids <- unique(names(probe2gene[probe2gene == gene]))
        if (length(probe.ids) == 1) {
          expr <- expr.residuals[, probe.ids[1]]
        } else {
          expr <- apply(expr.residuals[, probe.ids], 1, function(x) mean(x))
        }
        return(expr)
      } )
      path.gene.data <- data.frame(matrix(unlist(path.gene.data.list), byrow=FALSE, ncol = length(path.gene.data.list)))
      colnames(path.gene.data) <- extra.genes
      rownames(path.gene.data) <- rownames(expr.residuals)
      ggm.data <- cbind.data.frame(ggm.data, path.gene.data[rownames(ggm.data), , drop = F])
    }
    data.meta[[set.name]]$path.genes <- path.genes
    data.overview[set.name, 'total.entities'] <- data.overview[set.name, 'total.entities'] + length(extra.genes) 
    data.overview[set.name, 'path.genes'] <- length(path.genes)
    
    save(ggm.data, file = paste0(GGM.STRING.DIR, 'data/', set.name, '.RData'))
  }
  save(data.overview, file = paste0(GGM.STRING.DIR, 'data.overview', batch, '.RData'))
  save(data.meta, file = paste0(GGM.STRING.DIR, 'data.meta', batch, '.RData'))
}


enrich.ggm.data.with.string <- function(set = 'hervS1', filter = 'snp', seed = 'meqtl', flanking = 5e5, batch = NULL, batch.size = NULL) {
  GGM.NO.STRING.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb/')
  GGM.STRING.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb.string/')
  
  load(paste0(GGM.NO.STRING.DIR, 'snps.RData'))
  load(paste0(GGM.NO.STRING.DIR, 'data.meta.RData'))
  load(paste0(GGM.NO.STRING.DIR, 'data.overview.RData'))
  
  if(is.null(batch) | is.null(batch.size)) {
    range <- c(1:length(snps))
  } else {
    range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(snps)))    
  }
  
  load.string.db()
  
  for(snp in snps[range]) {
    cat(paste0('Processing snp: ', snp), fill = T)
    load(paste0(paste0(GGM.NO.STRING.DIR, 'data/', snp, '.RData')))
    meth.ids <- data.meta[[snp]]$meth.ids
    snp.genes <- data.meta[[snp]]$snp.genes
    
    tfbs.ann <- get.chipseq.context(meth.ids)
    cpgs.with.tfbs <- meth.ids[meth.ids %in% rownames(tfbs.ann[rowSums(tfbs.ann)>0,])]
    
    if(length(cpgs.with.tfbs) > 0) {
      snp.genes.in.string <- snp.genes[snp.genes %in% nodes(STRING.DB)]
      
      string.db <- STRING.DB
      string.db <- add.to.graphs(list(string.db), snp, snp.genes, cpgs.with.tfbs, tfbs.ann)[[1]]
      
      tfs = unique(unlist(adj(string.db, cpgs.with.tfbs)))
      
      nodeset = c(nodes(STRING.DB), setdiff(tfs, "KAP1"), snp.genes.in.string, cpgs.with.tfbs)
      string.db = subGraph(intersect(nodes(string.db), nodeset), string.db)
      
      path.genes <- get.string.shortest.paths(cis = cpgs.with.tfbs, 
                                              trans=unique(c(snp.genes.in.string, tfs)), 
                                              snp.genes=snp.genes.in.string,
                                              string.db)
      cat(paste(path.genes, collapse = ', '), fill = T)
    } else {
      cat('No cpgs.with.tfbs', fill = T)
      path.genes <- character(0)
    }
    
    
    extra.genes <- path.genes[!path.genes %in% c(snp.genes, data.meta[[snp]]$meth.genes, data.meta[[snp]]$tfbs.genes)]
    
    if(length(extra.genes) > 0 ){
      path.gene.data.list <- lapply(extra.genes, function(gene) {
        probe.ids <- unique(names(probe2gene[probe2gene == gene]))
        if (length(probe.ids) == 1) {
          expr <- expr.residuals[, probe.ids[1]]
        } else {
          expr <- apply(expr.residuals[, probe.ids], 1, function(x) mean(x))
        }
        return(expr)
      } )
      path.gene.data <- data.frame(matrix(unlist(path.gene.data.list), byrow=FALSE, ncol = length(path.gene.data.list)))
      colnames(path.gene.data) <- extra.genes
      rownames(path.gene.data) <- rownames(expr.residuals)
      ggm.data <- cbind.data.frame(ggm.data, path.gene.data[rownames(ggm.data), ])
    }
    data.meta[[snp]]$path.genes <- path.genes
    data.overview[snp, 'total.entities'] <- data.overview[snp, 'total.entities'] + length(extra.genes) 
    data.overview[snp, 'path.genes'] <- length(path.genes)
    
    save(ggm.data, file = paste0(GGM.STRING.DIR, 'data/', snp, '.RData'))
  }
  save(data.overview, file = paste0(GGM.STRING.DIR, 'data.overview', batch, '.RData'))
  save(data.meta, file = paste0(GGM.STRING.DIR, 'data.meta', batch, '.RData'))
}