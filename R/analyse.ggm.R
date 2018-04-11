source('Scripts/R/paths.R')
source('Scripts/R/lib.R')
source('Scripts/R/go.enrichment.R')

require(illuminaHumanv3.db)
require(BDgraph)
require(GenomicRanges)

set <- 'hervS2'
filter <- 'snp'
seed <- 'meqtl'
flanking <- 2.5e5
string  <- T

load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$HERV.METH.OVERLAP.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)
load(PATHS$EQTM.ME.DATA)
load(PATHS$MAF001.RES.ME.DATA)

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
if (filter == 'snp' ) {
  load(paste0(GGM.DIR, 'filtered.snps.RData'))
} else {
  load(paste0(GGM.DIR, 'cpgs.RData'))
}
load(paste0(GGM.DIR, 'data.meta.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))


if (filter == 'snp' ) {
  cis.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('cis.', filter)]][, 1:2]
  trans.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('trans.', filter)]][, 1:2]
  all.eqtl.pairs <- rbind.data.frame(cis.eqtl.pairs, trans.eqtl.pairs, stringsAsFactors = F)
} else {
  all.eqtl.pairs <- rbind.data.frame(eqtl.me$cis$eqtls[,1:2], eqtl.me$trans$eqtls[,1:2], stringsAsFactors = F)
  colnames(all.eqtl.pairs) <- c('snp', 'expr.id')
}

herv.cpgs <-  names(get(paste0(set, '.meth.overlap'))$meth.ranges)

probe2gene <- unlist(as.list(illuminaHumanv3SYMBOL))
probe2gene <- probe2gene[!is.na(probe2gene)]

herv.expr.ids <- names(get(paste0(set, '.expr.overlap'))$expr.ranges)
herv.genes <- unique(probe2gene[herv.expr.ids[herv.expr.ids %in% names(probe2gene)]])

trans.eqtm.pairs <- eqtm.me$trans$eqtls[, c('snps', 'gene')]
cis.eqtm.pairs <- eqtm.me$cis$eqtls[, c('snps', 'gene')]
all.eqtm.pairs <- rbind.data.frame(trans.eqtm.pairs, cis.eqtm.pairs, stringsAsFactors = F)
colnames(all.eqtm.pairs) <- c('cpg', 'expr.id')
all.eqtm.pairs$cpg <- as.character(all.eqtm.pairs$cpg)
all.eqtm.pairs$expr.id <- as.character(all.eqtm.pairs$expr.id)

cutoff <- 0.9



analyse.cpg.ggms <- function(cpg.sets, cutoff = 0.9) {
  ggm.meta <- list()
  set.names <- sapply(cpg.sets, paste, collapse = '|')
  ggm.overview <- data.frame(matrix(nrow = length(set.names), ncol = 36))
  rownames(ggm.overview) <- set.names
  colnames(ggm.overview) <- c('entities', 'edges', 'ccs', 'cc.sizes', 'seed.cpg.cc.sizes',
                              'seed.cpgs', 'con.seed.cpgs', 'extra.cpgs', 'extra.con.cpgs', 'snps', 'con.snps', 
                              'tfs', 'con.tfs', 'path.genes', 'con.path.genes',
                              'seed.cpg.genes', 'con.seed.cpg.genes', 'direct.seed.cpg.genes', 'con.direct.seed.cpg.genes' ,'herv.seed.cpg.genes',
                              'extra.cpg.genes', 'con.extra.cpg.genes', 'direct.extra.cpg.genes', 'con.direct.extra.cpg.genes' ,'herv.extra.cpg.genes',
                              'snp.genes', 'con.snp.genes', 'direct.snp.genes', 'herv.snp.genes', 'con.herv.snp.genes', 
                              'pos.eqtl', 'con.eqtl', 'direct.eqtl', 
                              'pos.eqtm', 'con.eqtm', 'direct.eqtm')
  
  ggms <- list()
  for (set.name in set.names) {
    cat(paste0('Processing cpg set: ', set.name), fill = T)
    load(paste0(GGM.DIR, 'ggm/', set.name, '.RData'))
    ggms[[set.name]] <- ggm
    meta <- data.meta[[set.name]]
    ggm.res <- graph.from.fit(ggm, cutoff)
    ccs <- connComp(ggm.res)
    cc.lengths <- unlist(lapply(ccs, length))
    biggest.cc <- ccs[[which.max(cc.lengths)]]
    
    entities <- nodes(ggm.res)
    edges <- numEdges(ggm.res)
    seed.cpgs <- data.frame(row.names = meta$seed.meth.ids)
    seed.cpgs$con <- rownames(seed.cpgs) %in% biggest.cc

    seed.cpg.ccs <- ccs[sapply(ccs, function(cc) {
      return(length(intersect(rownames(seed.cpgs), cc)) > 0) 
    })]
    seed.cpg.cc.lengths <- unlist(lapply(seed.cpg.ccs, length))
    
    extra.cpgs <- data.frame(row.names = meta$extra.meth.ids)
    extra.cpgs$con <- rownames(extra.cpgs) %in% biggest.cc
    
    seed.cpg.genes <- data.frame(row.names = c(meta$seed.meth.genes, meta$seed.meth.no.gene.probes))
    seed.cpg.genes$con <- rownames(seed.cpg.genes) %in% biggest.cc
    seed.cpg.genes$direct <- rownames(seed.cpg.genes) %in% unlist(adj(ggm.res, rownames(seed.cpgs)))
    seed.cpg.genes$herv <- rownames(seed.cpg.genes) %in% herv.expr.ids | rownames(seed.cpg.genes) %in% herv.genes
    
    extra.cpg.genes <- data.frame(row.names = c(meta$extra.meth.genes, meta$extra.meth.no.gene.probes))
    if(nrow(extra.cpgs) > 0) {
      extra.cpg.genes$con <- rownames(extra.cpg.genes) %in% biggest.cc
      extra.cpg.genes$direct <- rownames(extra.cpg.genes) %in% unlist(adj(ggm.res, rownames(extra.cpgs)))
      extra.cpg.genes$herv <- rownames(extra.cpg.genes) %in% herv.expr.ids | rownames(extra.cpg.genes) %in% herv.genes
    }
    snps <- data.frame(row.names = meta$best.snps)
    snps$con <- rownames(snps) %in% biggest.cc
    
    snp.genes <- data.frame(row.names = c(meta$snp.genes, meta$snp.no.gene.probes))
    snp.genes$con <- rownames(snp.genes) %in% biggest.cc
    snp.genes$direct <- rownames(snp.genes) %in% unlist(adj(ggm.res, rownames(snps)))
    snp.genes$herv <- rownames(snp.genes) %in% herv.expr.ids | rownames(snp.genes) %in% herv.genes
    
    tfs <- data.frame(row.names = unique(c(meta$seed.tfbs.genes, meta$extra.tfbs.genes)))
    tfs$con <- rownames(tfs) %in% biggest.cc
    
    path.genes <- data.frame(row.names = meta$path.genes)
    path.genes$con <- rownames(path.genes) %in% biggest.cc
    
    all.expr.ids <- unique(c(names(probe2gene[probe2gene %in% c(meta$seed.meth.genes, meta$extra.meth.genes,
                                                                meta$seed.tfbs.genes, meta$extra.tfbs.genes,
                                                                meta$snp.genes, meta$path.genes)]), 
                             meta$snp.no.gene.probes, meta$seed.meth.no.gene.probes, meta$extra.meth.no.gene.probes))
    
    eqtl.pairs <- all.eqtl.pairs[all.eqtl.pairs$snp %in% snps & all.eqtl.pairs$expr.id %in% all.expr.ids,]
    eqtl.pairs$gene <- probe2gene[eqtl.pairs$expr.id]
    
    eqtl.pairs$con <- eqtl.pairs$snp %in% biggest.cc & (eqtl.pairs$expr.id %in% biggest.cc | eqtl.pairs$gene %in% biggest.cc)
    eqtl.pairs$direct <- apply(eqtl.pairs, 1, function(pair) {
      snp.edges <- adj(ggm.res, pair[1])
      return(pair[2] %in% snp.edges | pair[3] %in% snp.edges)
    })
    
    
    eqtm.pairs <- all.eqtm.pairs[all.eqtm.pairs$cpg %in% c(rownames(seed.cpgs), rownames(extra.cpgs)) & all.eqtm.pairs$expr.id %in% all.expr.ids, ]
    eqtm.pairs$gene <- probe2gene[eqtm.pairs$expr.id]
    
    eqtm.pairs$con <- eqtm.pairs$cpg %in% biggest.cc & (eqtm.pairs$expr.id %in% biggest.cc | eqtm.pairs$gene %in% biggest.cc)
    
    eqtm.pairs$direct <- apply(eqtm.pairs, 1, function(row) {
      cpg.edges <- adj(ggm.res, row[1])
      return(row[2] %in% cpg.edges | row[3] %in% cpg.edges)
    })
    
    ggm.overview[set.name,] <- c(length(entities), edges, length(ccs), paste(sort(cc.lengths[cc.lengths > 1], decreasing = T), collapse = ', '), seed.cpg.cc.lengths,
                            nrow(seed.cpgs), sum(seed.cpgs$con), nrow(extra.cpgs), sum(extra.cpgs$con), nrow(snps), sum(snps$con), 
                            nrow(tfs), sum(tfs$con), nrow(path.genes), sum(path.genes$con),
                            nrow(seed.cpg.genes), sum(seed.cpg.genes$con), sum(seed.cpg.genes$direct), sum(seed.cpg.genes$direct & seed.cpg.genes$con), sum(seed.cpg.genes$herv), 
                            nrow(extra.cpg.genes), sum(extra.cpg.genes$con), sum(extra.cpg.genes$direct), sum(extra.cpg.genes$direct & extra.cpg.genes$con), sum(extra.cpg.genes$herv), 
                            nrow(snp.genes), sum(snp.genes$con), sum(snp.genes$direct), sum(snp.genes$herv), sum(snp.genes$herv & snp.genes$con),
                            nrow(eqtl.pairs), sum(eqtl.pairs$con), sum(eqtl.pairs$direct), 
                            dim(eqtm.pairs)[1], sum(eqtm.pairs$con), sum(eqtm.pairs$direct))
    
    ggm.meta[[set.name]] <- list(ccs = ccs, seed.cpgs = seed.cpgs, extra.cpgs = extra.cpgs, 
                            seed.cpg.genes = seed.cpg.genes, extra.cpg.genes = extra.cpg.genes, snps = snps,
                            snp.genes = snp.genes, tfs = tfs, path.genes = path.genes, eqtl.pairs = eqtl.pairs, eqtm.pairs = eqtm.pairs)
    
  }
  for(i in c(1:3, 5:36)) {
    ggm.overview[, i] <- as.integer(ggm.overview[, i])
  }
  
  save(ggm.overview, file = paste0(GGM.DIR, 'ggm.overview.', cutoff,'.RData'))
  save(ggm.meta, file = paste0(GGM.DIR, 'ggm.meta', cutoff, '.RData'))
}

###snp

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.snp.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')

analyse.snp.ggms <- function(snps, cutoff = 0.9) {
  ggm.meta <- list()
  ggm.overview <- data.frame(matrix(nrow = length(snps), ncol = 31))
  rownames(ggm.overview) <- snps
  colnames(ggm.overview) <- c('entities', 'edges', 'ccs', 'cc.sizes', 'con.snp', 'snp.cc.size',
                              'cpgs', 'con.cpgs', 'direct.cpgs', 
                              'herv.cpgs', 'con.herv.cpgs', 'cpg.genes', 'con.cpg.genes', 'direct.cpg.genes', 'con.direct.cpg.genes' ,'herv.cpg.genes',
                              'snp.genes', 'con.snp.genes', 'direct.snp.genes', 'herv.snp.genes', 'con.herv.snp.genes', 
                              'tfs', 'con.tfs', 'path.genes', 'con.path.genes',
                              'pos.eqtl', 'con.eqtl', 'direct.eqtl', 
                              'pos.eqtm', 'con.eqtm', 'direct.eqtm')
  
  ggms <- list()
  for (snp in snps) {
    cat(paste0('Processing snp: ', snp), fill = T)
    load(paste0(GGM.DIR, 'ggm/', snp, '.RData'))
    meta <- data.meta[[snp]]
    ggm.res <- graph.from.fit(ggm, cutoff)
    ccs <- connComp(ggm.res)
    cc.lengths <- unlist(lapply(ccs, length))
    biggest.cc <- ccs[[which.max(cc.lengths)]]
    
    con.snp <- sum(snp %in% biggest.cc)
    snp.cc <- ccs[[1]]
        
    entities <- nodes(ggm.res)
    edges <- numEdges(ggm.res)
    cpgs <- data.frame(row.names = meta$meth.ids)
    cpgs$con <- rownames(cpgs) %in% biggest.cc
    cpgs$direct <- rownames(cpgs) %in% adj(ggm.res, snp)[[1]]
    cpgs$herv <- rownames(cpgs) %in% herv.cpgs
    
    
    cpg.genes <- data.frame(row.names = c(meta$meth.genes, meta$meth.no.gene.probes))
    cpg.genes$con <- rownames(cpg.genes) %in% biggest.cc
    cpg.genes$direct <- rownames(cpg.genes) %in% unlist(adj(ggm.res, rownames(cpgs)))
    cpg.genes$herv <- rownames(cpg.genes) %in% herv.expr.ids | rownames(cpg.genes) %in% herv.genes
    
    snp.genes <- data.frame(row.names = c(meta$snp.genes, meta$snp.no.gene.probes))
    snp.genes$con <- rownames(snp.genes) %in% biggest.cc
    snp.genes$direct <- rownames(snp.genes) %in% adj(ggm.res, snp)[[1]]
    snp.genes$herv <- rownames(snp.genes) %in% herv.expr.ids | rownames(snp.genes) %in% herv.genes
    
    tfs <- data.frame(row.names = meta$tfbs.genes)
    tfs$con <- rownames(tfs) %in% biggest.cc
    
    path.genes <- data.frame(row.names = meta$path.genes)
    path.genes$con <- rownames(path.genes) %in% biggest.cc
    
  
    eqtl.expr.ids <- unique(as.character(all.eqtl.pairs[all.eqtl.pairs$snps == snp, 'gene']))
    eqtl.genes <- unique(probe2gene[eqtl.expr.ids[eqtl.expr.ids %in% names(probe2gene)]])
    
    eqtls <- data.frame(row.names = entities[entities %in% eqtl.expr.ids | entities %in% eqtl.genes])
    eqtls$con <- rownames(eqtls) %in% biggest.cc
    eqtls$direct <- rownames(eqtls) %in% adj(ggm.res, snp)[[1]]
    
    all.expr.ids <- unique(c(names(probe2gene[probe2gene %in% c(meta$tfbs.genes, meta$snp.genes, meta$meth.genes)]), meta$snp.no.gene.probes, meta$meth.no.gene.probes))
    
    eqtm.pairs <- all.eqtm.pairs[all.eqtm.pairs$cpg %in% rownames(cpgs) & all.eqtm.pairs$expr.id %in% all.expr.ids, ]
    eqtm.pairs$gene <- probe2gene[eqtm.pairs$expr.id]
    
    eqtm.pairs$con <- eqtm.pairs$cpg %in% biggest.cc & (eqtm.pairs$expr.id %in% biggest.cc | eqtm.pairs$gene %in% biggest.cc)
    
    eqtm.pairs$direct <- apply(eqtm.pairs, 1, function(row) {
      cpg.edges <- adj(ggm.res, row[1])
      return(row[2] %in% cpg.edges | row[3] %in% cpg.edges)
    })
    ggm.overview[snp,] <- c(length(entities), edges, length(ccs), paste(sort(cc.lengths[cc.lengths > 1], decreasing = T), collapse = ', '), con.snp, length(snp.cc),
                            dim(cpgs)[1], sum(cpgs$con), sum(cpgs$direct), sum(cpgs$herv), sum(cpgs$herv & cpgs$con), 
                            dim(cpg.genes)[1], sum(cpg.genes$con), sum(cpg.genes$direct), sum(cpg.genes$direct & cpg.genes$con), sum(cpg.genes$herv), 
                            dim(snp.genes)[1], sum(snp.genes$con), sum(snp.genes$direct), sum(snp.genes$herv), sum(snp.genes$herv & snp.genes$con),
                            nrow(tfs), sum(tfs$con), nrow(path.genes), sum(path.genes$con),
                            dim(eqtls)[1], sum(eqtls$con), sum(eqtls$direct), 
                            dim(eqtm.pairs)[1], sum(eqtm.pairs$con), sum(eqtm.pairs$direct))
    
    ggm.meta[[snp]] <- list(ccs = ccs, cpgs = cpgs, cpg.genes = cpg.genes, snp.genes = snp.genes, tfs = tfs, path.genes = path.genes, eqtls = eqtls, eqtm.pairs = eqtm.pairs)
    
  }
  for(i in c(1:3, 5:31)) {
    ggm.overview[, i] <- as.integer(ggm.overview[, i])
  }
  
  save(ggm.overview, file = paste0(GGM.DIR, 'ggm.overview.', cutoff,'.RData'))
  save(ggm.meta, file = paste0(GGM.DIR, 'ggm.meta', cutoff, '.RData'))
}

analyse.snp.ggms(filtered.snps)
