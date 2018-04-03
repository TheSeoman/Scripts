source('Scripts/R/paths.R')
source('Scripts/R/preprocess-ggm-data.methods.R')

require('GenomicRanges')
require('illuminaHumanv3.db')

require(Rsamtools)

load(PATHS$HERV.MEQTL.TRANS.OVERLAP.DATA)
load(PATHS$MEQTL.TRANS.PAIRS.DATA)
load(PATHS$METH.RESIDUALS.DATA)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$METH.TFBS.OVERLAP.DATA)
load(PATHS$HERV.METH.OVERLAP.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)

load(PATHS$EXPR.GENE.ANNOT.DATA)


load(PATHS$EXPR.RANGES.DATA)
load(PATHS$METH.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

load(PATHS$SNP.SAMPLES.DATA)

load(PATHS$ID.MAP.DATA)

enlarge.ranges <- function(ranges, flanking) {
  enlarged.ranges <- GRanges(
    seqnames = seqnames(ranges),
    ranges = IRanges(start = start(ranges) - flanking, end = end(ranges) + flanking),
    strand = strand(ranges),
    name = ranges$name,
    score = ranges$score
  )
  names(enlarged.ranges) <- names(ranges)
  return (enlarged.ranges)
}

get.snp.data <- function(snp.range, snp.samples) {
  data = scanTabix(PATHS$F.SNP, param=snp.range)
  snp.data.list <- lapply(data, function(x) strsplit(x, '\t'))
  snp.data.table <- data.frame(matrix(unlist(snp.data.list), nrow=length(snp.samples)+5, byrow=F), stringsAsFactors = FALSE)
  colnames(snp.data.table) <- snp.data.table[2, ]
  snp.data.table <- snp.data.table[-(1:5), names(snp.range), drop = FALSE]
  snp.data.table <- data.frame(data.matrix(snp.data.table))
  rownames(snp.data.table) <- snp.samples
  return(snp.data.table)
}

get.nearby.probes <- function(snp.range, expr.ranges, distance = 5e5, overlap.type = 'any') {
  area.range <- enlarge.ranges(snp.range, distance)
  overlap.hits <- findOverlaps(area.range, expr.ranges, type = overlap.type)
  expr.ids <- names(expr.ranges[unique(subjectHits(overlap.hits))])
  return(expr.ids)
}

get.neighbour.probes <- function(meth.ranges, expr.ranges, max.distance = 5e5) {
  overlap.hits <- findOverlaps(meth.ranges, expr.ranges, type = 'any')
  precede.indices <- precede(meth.ranges, expr.ranges, ignore.strand = T)
  follow.indices <- follow(meth.ranges, expr.ranges, ignore.strand = T)
  
  precede.ranges <- expr.ranges[precede.indices]
  follow.ranges <- expr.ranges[follow.indices]
  
  expr.ids <- unique(c(names(expr.ranges[subjectHits(overlap.hits)]),
                       names(precede.ranges[distance(meth.ranges, precede.ranges) < max.distance]),
                       names(follow.ranges[distance(meth.ranges, follow.ranges) < max.distance])))
  
  return(expr.ids)    
}

collect.pair.overlaps <- function(pairs) {
  gn <- graphNEL(nodes=unique(c(pairs[,1], pairs[,2])))
  gn <- addEdge(pairs[,1], pairs[,2], gn)
  ccs <- connectedComp(gn)
  pair.overlaps <- lapply(ccs, function(cc) {
    cpg.set <- cc[grep('cg', cc)]
    herv.set <- cc[grep('chr', cc)]
    return(cpg.set)
  })  
  return(pair.overlaps)
}

if(F){
  set <- 'hervS2'
  filter <- 'snp'
  seed <- 'meqtl'
  snp.count.threshold = 5
  flanking <- 2.5e5
  string <- T
}



prepare.ggm.data.cpg.herv <- function(set = 'hervS2', filter = 'meth', seed = 'meqtl', string = T, flanking = 2.5e5, batch = NULL, batch.size = NULL) {
  GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
  dir.create(paste0(GGM.DIR, 'data/'), showWarnings = F, recursive = T)
  data <- list()
  data.meta <- list()
  herv.meqtl.pairs <- get(paste0(set, '.meqtl.trans.overlap'))[[filter]]
  herv.meqtl.pairs <- herv.meqtl.pairs[herv.meqtl.pairs$snp %in% names(snp.ranges), ]
  herv.expr.pairs <- get(paste0(set, '.expr.overlap'))$pairs
  
  if(!file.exists(paste0(GGM.DIR, 'cpgs.RData'))) {
    cpgs <- as.character(unique(herv.meqtl.pairs$cpg))
    herv.meth.pairs <- get(paste0(set, '.meth.overlap'))$pairs
    herv.meth.pairs <- herv.meth.pairs[herv.meth.pairs$meth.id %in% cpgs, ]
    cpg.sets <- collect.pair.overlaps(herv.meth.pairs)
    save(cpg.sets, file = paste0(GGM.DIR, 'cpgs.RData'))
  } else {
    load(paste0(GGM.DIR, 'cpgs.RData'))
  }
  
  data.overview <- data.frame(matrix(ncol = ifelse(string, 13, 12), nrow = length(cpg.sets)))
  if (string) {
    colnames(data.overview) <- c('total.entities', 'seed.cpgs', 'extra.cpgs', 'best.snps',  'seed.meth.tfs', 'extra.meth.tfs', 'snp.genes', 'snp.no.gene.probes', 'seed.meth.genes', 'seed.meth.no.gene.probes', 'extra.meth.genes', 'extra.meth.no.gene.probes', 'path.genes')
  } else {
    colnames(data.overview) <- c('total.entities', 'seed.cpgs', 'extra.cpgs', 'best.snps',  'seed.meth.tfs', 'extra.meth.tfs', 'snp.genes', 'snp.no.gene.probes', 'seed.meth.genes', 'seed.meth.no.gene.probes', 'extra.meth.genes', 'extra.meth.no.gene.probes')
  }
  rownames(data.overview) <- unlist(lapply(cpg.sets, paste, collapse = '|')) 
  
  
  if(string) {
    load.string.db()
  }
  
  if(is.null(batch) | is.null(batch.size)) {
    range <- c(1:length(cpg.sets))
  } else {
    range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(cpgs.sets)))    
  }
  
  for (seed.meth.ids in cpg.sets[range]) {
    set.name <- paste(seed.meth.ids, collapse = '|')
    cat(paste0('Processing cpg-set: ', set.name), fill = T)
    best.snps <- unique(sapply(seed.meth.ids, function(meth.id) {
      pairs <- herv.meqtl.pairs[herv.meqtl.pairs$cpg == meth.id, ]
      return(as.character(pairs[which.min(pairs$p.comb)[1], 'snp']))  
    }))
    
    extra.meth.ids <- unique(as.character(meqtl.trans.pairs[meqtl.trans.pairs$snp %in% best.snps, 'cpg']))
    extra.meth.ids <- extra.meth.ids[!extra.meth.ids %in% seed.meth.ids]
    
    best.snp.ranges <- snp.ranges[best.snps]
    
    snp.data <- get.snp.data(best.snp.ranges, snp.samples)
    
    snp.expr.ids <- get.nearby.probes(best.snp.ranges, expr.ranges, flanking)
    snp.expr.no.gene.ids <- snp.expr.ids[!snp.expr.ids %in% names(probe2gene)]
    snp.expr.with.gene.ids <- snp.expr.ids[snp.expr.ids %in% names(probe2gene)]
    snp.genes <- unique(probe2gene[snp.expr.with.gene.ids])
    
    seed.meth.expr.ids <- get.neighbour.probes(meth.ranges[seed.meth.ids], expr.ranges, flanking)
    seed.meth.expr.no.gene.ids <- seed.meth.expr.ids[!seed.meth.expr.ids %in% names(probe2gene)]
    seed.meth.expr.with.gene.ids <- seed.meth.expr.ids[seed.meth.expr.ids %in% names(probe2gene)]
    seed.meth.genes <- unique(probe2gene[seed.meth.expr.with.gene.ids])
    
    if (length(extra.meth.ids) > 0) {
      extra.meth.expr.ids <- get.neighbour.probes(meth.ranges[extra.meth.ids], expr.ranges, flanking)
      extra.meth.expr.no.gene.ids <- extra.meth.expr.ids[!extra.meth.expr.ids %in% names(probe2gene)]
      extra.meth.expr.with.gene.ids <- extra.meth.expr.ids[extra.meth.expr.ids %in% names(probe2gene)]
      extra.meth.genes <- unique(probe2gene[extra.meth.expr.with.gene.ids])
    } else {
      extra.meth.expr.no.gene.ids <- character(0)
      extra.meth.genes <- character(0)
    }

    expr.no.gene.data <- expr.residuals[, unique(c(snp.expr.no.gene.ids, seed.meth.expr.no.gene.ids, extra.meth.expr.no.gene.ids)), drop = F]
    
    meth.data <- meth.residuals[, unique(c(seed.meth.ids, extra.meth.ids)), drop = F]
    
    seed.tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% seed.meth.ids, 'tfbs.id'])
    seed.tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[seed.tfbs.ids]$TF)
    extra.tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% extra.meth.ids, 'tfbs.id'])
    extra.tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[extra.tfbs.ids]$TF)
    
    data.meta[[set.name]] <- list(seed.meth.ids = seed.meth.ids, 
                                  extra.meth.ids = extra.meth.ids, 
                                  best.snps = best.snps, 
                                  seed.tfbs.genes = seed.tfbs.genes, 
                                  extra.tfbs.genes = extra.tfbs.genes, 
                                  snp.genes = snp.genes, 
                                  snp.no.gene.probes = snp.expr.no.gene.ids, 
                                  seed.meth.genes = seed.meth.genes, 
                                  seed.meth.no.gene.probes = seed.meth.expr.no.gene.ids,
                                  extra.meth.genes = extra.meth.genes, 
                                  extra.meth.no.gene.probes = extra.meth.expr.no.gene.ids)
    
    total.genes <- unique(c(snp.genes, seed.meth.genes, extra.meth.genes, seed.tfbs.genes, extra.tfbs.genes))
    
    overview <-  c(length(seed.meth.ids), length(extra.meth.ids), length(best.snps), length(seed.tfbs.genes), length(extra.tfbs.genes), length(snp.genes), length(snp.expr.no.gene.ids), 
                   length(seed.meth.genes), length(seed.meth.expr.no.gene.ids), length(extra.meth.genes), length(extra.meth.expr.no.gene.ids))
    
    if(string) {
      meth.ids <- unique(c(seed.meth.ids, extra.meth.ids))
      tfbs.ann <- get.chipseq.context(meth.ids)
      cpgs.with.tfbs <- meth.ids[meth.ids %in% rownames(tfbs.ann[rowSums(tfbs.ann)>0,])]
      if (length(cpgs.with.tfbs) > 0) {
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
      } else {
        path.genes <- character(0)
      }
      
      total.genes <- unique(c(total.genes, path.genes))
      
      data.meta[[set.name]][['path.genes']] <- path.genes
      overview <- c(overview, length(path.genes))
    }
    
    if(length(total.genes) > 0 ){
      expr.gene.data.list <- lapply(total.genes, function(gene) {
        probe.ids <- unique(names(probe2gene[probe2gene == gene]))
        if (length(probe.ids) == 1) {
          expr <- expr.residuals[, probe.ids[1]]
        } else {
          expr <- apply(expr.residuals[, probe.ids], 1, function(x) mean(x))
        }
        return(expr)
      } )
      expr.gene.data <- data.frame(matrix(unlist(expr.gene.data.list), byrow=FALSE, ncol = length(expr.gene.data.list)))
      colnames(expr.gene.data) <- total.genes
      rownames(expr.gene.data) <- rownames(expr.residuals)
    } else {
      expr.gene.data <- data.frame(matrix(ncol = 0, nrow = nrow(expr.residuals)))
    }
    
    ggm.data <- cbind.data.frame(snp.data[id.map$axio_s4f4, , drop=F], meth.data[id.map$meth_f4, , drop=F], expr.no.gene.data[id.map$expr_s4f4ogtt, , drop=F], expr.gene.data[id.map$expr_s4f4ogtt,])
    rownames(ggm.data) <- id.map$expr_s4f4ogtt
    
    overview <- c(dim(ggm.data)[2], overview)
    
    data.overview[set.name,] <- overview
    
    save(ggm.data, file = paste0(GGM.DIR, 'data/', set.name, '.RData'))
  }
  save(data.overview, file = paste0(GGM.DIR, 'data.overview', batch, '.RData'))
  save(data.meta, file = paste0(GGM.DIR, 'data.meta', batch, '.RData'))
}

prepare.ggm.data.snp.herv <- function(set = 'hervS2', filter = 'snp', seed = 'meqtl', string = T, snp.count.threshold = 5, flanking = 2.5e5, batch = NULL, batch.size = NULL) {
  GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
  dir.create(paste0(GGM.DIR, 'data/'), showWarnings = F, recursive = T)
  data <- list()
  data.meta <- list()
  meqtl.pairs <- get(paste0(set, '.meqtl.trans.overlap'))[[filter]]
  
  if(!file.exists(paste0(GGM.DIR, 'snps.RData'))) {
    meqtl.count <- table(meqtl.pairs$snp)[table(meqtl.pairs$snp) > 0]
    snps <- names(meqtl.count[meqtl.count >= snp.count.threshold])
    save(snps, file = paste0(GGM.DIR, 'snps.RData'))
  } else {
    load(paste0(GGM.DIR, 'snps.RData'))
  }
  
  data.overview <- data.frame(matrix(ncol = ifelse(string, 8, 7), nrow = length(snps)))
  if (string) {
    colnames(data.overview) <- c('cpgs', 'TFs', 'snp.genes', 'snp.no.gene.probes', 'meth.genes', 'meth.no.gene.probes', 'path.genes', 'total.entities')
  } else {
    colnames(data.overview) <- c('cpgs', 'TFs', 'snp.genes', 'snp.no.gene.probes', 'meth.genes', 'meth.no.gene.probes', 'total.entities')
  }
  rownames(data.overview) <- snps
  
  if(string) {
    load.string.db()
  }
  
  if(is.null(batch) | is.null(batch.size)) {
    range <- c(1:length(snps))
  } else {
    range <- c(((batch-1)*batch.size+1):min(batch*batch.size, length(snps)))   
  }
  
  for (snp in snps[range]) {
    cat(paste0('Processing snp: ', snp), fill = T)
    snp.range <- snp.ranges[snp]
    snp.expr.ids <- get.nearby.probes(snp.range, expr.ranges, flanking)
    snp.expr.no.gene.ids <- snp.expr.ids[!snp.expr.ids %in% names(probe2gene)]
    snp.expr.with.gene.ids <- snp.expr.ids[snp.expr.ids %in% names(probe2gene)]
    snp.genes <- unique(probe2gene[snp.expr.with.gene.ids])
    
    meth.ids <- as.character(meqtl.pairs[meqtl.pairs$snp == snp, 'cpg'])
    meth.ids <- meth.ids[meth.ids %in% colnames(meth.residuals) & meth.ids %in% names(meth.ranges)]
    meth.expr.ids <- get.neighbour.probes(meth.ranges[meth.ids], expr.ranges, flanking)
    meth.expr.no.gene.ids <- meth.expr.ids[!meth.expr.ids %in% names(probe2gene)]
    meth.expr.with.gene.ids <- meth.expr.ids[meth.expr.ids %in% names(probe2gene)]
    meth.genes <- unique(probe2gene[meth.expr.with.gene.ids])
    expr.no.gene.data <- expr.residuals[, unique(c(snp.expr.no.gene.ids, meth.expr.no.gene.ids)), drop = F]
    
    meth.data <- meth.residuals[, meth.ids]
    tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% meth.ids, 'tfbs.id'])
    tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[tfbs.ids]$TF)
    
    data.meta[[snp]] <- list(meth.ids = meth.ids, tfbs.genes = tfbs.genes, snp.genes = snp.genes, snp.no.gene.probes = snp.expr.no.gene.ids, 
                             meth.genes = meth.genes, meth.no.gene.probes = meth.expr.no.gene.ids)
    
    total.genes <- unique(c(snp.genes, meth.genes, tfbs.genes))
    
    overview <-  c(length(meth.ids), length(tfbs.genes), length(snp.genes), length(snp.expr.no.gene.ids), 
                   length(meth.genes), length(meth.expr.no.gene.ids))
    
    if(string) {
      tfbs.ann <- get.chipseq.context(meth.ids)
      cpgs.with.tfbs <- meth.ids[meth.ids %in% rownames(tfbs.ann[rowSums(tfbs.ann)>0,])]
      if (length(cpgs.with.tfbs) > 0) {
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
      } else {
        path.genes <- character(0)
      }
      
      total.genes <- unique(c(total.genes, path.genes))
      
      data.meta[[snp]][['path.genes']] <- path.genes
      overview <- c(overview, length(path.genes))
    }
    
    expr.gene.data.list <- lapply(total.genes, function(gene) {
      probe.ids <- unique(names(probe2gene[probe2gene == gene]))
      if (length(probe.ids) == 1) {
        expr <- expr.residuals[, probe.ids[1]]
      } else {
        expr <- apply(expr.residuals[, probe.ids], 1, function(x) mean(x))
      }
      return(expr)
    } )
    expr.gene.data <- data.frame(matrix(unlist(expr.gene.data.list), byrow=FALSE, ncol = length(expr.gene.data.list)))
    colnames(expr.gene.data) <- total.genes
    rownames(expr.gene.data) <- rownames(expr.residuals)
    
    snp.data <- get.snp.data(snp.range, snp.samples)
    
    ggm.data <- cbind.data.frame(snp.data[id.map$axio_s4f4, , drop=F], meth.data[id.map$meth_f4, , drop=F], expr.no.gene.data[id.map$expr_s4f4ogtt, , drop=F], expr.gene.data[id.map$expr_s4f4ogtt,])
    rownames(ggm.data) <- id.map$expr_s4f4ogtt
    
    overview <- c(overview, dim(ggm.data)[2])
    
    data.overview[snp,] <- overview
    
    save(ggm.data, file = paste0(GGM.DIR, 'data/', snp, '.RData'))
  }
  save(data.overview, file = paste0(GGM.DIR, 'data.overview', batch, '.RData'))
  save(data.meta, file = paste0(GGM.DIR, 'data.meta', batch, '.RData'))
}