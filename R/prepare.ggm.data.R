source('Scripts/R/paths.R')
source('Scripts/R/residuals.R')
source('Scripts/R/preprocess-ggm-data.methods.R')

require('GenomicRanges')
require('illuminaHumanv3.db')

require(Rsamtools)

load(PATHS$HERV.MEQTL.TRANS.OVERLAP.DATA)
load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$METH.COV.MATRIX.DATA)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$METH.TFBS.OVERLAP.DATA)


load(PATHS$EXPR.RANGES.DATA)
load(PATHS$METH.RANGES.DATA)
load(PATHS$SNP.RANGES.DATA)

snp.samples <- scan(PATHS$F.SNP.SAMPLES)

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
  snp.data.table[, names(snp.range)] <- as.numeric(snp.data.table[, names(snp.range)])
  rownames(snp.data.table) <- snp.samples
  return(snp.data.table)
}

get.nearby.probes <- function(snp.range, expr.ranges, distance = 5e5, overlap.type = 'any') {
  area.range <- enlarge.ranges(snp.range, distance)
  overlap.hits <- findOverlaps(area.range, expr.ranges, type = overlap.type)
  expr.ids <- names(expr.ranges[subjectHits(overlap.hits)])
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

covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)
id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                         & covariates.all$axio_s4f4 %in% snp.samples 
                         & covariates.all$meth_f4 %in% rownames(meth.matrix), c('expr_s4f4ogtt', 'axio_s4f4', 'meth_f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]
id.map$expr_s4f4ogtt <- as.character(id.map$expr_s4f4ogtt)
id.map$axio_s4f4 <- as.character(id.map$axio_s4f4)
id.map$meth_f4 <- as.character(id.map$meth_f4)

probe2gene <- unlist(as.list(illuminaHumanv3SYMBOL))
probe2gene <- probe2gene[!is.na(probe2gene)]

if(F){
  set <- 'hervS1'
  filter <- 'snp'
  seed <- 'meqtl'
  snp.count.threshold = 5
  flanking <- 5e5
  string <- F
}

prepare.ggm.data <- function(set = 'hervS1', filter = 'snp', seed = 'meqtl', string = F, snp.count.threshold = 5, flanking = 5e5) {
  GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
  dir.create(paste0(GGM.DIR, 'data/'), showWarnings = F, recursive = T)
  
  meqtl.pairs <- get(paste0(set, '.meqtl.trans.overlap'))[[filter]]
  meqtl.count <- table(meqtl.pairs$snp)[table(meqtl.pairs$snp) > 0]
  snps <- names(meqtl.count[meqtl.count >= snp.count.threshold])
  
  data <- list()
  
  data.meta <- list()
  
  data.overview <- data.frame(matrix(ncol = ifelse(string, 8, 7), nrow = length(snps)))
  if (string) {
    colnames(data.overview) <- c('cpgs', 'TFs', 'snp.genes', 'snp.no.gene.probes', 'meth.genes', 'meth.no.gene.probes', 'path.genes', 'total.entities')
  } else {
    colnames(data.overview) <- c('cpgs', 'TFs', 'snp.genes', 'snp.no.gene.probes', 'meth.genes', 'meth.no.gene.probes', 'total.entities')
  }
  rownames(data.overview) <- snps
  
  save(snps, file = paste0(GGM.DIR, 'snps.RData'))
  
  if(string) {
    load.string.db()
  }
  
  for (snp in snps) {
    cat(paste0('Processing snp: ', snp), fill = T)
    snp.range <- snp.ranges[snp]
    snp.expr.ids <- get.nearby.probes(snp.range, expr.ranges, flanking)
    snp.expr.no.gene.ids <- snp.expr.ids[!snp.expr.ids %in% names(probe2gene)]
    snp.expr.with.gene.ids <- snp.expr.ids[snp.expr.ids %in% names(probe2gene)]
    snp.genes <- unique(probe2gene[snp.expr.with.gene.ids])
    
    meth.ids <- as.character(meqtl.pairs[meqtl.pairs$snp == snp, 'cpg'])
    meth.expr.ids <- get.neighbour.probes(meth.ranges[meth.ids], expr.ranges, flanking)
    meth.expr.no.gene.ids <- meth.expr.ids[!meth.expr.ids %in% names(probe2gene)]
    meth.expr.with.gene.ids <- meth.expr.ids[meth.expr.ids %in% names(probe2gene)]
    meth.genes <- unique(probe2gene[meth.expr.with.gene.ids])
    
    expr.no.gene.data <- expr.residuals[, unique(c(snp.expr.no.gene.ids, meth.expr.no.gene.ids)), drop = F]
    
    meth.data <- get.residuals(meth.matrix[, c(1:25, which(colnames(meth.matrix) %in% meth.ids))], 'meth', meth.ids)
    tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% meth.ids, 'tfbs.id'])
    tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[tfbs.ids]$TF)
    
    data.meta[[snp]] <- list(meth.ids = meth.ids, tfbs.genes = tfbs.genes, snp.genes = snp.genes, snp.no.gene.probes = snp.expr.no.gene.ids, 
                             meth.genes = meth.genes, meth.no.gene.probes = meth.expr.no.gene.ids)
    
    total.genes <- unique(c(snp.genes, meth.genes, tfbs.genes))
    
    overview <-  c(length(meth.ids), length(tfbs.genes), length(snp.genes), length(snp.expr.no.gene.ids), 
                   length(meth.genes), length(meth.expr.no.gene.ids))
    
      tfbs.ann <- get.chipseq.context(meth.ids)
      cpgs.with.tfbs <- meth.ids[meth.ids %in% rownames(tfbs.ann[rowSums(tfbs.ann)>0,])]
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
    
    snp.data <- get.snp.data(snp.range)
    
    ggm.data <- cbind.data.frame(snp.data[id.map$axio_s4f4, , drop=F], meth.data[id.map$meth_f4, , drop=F], expr.no.gene.data[id.map$expr_s4f4ogtt, , drop=F], expr.gene.data[id.map$expr_s4f4ogtt,])
    rownames(ggm.data) <- id.map$expr_s4f4ogtt
    
    overview <- c(overview, dim(ggm.data)[2])
    
    data.overview[snp,] <- overview
    
    save(ggm.data, file = paste0(GGM.DIR, 'data/', snp, '.RData'))
  }
  save(data.overview, file = paste0(GGM.DIR, 'data.overview.RData'))
  save(data.meta, file = paste0(GGM.DIR, 'data.meta.RData'))
}

prepare.ggm.data(set = 'hervS1', filter = 'snp', seed = 'meqtl', string = T, snp.count.threshold = 5, flanking = 2.5e5)
