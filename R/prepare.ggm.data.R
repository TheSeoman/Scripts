source('Scripts/R/paths.R')
source('Scripts/R/residuals.R')

require('GenomicRanges')
require('illuminaHumanv3.db')

require(Rsamtools)

require(BDgraph)
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

get.snp.data <- function(snp.range) {
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

genes <- unlist(as.list(illuminaHumanv3SYMBOL))
genes <- genes[!is.na(genes)]

set <- 'hervS1'
filter <- 'meth'
seed <- 'meqtl'

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '/')

snp.count.threshold <- 5
meqtl.pairs <- get(paste0(set, '.meqtl.trans.overlap'))[[filter]]
meqtl.count <- table(meqtl.pairs$cpg)[table(meqtl.pairs$cpg) > 0]
snps <- names(meqtl.count[meqtl.count >= snp.count.threshold])

data <- list()

data.overview <- data.frame(matrix(ncol = 7, nrow = length(snps)))
colnames(data.overview) <- c('cpgs', 'TFs', 'snp.expr.genes', 'snp.no.gene.probes', 
                             'meth.expr.genes', 'meth.no.gene.probes', 'total.entities')
rownames(data.overview) <- snps

save(snps, file = paste0(GGM.DIR, 'snps.RData'))

for (snp in snps) {
  cat(paste0('Processing snp: ', snp), fill = T)
  snp.range <- snp.ranges[snp]
  snp.expr.ids <- get.nearby.probes(snp.range, expr.ranges)
  snp.expr.no.gene.ids <- snp.expr.ids[!snp.expr.ids %in% names(genes)]
  snp.expr.with.gene.ids <- snp.expr.ids[snp.expr.ids %in% names(genes)]
  snp.expr.genes <- unique(genes[snp.expr.with.gene.ids])
  
  meth.ids <- as.character(meqtl.pairs[meqtl.pairs$snp == snp, 'cpg'])
  meth.expr.ids <- get.neighbour.probes(meth.ranges[meth.ids], expr.ranges)
  meth.expr.no.gene.ids <- meth.expr.ids[!meth.expr.ids %in% names(genes)]
  meth.expr.with.gene.ids <- meth.expr.ids[meth.expr.ids %in% names(genes)]
  meth.expr.genes <- unique(genes[meth.expr.with.gene.ids])
  
  expr.no.gene.data <- expr.residuals[, unique(c(snp.expr.no.gene.ids, meth.expr.no.gene.ids)), drop = F]
  
  meth.data <- get.residuals(meth.matrix[, c(1:25, which(colnames(meth.matrix) %in% meth.ids))], 'meth', meth.ids)
  meth.tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% meth.ids, 'tfbs.id'])
  meth.tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[meth.tfbs.ids]$TF)
  
  total.genes <- unique(c(snp.expr.genes, meth.expr.genes, meth.tfbs.genes))
  
  expr.gene.data.list <- lapply(total.genes, function(gene) {
    probe.ids <- unique(names(genes[genes == gene]))
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
  
  
  data.overview[snp,] <- c(length(meth.ids), length(meth.tfbs.genes), length(snp.expr.genes), length(snp.expr.no.gene.ids), 
                           length(meth.expr.genes), length(meth.expr.no.gene.ids), 10) #, dim(ggm.data)[2])

  # save(ggm.data, file = paste0(GGM.DIR, 'data/', snp, '.RData'))
    
  data[[snp]] <- ggm.data
}

save(data.overview, file = paste0(GGM.DIR, 'data.overview.RData'))




# # snps in meQTLs with snp and meth probe in herv
# snps <- intersect(as.character(unique(meqtl.pairs$snp)), as.character(unique(eqtl.pairs$snps)))
# eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('cis.', filter)]]
# for(snp in snps[1:5]) {
#   cat(paste0('Processing snp: ', snp), fill = T)
#   snp.tfbs.ids <- snp.tfbs.overlap$pairs[snp.tfbs.overlap$pairs$snp.id == snp, 'tfbs.id']
#   snp.tfbs.genes <- unique(snp.tfbs.overlap$tfbs.ranges[snp.tfbs.ids]$TF)
#   cat(paste0('tfs overlapping with snp: ', paste(snp.tfbs.genes, collapse=', ') ), fill = T)
#   
#   meth.ids <- as.character(meqtl.pairs[meqtl.pairs$snp == snp, 'cpg'])
#   cat(paste0('cpgs associated with snp: ', paste(meth.ids, collapse=', ') ), fill = T)
#   meth.tfbs.ids <- unique(meth.tfbs.overlap$pairs[meth.tfbs.overlap$pairs$meth.id %in% meth.ids, 'tfbs.id'])
#   meth.tfbs.genes <- unique(meth.tfbs.overlap$tfbs.ranges[meth.tfbs.ids]$TF)
#   cat(paste0('tfs overlapping with cpgs: ', paste(meth.tfbs.genes, collapse=', ') ), fill = T)
#   
#   expr.ids <- as.character(eqtl.pairs[eqtl.pairs$snps == snp, 'gene'])
#   cat(paste0('expr probes associated with snp: ', paste(expr.ids, sep=', ') ), fill = T)
#   expr.values <- expr.residuals[, expr.ids, drop = F]
#   expr.tfbs.ids <- unique(expr.tfbs.overlap$pairs[expr.tfbs.overlap$pairs$expr.id %in% expr.ids, 'tfbs.id'])
#   expr.tfbs.genes <- unique(expr.tfbs.overlap$tfbs.ranges[expr.tfbs.ids]$TF)
#   cat(paste0('tfs overlapping with expr probes: ', paste(expr.tfbs.genes, collapse=', ') ), fill = T)
#   
#   tfbs.genes <- unique(c(snp.tfbs.genes, meth.tfbs.genes, expr.tfbs.genes))
#   
#   tfbs.expr.values <- lapply(tfbs.genes, function(gene) {
#     probe.ids <- unique(names(genes[genes == gene]))
#     if (length(probe.ids) == 1) {
#       expr <- expr.residuals[, probe.ids[1]]
#     } else {
#       expr <- apply(expr.residuals[, probe.ids], 1, function(x) mean(x))
#     }
#     return(expr)
#   } )
#   tfbs.expr.values <- matrix(unlist(tfbs.expr.values), byrow=FALSE, ncol = length(tfbs.expr.values))
#   colnames(tfbs.expr.values) <- tfbs.genes
#   rownames(tfbs.expr.values) <- rownames(expr.residuals)
#   
#   meth.values <- get.residuals(meth.matrix, 'meth', meth.ids)
#   
#   snp.values <- get.snp.data(snp.ranges[snp])
#   
#   data.matrix <- cbind.data.frame(snp.values[id.map$axio_s4f4, , drop=F], meth.values[id.map$meth_f4, , drop=F], expr.values[id.map$expr_s4f4ogtt, , drop=F], tfbs.expr.values[id.map$expr_s4f4ogtt,])
#   rownames(data.matrix) <- id.map$expr_s4f4ogtt
#   
#   save(data.matrix, file = paste0(PATHS$DATA.DIR, 'ggm/', set, '/', snp, '.RData'))
# }
