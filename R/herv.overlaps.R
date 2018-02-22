source('Scripts/R/paths.R')
source('Scripts/R/go.enrichment.R')

require(GenomicRanges)
require(illuminaHumanv3.db)
require(FDb.InfiniumMethylation.hg19)
require(rtracklayer)
require(data.table)

if(!file.exists(PATHS$EXPR.RANGES.DATA)) {
  get.expression.ranges <- function () {
    require(illuminaHumanv3.db)
    allLocs <- unlist(as.list(illuminaHumanv3GENOMICLOCATION))
    start <- as.numeric(unlist(lapply(allLocs , function(x)
      strsplit(as.character(x),":")[[1]][2])));
    validPos <- !is.na(start);
    start <- start[validPos];
    
    chrs <- unlist(lapply(allLocs, function(x)
      strsplit(as.character(x),":")[[1]][1]))[validPos];
    end <- as.numeric(unlist(lapply(allLocs , function(x)
      strsplit(as.character(x),":")[[1]][3])))[validPos];
    strand <- substr(unlist(lapply(allLocs , function(x)
      strsplit(as.character(x),":")[[1]][4])), 1, 1)[validPos];
    ids <- names(allLocs)[validPos];
    gr <- GRanges(chrs, ranges=IRanges(start,end), strand=strand)
    return(gr)
  }
  expr.ranges <- get.expression.ranges()
  save(expr.ranges, file = PATHS$EXPR.RANGES.DATA)
} else {
  load(PATHS$EXPR.RANGES)
}

if (!file.exists(PATHS$METH.RANGES.DATA)) {
  meth.ranges <- getPlatform(platform = 'HM450', genome = 'hg19')
  meth.ranges <- meth.ranges[grep('cg|ch', meth.ranges$probeType)]
  save(meth.ranges, file = PATHS$METH.RANGES.DATA)
} else {
  load(PATHS$METH.RANGES.DATA)
}

if(!file.exists(PATHS$SNP.RANGES.DATA)) {
  snp.info <- fread(PATHS$F.SNP.POS, sep = '\t', header = T)
  snp.ranges <- GRanges(seqnames=snp.info$chr, ranges=IRanges(start=snp.info$pos, width=1), from=snp.info$from, to=snp.info$to)
  names(snp.ranges) <- snp.info$snp
  save(snp.ranges, file = PATHS$SNP.RANGES.DATA)
} else {
  load(PATHS$SNP.RANGES.DATA)
}
load(PATHS$TFBS.RANGES.DATA)

calc.overlap.data <- function (herv.ranges, essay.ranges, data.type) {
  overlap.hits <- findOverlaps(herv.ranges, essay.ranges, type = 'any')
  herv.overlap.ranges <- herv.ranges[unique(queryHits(overlap.hits))]
  essay.overlap.ranges <- essay.ranges[unique(subjectHits(overlap.hits))]
  out <- list()
  out$pairs <- cbind.data.frame(names(herv.ranges[queryHits(overlap.hits)]), names(essay.ranges[subjectHits(overlap.hits)]), stringsAsFactors = FALSE)
  colnames(out$pairs) <- c('herv.id', paste0(data.type, '.id'))
  out$herv.ranges <- herv.overlap.ranges
  out[[paste0(data.type, '.ranges')]] <- essay.overlap.ranges
  return(out)
}

combine.overlaps <- function (overlap1, overlap2, overlap1.type, overlap2.type) {
  herv.ids <- intersect(overlap1$pairs$herv.id, overlap2$pairs$herv.id)
  out <- list()
  out$herv.ranges <- overlap1$herv.ranges[herv.ids]
  overlap1.pairs <- overlap1$pairs[overlap1$pairs$herv.id %in% herv.ids, ]
  overlap2.pairs <- overlap2$pairs[overlap2$pairs$herv.id %in% herv.ids, ]
  out$pairs <- data.frame(matrix(nrow = 0, ncol = 3))
  for(herv.id in herv.ids) {
    for(essay1.id in overlap1.pairs[overlap1.pairs$herv.id == herv.id, 2]) {
      for(essay2.id in overlap2.pairs[overlap2.pairs$herv.id == herv.id, 2]) {
        out$pairs <- rbind.data.frame(out$pairs, c(herv.id, essay1.id, essay2.id), stringsAsFactors = FALSE)
      }
    }
  }
  colnames(out$pairs) <- c('herv.id', paste0(overlap1.type, '.id'), paste0(overlap2.type, '.id'))
  
  out[[paste0(overlap1.type, '.ranges')]] <- overlap1[[paste0(overlap1.type, '.ranges')]][unique(overlap1.pairs[, 2])]
  out[[paste0(overlap2.type, '.ranges')]] <- overlap1[[paste0(overlap1.type, '.ranges')]][unique(overlap1.pairs[, 2])]
  return (out)
}

print.overlap.info <- function(overlap, overlap.type) {
  cat(paste0("Overlap info:\n# Overlaps: ", dim(overlap$pairs)[1], 
             "\n# hERVs: ", length(overlap$herv.ranges),
             "\n# probes: ", length(overlap[[paste0(overlap.type, '.ranges')]])), fill = T)
}


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

name.ranges.by.coordinates <- function(ranges) {
  names <- paste0(as.character(seqnames(ranges)), ':', start(ranges), '-', end(ranges))
  names(ranges) <- names
  return (ranges)
}

merge.overlapping.ranges <- function(ranges) {
  merged.ranges <- reduce(ranges)
  merged.ranges <- name.ranges.by.coordinates(merged.ranges)
  merged.ranges <- sort(merged.ranges, ignore.strand = T)
  
  overlaps <- findOverlaps(merged.ranges, ranges)
  
  df <- cbind.data.frame(queryHits(overlaps), ranges$name, ranges$score, stringsAsFactors = F)
  colnames(df) <- c('index', 'name', 'score')
  dt <- data.table(df, key='index')
  dt2 <- dt[, list(name=paste(unique(name), collapse='|'),
                   score=sum(score)), by=key(dt)]
  
  merged.ranges$name <- dt2[, name]
  merged.ranges$score <- dt2[, score]
  
  return(merged.ranges)
}

if(!file.exists(PATHS$HERV.RANGES.DATA)) {
  hervS1.ranges <- import(PATHS$HERVS1.ANNOT, format = 'BED')
  strand(hervS1.ranges) <- '*'
  hervS1.ranges <- name.ranges.by.coordinates(hervS1.ranges)
  hervS2.ranges <- import(PATHS$HERVS2.ANNOT, format = 'BED')
  strand(hervS2.ranges) <- '*'
  hervS2.ranges <- name.ranges.by.coordinates(hervS2.ranges)
  hervS3.ranges <- import(PATHS$HERVS3.ANNOT, format = 'BED')
  strand(hervS3.ranges) <- '*'
  hervS3.ranges <- name.ranges.by.coordinates(hervS3.ranges)
  
  save(
    hervS1.ranges,
    hervS2.ranges,
    hervS3.ranges,
    file = PATHS$HERV.RAW.RANGES.DATA)
  
  hervS1.ranges <- merge.overlapping.ranges(hervS1.ranges)
  hervS2.ranges <- merge.overlapping.ranges(hervS2.ranges)
  hervS3.ranges <- merge.overlapping.ranges(hervS3.ranges)
  
  hervS1.1kb.ranges <- enlarge.ranges(hervS1.ranges, 1000)
  hervS2.1kb.ranges <- enlarge.ranges(hervS2.ranges, 1000)
  hervS3.1kb.ranges <- enlarge.ranges(hervS3.ranges, 1000)
  
  hervS1.2kb.ranges <- enlarge.ranges(hervS1.ranges, 2000)
  hervS2.2kb.ranges <- enlarge.ranges(hervS2.ranges, 2000)
  hervS3.2kb.ranges <- enlarge.ranges(hervS3.ranges, 2000)
  
  save(
    hervS1.ranges,
    hervS2.ranges,
    hervS3.ranges,
    hervS1.1kb.ranges,
    hervS2.1kb.ranges,
    hervS3.1kb.ranges,
    hervS1.2kb.ranges,
    hervS2.2kb.ranges,
    hervS3.2kb.ranges,
    file = PATHS$HERV.RANGES.DATA
  )
  save(hervS2.ranges, file = PATHS$HERVS2.RANGES.DATA)
  save(hervS2.2kb.ranges, file = PATHS$HERVS2.2KB.RANGES.DATA)
} else {
  load(PATHS$HERV.RANGES.DATA)
}

load(PATHS$EXPR.DATA)
expr.data <- f4.norm

load(PATHS$METH.DATA)
meth.data <- data.frame(beta)
rm(beta)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    expr.overlap.name <- paste0('herv', set, flanking, '.expr.overlap')
    assign(expr.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), expr.ranges, 'expr'))
    meth.overlap.name <- paste0('herv', set, flanking, '.meth.overlap')
    assign(meth.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), meth.ranges, 'meth'))
    snp.overlap.name <- paste0('herv', set, flanking, '.snp.overlap')
    assign(snp.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), snp.ranges, 'snp'))
    tfbs.overlap.name <- paste0('herv', set, flanking, '.tfbs.overlap')
    assign(tfbs.overlap.name, calc.overlap.data(get(paste0('herv', set, flanking, '.ranges')), blood.tfbs.ranges, 'tfbs'))
  }
}

if(!file.exists(PATHS$EXPR.GENE.ANNOT.DATA)) {
  probe2gene <- unlist(as.list(illuminaHumanv3SYMBOL))
  probe2gene <- probe2gene[names(expr.ranges)]
  probe2gene <- probe2gene[!is.na(probe2gene)]
  save(probe2gene, file = PATHS$EXPR.GENE.ANNOT.DATA)
} else {
  load(PATHS$EXPR.GENE.ANNOT.DATA)
}

gene.universe <- unique(probe2gene)

expr.overlap.overview <- data.frame(matrix(nrow = 4, ncol = 10))
rownames(expr.overlap.overview) <- c('Pairs', 'HERVs', 'Probes', 'Genes')
sets <- c('S1', 'S2', 'S3')
flanking.distances <- c('', '.1kb', '.2kb')
colnames(expr.overlap.overview) <- c('Set', as.vector(t(outer(sets, flanking.distances, paste0))))
expr.overlap.overview$Set <- rownames(expr.overlap.overview)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    expr.overlap.name <- paste0('herv', set, flanking, '.expr.overlap')
    cat(paste0('Processing ', expr.overlap.name), fill = T)
    expr.overlap <- get(expr.overlap.name)
    expr.overlap.overview['Pairs', paste0(set, flanking)] <- dim(expr.overlap$pairs)[1]
    expr.overlap.overview['HERVs', paste0(set, flanking)] <- length(expr.overlap$herv.ranges)
    expr.overlap.overview['Probes', paste0(set, flanking)] <- length(expr.overlap$expr.ranges)
    expr.overlap.genes <- unique(probe2gene[names(expr.overlap$expr.ranges)[names(expr.overlap$expr.ranges) %in% names(probe2gene)]])
    expr.overlap.overview['Genes', paste0(set, flanking)] <- length(expr.overlap.genes)
    assign(paste0('herv', set, flanking, '.expr.enrichment'), go.enrichment(expr.overlap.genes, gene.universe, gsc, c('BP')))
  }
}

write.table(expr.overlap.overview, file = paste0(PATHS$TABLE.DIR, 'expr.overlap.overview.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)

hervS2.2kb.expr.enrichment.signigicant <- hervS2.2kb.expr.enrichment[1:6, c(2, 8, 3, 9)]
colnames(hervS2.2kb.expr.enrichment.signigicant) <- c('Term ID', 'Term', 'p', 'fdr')
write.table(hervS2.2kb.expr.enrichment.signigicant, file = paste0(PATHS$TABLE.DIR, 'hervS2.2kb.expr.enrichment.tsv'), 
              quote = F, sep = '\t', row.names = F, col.names = T)

save(
  hervS1.expr.enrichment,
  hervS2.expr.enrichment,
  hervS3.expr.enrichment,
  hervS1.1kb.expr.enrichment,
  hervS2.1kb.expr.enrichment,
  hervS3.1kb.expr.enrichment,
  hervS1.2kb.expr.enrichment,
  hervS2.2kb.expr.enrichment,
  hervS3.2kb.expr.enrichment,
  file = PATHS$HERV.EXPR.ENRICHMENT.DATA
)

meth.overlap.overview <- data.frame(matrix(nrow = 3, ncol = 10))
rownames(meth.overlap.overview) <- c('Pairs', 'HERVs', 'CpGs')
sets <- c('S1', 'S2', 'S3')
flanking.distances <- c('', '.1kb', '.2kb')
colnames(meth.overlap.overview) <- c('Set', as.vector(t(outer(sets, flanking.distances, paste0))))
meth.overlap.overview$Set <- rownames(meth.overlap.overview)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    meth.overlap.name <- paste0('herv', set, flanking, '.meth.overlap')
    cat(paste0('Processing ', meth.overlap.name), fill = T)
    meth.overlap <- get(meth.overlap.name)
    meth.overlap.overview['Pairs', paste0(set, flanking)] <- dim(meth.overlap$pairs)[1]
    meth.overlap.overview['HERVs', paste0(set, flanking)] <- length(meth.overlap$herv.ranges)
    meth.overlap.overview['CpGs', paste0(set, flanking)] <- length(meth.overlap$meth.ranges)
  }
}

write.table(meth.overlap.overview, file = paste0(PATHS$TABLE.DIR, 'meth.overlap.overview.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)

snp.overlap.overview <- data.frame(matrix(nrow = 3, ncol = 10))
rownames(snp.overlap.overview) <- c('Pairs', 'HERVs', 'SNPs')
sets <- c('S1', 'S2', 'S3')
flanking.distances <- c('', '.1kb', '.2kb')
colnames(snp.overlap.overview) <- c('Set', as.vector(t(outer(sets, flanking.distances, paste0))))
snp.overlap.overview$Set <- rownames(snp.overlap.overview)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    snp.overlap.name <- paste0('herv', set, flanking, '.snp.overlap')
    cat(paste0('Processing ', snp.overlap.name), fill = T)
    snp.overlap <- get(snp.overlap.name)
    snp.overlap.overview['Pairs', paste0(set, flanking)] <- dim(snp.overlap$pairs)[1]
    snp.overlap.overview['HERVs', paste0(set, flanking)] <- length(snp.overlap$herv.ranges)
    snp.overlap.overview['SNPs', paste0(set, flanking)] <- length(snp.overlap$snp.ranges)
  }
}

write.table(snp.overlap.overview, file = paste0(PATHS$TABLE.DIR, 'snp.overlap.overview.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)




for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    for (type in c('expr', 'meth', 'snp', 'tfbs')) {
      overlap.name <- paste0('herv', set, flanking, '.', type, '.overlap')
      cat(overlap.name, fill = T)
      print.overlap.info(get(overlap.name), type)
    }
  }
}

save(
  hervS1.expr.overlap,
  hervS2.expr.overlap,
  hervS3.expr.overlap,
  hervS1.1kb.expr.overlap,
  hervS2.1kb.expr.overlap,
  hervS3.1kb.expr.overlap,
  hervS1.2kb.expr.overlap,
  hervS2.2kb.expr.overlap,
  hervS3.2kb.expr.overlap,
  file = PATHS$HERV.EXPR.OVERLAP.DATA
)

save(
  hervS1.meth.overlap,
  hervS2.meth.overlap,
  hervS3.meth.overlap,
  hervS1.1kb.meth.overlap,
  hervS2.1kb.meth.overlap,
  hervS3.1kb.meth.overlap,
  hervS1.2kb.meth.overlap,
  hervS2.2kb.meth.overlap,
  hervS3.2kb.meth.overlap,
  file = PATHS$HERV.METH.OVERLAP.DATA
)

save(
  hervS1.snp.overlap,
  hervS2.snp.overlap,
  hervS3.snp.overlap,
  hervS1.1kb.snp.overlap,
  hervS2.1kb.snp.overlap,
  hervS3.1kb.snp.overlap,
  hervS1.2kb.snp.overlap,
  hervS2.2kb.snp.overlap,
  hervS3.2kb.snp.overlap,
  file = PATHS$HERV.SNP.OVERLAP.DATA
)

save(
  hervS1.tfbs.overlap,
  hervS2.tfbs.overlap,
  hervS3.tfbs.overlap,
  hervS1.1kb.tfbs.overlap,
  hervS2.1kb.tfbs.overlap,
  hervS3.1kb.tfbs.overlap,
  hervS1.2kb.tfbs.overlap,
  hervS2.2kb.tfbs.overlap,
  hervS3.2kb.tfbs.overlap,
  file = PATHS$HERV.TFBS.OVERLAP.DATA
)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    overlap.name <- paste0('herv', set, flanking, '.both.overlap') 
    cat(overlap.name, fill = TRUE)
    assign(overlap.name, combine.overlaps(get(paste0('herv', set, flanking, '.expr.overlap')),
                                          get(paste0('herv', set, flanking, '.meth.overlap')),
                                          'expr', 'meth'))
  }
}
