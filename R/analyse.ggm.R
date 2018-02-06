source('Scripts/R/paths.R')
source('Scripts/R/lib.R')

require(illuminaHumanv3.db)
require(BDgraph)

set <- 'hervS1'
filter <- 'snp'
seed <- 'meqtl'
flanking <- 2.5e5
string  <- F

load(PATHS$HERV.EQTL.OVERLAP.DATA)
load(PATHS$HERV.METH.OVERLAP.DATA)
load(PATHS$HERV.EXPR.OVERLAP.DATA)
load(PATHS$EQTM.ME.DATA)

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
load(paste0(GGM.DIR, 'snps.RData'))
load(paste0(GGM.DIR, 'data.meta.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))


cis.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('cis.', filter)]][, 1:2]
trans.eqtl.pairs <- get(paste0(set, '.eqtl.overlap'))[[paste0('trans.', filter)]][, 1:2]

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



snp <- snps[1]



ggm.overview <- matrix(nrow = length(snps), ncol = 23)
ggm.meta <- list()
rownames(ggm.overview) <- snps
colnames(ggm.overview) <- c('entities', 'edges', 'ccs', 'cc.sizes', 'snp.cc.size',
                            'cpgs', 'con.cpgs', 'direct.cpgs', 
                            'herv.cpgs', 'con.herv.cpgs', 'cpg.genes', 'con.cpg.genes', 'herv.cpg.genes',
                            'snp.genes', 'con.snp.genes', 'herv.snp.genes', 'con.herv.snp.genes', 
                            'pos.eqtl', 'con.eqtl', 'direct.eqtl', 
                            'pos.eqtm', 'con.eqtm', 'direct.eqtm')

ggms <- list()
for (snp in snps) {
    cat(paste0('Processing snp: ', snp), fill = T)
    load(paste0(GGM.DIR, 'ggm/', snp, '.RData'))
    ggms[[snp]] <- ggm
    meta <- data.meta[[snp]]
    ggm.res <- graph.from.fit(ggm, NULL)
    ccs <- connComp(ggm.res)
    cc.lengths <- unlist(lapply(ccs, length))
    biggest.cc <- ccs[[which.max(cc.lengths)]]
    if (!snp %in% biggest.cc) {
      cat(paste0('SNP ', snp, ' not in biggest cc of ggm.'), fill = TRUE)
    }
    snp.cc <- ccs[[1]]
    
    entities <- nodes(ggm.res)
    edges <- numEdges(ggm.res)
    cpgs <- data.frame(row.names = meta$meth.ids)
    cpgs$con <- rownames(cpgs) %in% snp.cc
    cpgs$direct <- rownames(cpgs) %in% adj(ggm.res, snp)[[1]]
    cpgs$herv <- rownames(cpgs) %in% herv.cpgs
    
    
    cpg.genes <- data.frame(row.names = c(meta$meth.genes, meta$meth.no.gene.probes))
    cpg.genes$con <- rownames(cpg.genes) %in% snp.cc
    cpg.genes$herv <- rownames(cpg.genes) %in% herv.expr.ids | rownames(cpg.genes) %in% herv.genes
    
    snp.genes <- data.frame(row.names = c(meta$snp.genes, meta$snp.no.gene.probes))
    snp.genes$con <- rownames(snp.genes) %in% snp.cc
    snp.genes$herv <- rownames(snp.genes) %in% herv.expr.ids | rownames(snp.genes) %in% herv.genes
    
    eqtl.expr.ids <- unique(as.character(cis.eqtl.pairs[cis.eqtl.pairs$snps == snp, 'gene']))
    eqtl.genes <- unique(probe2gene[eqtl.expr.ids[eqtl.expr.ids %in% names(probe2gene)]])
    
    eqtls <- data.frame(row.names = rownames(snp.genes)[rownames(snp.genes) %in% eqtl.genes | rownames(snp.genes) %in% eqtl.expr.ids])
    eqtls$con <- rownames(eqtls) %in% snp.cc
    eqtls$direct <- rownames(eqtls) %in% adj(ggm.res, snp)[[1]]

    all.expr.ids <- unique(c(names(probe2gene[probe2gene %in% c(meta$tfbs.genes, meta$snp.genes, meta$meth.genes)]), meta$snp.no.gene.probes, meta$meth.no.gene.probes))
    
    eqtm.pairs <- all.eqtm.pairs[all.eqtm.pairs$cpg %in% rownames(cpgs) & all.eqtm.pairs$expr.id %in% all.expr.ids, ]
    eqtm.pairs$gene <- probe2gene[eqtm.pairs$expr.id]
    
    eqtm.pairs$con <- eqtm.pairs$cpg %in% snp.cc & (eqtm.pairs$expr.id %in% snp.cc | eqtm.pairs$gene %in% snp.cc)
    
    eqtm.pairs$direct <- apply(eqtm.pairs, 1, function(row) {
      cpg.edges <- adj(ggm.res, row[1])
      return(row[2] %in% cpg.edges | row[3] %in% cpg.edges)
    })
    #                        entities         edges0.5   edges0.9    ccs                   
    ggm.overview[snp,] <- c(length(entities), edges, length(ccs), paste(sort(cc.lengths[cc.lengths > 1], decreasing = T), collapse = ', '), length(snp.cc),
    #                       cpgs         con.cpgs           direct.cpgs 
                            dim(cpgs)[1], sum(cpgs$con), sum(cpgs$direct), sum(cpgs$herv), sum(cpgs$herv & cpgs$con), 
                            dim(cpg.genes)[1], sum(cpg.genes$con), sum(cpg.genes$herv), 
    #                       'snp.genes', 'con.snp.genes', 
                            dim(snp.genes)[1], sum(snp.genes$con), sum(snp.genes$herv), sum(snp.genes$herv & snp.genes$con),
    #                       'pos.eqtl', 'con.eqtl', 'direct.eqtl', 
                            dim(eqtls)[1], sum(eqtls$con), sum(eqtls$con), 
    #                       'pos.eqtm', 'con.eqtm', 'direct.eqtm'
                            dim(eqtm.pairs)[1], sum(eqtm.pairs$con), sum(eqtm.pairs$direct))
    
    ggm.meta[[snp]] <- list(ccs = ccs, cpgs = cpgs, cpg.genes = cpg.genes, snp.genes = snp.genes, eqtls = eqtls, eqtm.pairs = eqtm.pairs)
    
}

save(ggm.overview, file = paste0(GGM.DIR, 'ggm.overview.nocut.RData'))
save(ggm.meta, file = paste0(GGM.DIR, 'ggm.meta.nocut.RData'))