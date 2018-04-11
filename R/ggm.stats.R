source('Scripts/R/paths.R')
require(ggplot2)
require(gridExtra)

set <- 'hervS2'
filter <- 'meth'
seed <- 'meqtl'
flanking <- 2.5e5
string <- T
cutoff <- 0.9

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
load(paste0(GGM.DIR, 'ggm.overview.', cutoff,'.RData'))
if(filter == 'snp') {
  load(paste0(GGM.DIR, 'snps.RData'))
} else {
  load(paste0(GGM.DIR, 'cpgs.RData'))
}



plot.overview.bar <- function(overview.data, title, xlab, show.legend = T) {
  overview.flat <- melt(as.matrix(t(overview.data)))
  colnames(overview.flat) <- c('type', 'cpg.set', 'count')
  xtext <- gsub('\\|', '\n', rownames(overview.data))
  overview.bar <- ggplot(data=overview.flat, aes(x=cpg.set, y=count, fill=type)) + geom_bar(stat="identity")
  overview.bar <- overview.bar + labs(x = xlab, y = 'Entities', title=title)
  overview.bar <- overview.bar + scale_x_discrete(breaks=rownames(overview.data), labels=xtext)
  overview.bar <- overview.bar + guides(fill = ifelse(show.legend, 'legend', 'none'))
  overview.bar <- overview.bar + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
  return(overview.bar)
}

plot.blank.text <- function(text, title) {
  blank <- ggplot() + annotate(geom = "text", x = 50 , y = 50, label = text) + geom_point() + xlim(0, 100) + ylim(0, 100) 
  blank <- blank + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(Title=title)
  return(blank)
}

### ggm data overview

### cpg data
GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.meth.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')

load(paste0(GGM.DIR, 'data.meta.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))

filtered.cpg.data.overview <- data.overview[(data.overview$seed.meth.tfs > 0 | data.overview$extra.meth.tfs | data.overview$path.genes)
                                        & (data.overview$seed.meth.genes + data.overview$seed.meth.no.gene.probes + data.overview$extra.meth.genes + data.overview$extra.meth.no.gene.probes) > 0
                                        & (data.overview$snp.genes + data.overview$snp.no.gene.probes) > 0
                                        ,]
filtered.cpg.data.overview <- filtered.cpg.data.overview[order(filtered.cpg.data.overview$total.entities, decreasing = T),]

filtered.cpg.sets <- rownames(filtered.cpg.data.overview)

save(filtered.cpg.sets, file = paste0(GGM.DIR, 'filtered.cpgs.RData'))
save(filtered.cpg.data.overview, file = paste0(GGM.DIR, 'data.overview.filtered.RData'))

get.cpg.data.overview.reduced <- function(data.overview) {
  data.overview.rowSums <- rowSums(data.overview[, -1])
  
  data.overview.reduced <- data.frame(matrix(nrow = nrow(data.overview), ncol = 9))
  rownames(data.overview.reduced) <- rownames(data.overview)
  colnames(data.overview.reduced) <- c('total.entites', 'cpgs', 'snps', 'tfs', 'snp.genes', 'snp.probes', 'meth.genes', 'meth.probes', 'path.genes')
  
  for(i in 1:nrow(data.overview)) {
    if(data.overview[i , 'total.entities'] == data.overview.rowSums[i]) {
      tfs <- data.overview$seed.meth.tfs[i] +  data.overview$extra.meth.tfs[i]
      meth.genes <- data.overview$seed.meth.genes[i] + data.overview$extra.meth.genes[i]
      meth.probes <- data.overview$seed.meth.no.gene.probes[i] + data.overview$extra.meth.no.gene.probes[i]
      snp.genes <- data.overview$snp.genes[i]
      snp.probes <- data.overview$snp.no.gene.probes[i]
      path.genes <- data.overview$path.genes[i]
    } else {  
      meta <- data.meta[[rownames(data.overview.reduced)[i]]]
      meth.gene.ids <- unique(c(meta$seed.meth.genes, meta$extra.meth.genes))
      meth.genes <- length(meth.gene.ids)
      meth.probe.ids <- unique(c(meta$seed.meth.no.gene.probes, meta$extra.meth.no.gene.probes))
      meth.probes <- length(meth.probe.ids)
      snp.gene.ids <- meta$snp.genes[!meta$snp.genes %in% meth.gene.ids]
      snp.genes <- length(snp.gene.ids)
      snp.probes <- length(meta$snp.no.gene.probes[!meta$snp.no.gene.probes %in% meth.probe.ids])
      tf.gene.ids <- unique(c(meta$seed.tfbs.genes, meta$extra.tfbs.genes))
      tf.gene.ids <- tf.gene.ids[!tf.gene.ids %in% c(meth.gene.ids, snp.gene.ids)]
      tfs <- length(tf.gene.ids)
      path.gene.ids <- meta$path.genes[meta$path.genes %in% c(meth.gene.ids, snp.gene.ids, tf.gene.ids)]
      path.genes <- length(path.gene.ids)
    }
    data.overview.reduced[i, ] <- c(data.overview$total.entities[i], 
                                    data.overview$seed.cpgs[i] + data.overview$extra.cpgs[i], 
                                    data.overview$best.snps[i],
                                    tfs, 
                                    snp.genes,
                                    snp.probes, 
                                    meth.genes,
                                    meth.probes,
                                    path.genes)
    
  }
  
  data.overview.reduced <- data.overview.reduced[order(data.overview.reduced[, 'total.entites'], decreasing = T),]
  
  return(data.overview.reduced)
}

cpg.data.overview.reduced <- get.cpg.data.overview.reduced(filtered.cpg.data.overview)
View(filtered.cpg.ggm.overview[, c(1:2, 4:18, 21:23, 26:28)])


### snp data
GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', 'snp', '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')

load(paste0(GGM.DIR, 'data.meta.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))


filtered.snp.data.overview <- data.overview[(data.overview$meth.genes + data.overview$meth.no.gene.probes) > 0
                                            & (data.overview$snp.genes + data.overview$snp.no.gene.probes) > 0
                                            & (data.overview$TFs + data.overview$path.genes) > 0,]

filtered.snps <- rownames(filtered.snp.data.overview)
save(filtered.snps, file = paste0(GGM.DIR, 'filtered.snps.RData'))

filtered.snp.data.overview <- filtered.snp.data.overview[order(filtered.snp.data.overview$total.entities, decreasing = T), c(8, 1:7)]
save(filtered.snp.data.overview, file = paste0(GGM.DIR, 'data.overview.filtered.RData'))

get.snp.data.overview.reduced <- function(data.overview) {
  data.overview.rowSums <- rowSums(data.overview[, -1])
  
  data.overview.reduced <- data.frame(matrix(nrow = nrow(data.overview), ncol = 9))
  rownames(data.overview.reduced) <- rownames(data.overview)
  colnames(data.overview.reduced) <- c('total.entites', 'cpgs', 'snps', 'tfs', 'snp.genes', 'snp.probes', 'meth.genes', 'meth.probes', 'path.genes')
  
  for(snp in rownames(data.overview)) {
    row <- data.overview[snp, , drop = F]
    if(row$total.entities == data.overview.rowSums[snp]) {
      tfs <- row$TFs
      meth.genes <- row$meth.genes
      meth.probes <- row$meth.no.gene.probes
      snp.genes <- row$snp.genes
      snp.probes <- row$snp.no.gene.probes
      path.genes <- row$path.genes
    } else {  
      meta <- data.meta[[snp]]
      meth.genes <- length(meta$meth.genes)
      meth.probes <- length(meta$meth.no.gene.probes)
      snp.genes <- length(meta$snp.genes[!meta$snp.genes %in% meta$meth.genes])
      snp.probes <- length(meta$snp.no.gene.probes[!meta$snp.no.gene.probes %in% meta$meth.no.gene.probes])
      tfs <- length(meta$tfbs.genes[!meta$tfbs.genes %in% c( meta$meth.genes, meta$snp.genes)])
      path.genes <- length(meta$path.genes)
    }
    data.overview.reduced[snp, ] <- c(row$total.entities, 
                                    row$cpgs, 
                                    1,
                                    tfs, 
                                    snp.genes,
                                    snp.probes, 
                                    meth.genes,
                                    meth.probes,
                                    path.genes)
    
  }
  
  
  return(data.overview.reduced)
}

snp.data.overview.reduced <- get.snp.data.overview.reduced(filtered.snp.data.overview)



snp.ggm.data.top.bar <- plot.overview.bar(snp.data.overview.reduced[1:10, -1], 'A', 'Seed SNP', F)
snp.ggm.data.bot.bar <- plot.overview.bar(snp.data.overview.reduced[(nrow(snp.data.overview.reduced)-9):nrow(snp.data.overview.reduced), -1], 'B', 'Seed SNP')
cpg.ggm.data.top.bar <- plot.overview.bar(cpg.data.overview.reduced[1:10, -1], 'C', 'Seed CpG set', F)
cpg.ggm.data.bot.bar <- plot.overview.bar(cpg.data.overview.reduced[(nrow(cpg.data.overview.reduced)-9):nrow(cpg.data.overview.reduced), -1], 'D', 'Seed CpG set', T)

pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_ggm_data_entity_bar.pdf'), width = 6.3, height = 6)
grid.arrange(snp.ggm.data.top.bar, snp.ggm.data.bot.bar, cpg.ggm.data.top.bar, cpg.ggm.data.bot.bar, ncol = 2, widths = c(0.41, 0.59))
dev.off()




### calculated ggms
### cpg data
GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/hervS2.meqtl.meth.250kb.string/')

load(paste0(GGM.DIR, 'ggm.overview.', cutoff, '.RData'))
load(paste0(GGM.DIR, 'ggm.meta', cutoff, '.RData'))


filtered.cpg.ggm.overview <- ggm.overview[ggm.overview$con.seed.cpgs > 0
                                      & ggm.overview$con.snps > 0
                                      & (ggm.overview$con.seed.cpg.genes + ggm.overview$con.extra.cpg.genes) > 0
                                      & ggm.overview$con.snp.genes > 0
                                      & (ggm.overview$con.tfs > 0 | ggm.overview$con.path.genes > 0)
                                      ,]

save(filtered.cpg.ggm.overview, file = paste0(GGM.DIR, 'filtered.ggm.overview.', cutoff, '.RData'))

get.cpg.ggm.overview.reduced <- function(ggm.overview) {
  ggm.overview.reduced <- data.frame(matrix(nrow = nrow(ggm.overview), ncol = 8))
  rownames(ggm.overview.reduced) <- rownames(ggm.overview)
  colnames(ggm.overview.reduced) <- c('nodes', 'edges', 'cpgs', 'snps', 'cpg.genes', 'snp.genes', 'tfs', 'path.genes')
  
  biggest.cc.sum <- with(ggm.overview, con.seed.cpgs + extra.con.cpgs + con.snps + con.tfs + con.path.genes + con.seed.cpg.genes + con.extra.cpg.genes + con.snp.genes)
  names(biggest.cc.sum) <- rownames(ggm.overview)
  
  for(cpg.set in rownames(ggm.overview)) {
    cat(cpg.set, fill = T)
    row <- ggm.overview[cpg.set,]
    biggest.cc.size <- as.numeric(strsplit(row$cc.sizes, split = ', ')[[1]][1])
    if(biggest.cc.size == biggest.cc.sum[cpg.set]) {
      cpg.genes <- row$con.seed.cpg.genes + row$con.extra.cpg.genes
      snp.genes <- row$con.snp.genes
      tfs <- row$con.tfs
      path.genes <- row$con.path.genes
    } else {
      meta <- ggm.meta[[cpg.set]]
      
      cpg.gene.ids <- unique(c(rownames(meta$seed.cpg.genes)[meta$seed.cpg.genes$con], rownames(meta$extra.cpg.genes)[meta$extra.cpg.genes$con]))
      cpg.genes <- length(cpg.gene.ids)
      snp.gene.ids <- rownames(meta$snp.genes)[meta$snp.genes$con]
      snp.genes <- length(snp.gene.ids[!snp.gene.ids %in% cpg.gene.ids])
      tf.gene.ids <- rownames(meta$tfs)[meta$tfs$con]
      tf.gene.ids <- tf.gene.ids[!tf.gene.ids %in% c(cpg.gene.ids, snp.gene.ids)]
      tfs <- length(tf.gene.ids)
      path.gene.ids <- rownames(meta$path.genes)[meta$path.genes$con]
      path.genes <- length(path.gene.ids[!path.gene.ids %in% c(cpg.gene.ids, snp.gene.ids, tf.gene.ids)])
    }
    ggm.overview.reduced[cpg.set, ] <- c(biggest.cc.size,
                                                  row$edges,
                                                  row$con.seed.cpgs + row$extra.con.cpgs,
                                                  row$con.snps,
                                                  cpg.genes,
                                                  snp.genes,
                                                  tfs,
                                                  path.genes)
  }    
  return(ggm.overview.reduced)
}

cpg.ggm.overview.reduced <- get.cpg.ggm.overview.reduced(filtered.cpg.ggm.overview)


### snp ggms
GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', 'snp', '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')

load(paste0(GGM.DIR, 'ggm.overview.', cutoff, '.RData'))


filtered.snp.ggm.overview <- ggm.overview[ggm.overview$con.cpgs > 0
                                          & ggm.overview$con.snp > 0
                                          & ggm.overview$con.cpg.genes > 0
                                          & ggm.overview$con.snp.genes > 0
                                          & (ggm.overview$con.tfs > 0 | ggm.overview$con.path.genes > 0)
                                          ,]
filtered.snp.ggm.overview <- filtered.snp.ggm.overview[order(filtered.snp.ggm.overview$entities, decreasing = T), ]
                                                       
save(filtered.snp.ggm.overview, file = paste0(GGM.DIR, 'filtered.ggm.overview.', cutoff, '.RData'))
                                                       

get.snp.ggm.overview.reduced <- function(ggm.overview) {
  ggm.overview.reduced <- data.frame(matrix(nrow = nrow(ggm.overview), ncol = 8))
  rownames(ggm.overview.reduced) <- rownames(ggm.overview)
  colnames(ggm.overview.reduced) <- c('nodes', 'edges', 'cpgs', 'snps', 'cpg.genes', 'snp.genes', 'tfs', 'path.genes')
  
  con.sum <- with(ggm.overview, con.snp + con.cpgs + con.snp.genes + con.cpg.genes + con.tfs + con.path.genes)
  names(con.sum) <- rownames(ggm.overview)
  
  for(snp in rownames(ggm.overview)) {
    cat(snp, fill = T)
    row <- ggm.overview[snp,]
    biggest.cc.size <- as.numeric(strsplit(row$cc.sizes, split = ', ')[[1]][1])
    if(biggest.cc.size == con.sum[snp]) {
      cpg.genes <- row$con.cpg.genes
      snp.genes <- row$con.snp.genes
      tfs <- row$con.tfs
      path.genes <- row$con.path.genes
    } else {
      meta <- ggm.meta[[snp]]
      
      cpg.gene.ids <- rownames(meta$cpg.genes)[meta$cpg.genes$con]
      cpg.genes <- length(cpg.gene.ids)
      snp.gene.ids <- rownames(meta$snp.genes)[meta$snp.genes$con]
      snp.genes <- length(snp.gene.ids[!snp.gene.ids %in% cpg.gene.ids])
      tf.gene.ids <- rownames(meta$tfs)[meta$tfs$con]
      tf.gene.ids <- tf.gene.ids[!tf.gene.ids %in% c(cpg.gene.ids, snp.gene.ids)]
      tfs <- length(tf.gene.ids)
      path.gene.ids <- rownames(meta$path.genes)[meta$path.genes$con]
      path.genes <- length(path.gene.ids[!path.gene.ids %in% c(cpg.gene.ids, snp.gene.ids, tf.gene.ids)])
    }
    ggm.overview.reduced[snp, ] <- c(biggest.cc.size,
                                         row$edges,
                                         row$con.cpgs,
                                         row$con.snp,
                                         cpg.genes,
                                         snp.genes,
                                         tfs,
                                         path.genes)
  }    
  ggm.overview.reduced <- ggm.overview.reduced[order(ggm.overview.reduced$nodes, decreasing = T),]
  return(ggm.overview.reduced)
}

snp.ggm.overview.reduced <- get.snp.ggm.overview.reduced(filtered.snp.ggm.overview)


snp.ggm.top.entity.bar <- plot.overview.bar(snp.ggm.overview.reduced[1:10, -(1:2)], 'A', 'Seed SNP', F)
snp.ggm.bot.entity.bar <- plot.overview.bar(snp.ggm.overview.reduced[(nrow(snp.ggm.overview.reduced)-9):nrow(snp.ggm.overview.reduced), -(1:2)], 'B', 'Seed SNP', T)
cpg.ggm.entity.bar <- plot.overview.bar(cpg.ggm.overview.reduced[,-(1:2)], 'C', 'Seed CpG set')

pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_ggm_cc_nodes_bar.pdf'), width = 6.3, height = 6)
grid.arrange(snp.ggm.top.entity.bar, snp.ggm.bot.entity.bar, cpg.ggm.entity.bar, layout_matrix = rbind(c(1, 2), c(3, 3)), widths = c(0.41, 0.59))
dev.off()

plot.entity.cc.size <- function(ggm.overview, filtered.ggm.overview) {
  df <- data.frame(entities = ggm.overview$entities, biggest.cc = as.numeric(sapply(ggm.overview$cc.sizes, function(x) return(strsplit(x, split = ',')[[1]][1]))), retained = 'no')
  df[rownames(filtered.ggm.overview), 'retained'] <- 'yes'
  df[is.na(df$biggest.cc), 'biggest.cc'] <- 1
  g <- ggplot(df, aes(x=entities, y=biggest.cc, group = retained)) + geom_point(aes(color = retained), size = 0.5)
}

get.ggm.overview.table <- function(ggm.overview.reduced) {
  table <- data.frame(matrix(nrow = 8, ncol =4))
  rownames(table) <- colnames(ggm.overview.reduced)
  colnames(table) <- c('Min', 'Max', 'Median', 'Mean')
  for(x in rownames(table)) {
    table[x,] <- c(min(ggm.overview.reduced[, x]),
                   max(ggm.overview.reduced[, x]),
                   median(ggm.overview.reduced[, x]),
                   round(mean(ggm.overview.reduced[, x]), digits = 2))
  }
  table$Type <- c('Nodes', 'Edges', 'CpGs', 'SNPs', 'CpG genes', 'SNP genes', 'TFs', 'Path genes')
  table <- table[,c(5, 1:4)]
  return(table)
}

snp.ggm.overview.table <- get.ggm.overview.table(snp.ggm.overview.reduced)
cpg.ggm.overview.table <- get.ggm.overview.table(cpg.ggm.overview.reduced)

write.table(snp.ggm.overview.table, file = paste0(PATHS$TABLE.DIR, 'snp.ggm.overview.table.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)
write.table(cpg.ggm.overview.table, file = paste0(PATHS$TABLE.DIR, 'cpg.ggm.overview.table.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)
