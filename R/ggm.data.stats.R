source('Scripts/R/paths.R')

require(ggplot2)
require(gridExtra)
require(reshape2)

set <- 'hervS2'
filter <- 'snp'
seed <- 'meqtl'
flanking <- 2.5e5
string <- T

GGM.DIR <- paste0(PATHS$DATA.DIR, 'ggm/', set, '.', seed, '.', filter, '.', flanking/1000, 'kb', ifelse(string, '.string', ''), '/')
load(paste0(GGM.DIR, 'cpgs.RData'))

if(filter == 'snp') {
  load(paste0(GGM.DIR, 'snps.RData'))
} else {
  load(paste0(GGM.DIR, 'cpgs.RData'))
}

load(paste0(GGM.DIR, 'data.meta.RData'))
load(paste0(GGM.DIR, 'data.overview.RData'))

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
    } else {  
      meta <- data.meta[[rownames(data.overview.reduced)[i]]]
      meth.gene.ids <- unique(c(meta$seed.meth.genes, meta$extra.meth.genes))
      meth.genes <- length(meth.gene.ids)
      meth.probe.ids <- unique(c(meta$seed.meth.no.gene.probes, meta$extra.meth.no.gene.probes))
      meth.probes <- length(meth.probe.ids)
      snp.genes <- length(meta$snp.genes[!meta$snp.genes %in% meth.gene.ids])
      snp.probes <- length(meta$snp.no.gene.probes[!meta$snp.no.gene.probes %in% meth.probe.ids])
      tf.gene.ids <- unique(c(meta$seed.tfbs.genes, meta$extra.tfbs.genes))
      tf.gene.ids <- tf.gene.ids[!tf.gene.ids %in% c(meth.gene.ids, meta$snp.genes)]
      tfs <- length(tf.gene.ids)
    }
      data.overview.reduced[i, ] <- c(data.overview$total.entities[i], 
                                      data.overview$seed.cpgs[i] + data.overview$extra.cpgs[i], 
                                      data.overview$best.snps[i],
                                      tfs, 
                                      snp.genes,
                                      snp.probes, 
                                      meth.genes,
                                      meth.probes,
                                      data.overview$path.genes[i])
    
  }
  
  data.overview.reduced <- data.overview.reduced[order(data.overview.reduced[, 'total.entites'], decreasing = T),]
  
  data.overview.reduced <- data.overview.reduced[(data.overview.reduced$tfs > 0 | data.overview.reduced$path.genes)
                                                 & (data.overview.reduced$meth.genes + data.overview.reduced$meth.probes) > 0
                                                 & (data.overview.reduced$snp.genes + data.overview.reduced$snp.probes) > 0
                                                 ,]
  
  
  return(data.overview.reduced)
}

cpg.sets <- rownames(data.overview.reduced)
save(cpg.sets, file = paste0(GGM.DIR, 'cpgs.RData'))

cpg.data.overview.reduced <- get.cpg.data.overview.reduced(data.overview) 



snp.top.bar <- plot.blank.text('Same Plot for biggest SNP-GGMs', 'A')
snp.bot.bar <- plot.blank.text('Same Plot for smallest SNP-GGMs', 'B')
cpg.top.bar <- plot.overview.bar(data.overview.reduced[1:10, -1], 'C', 'Seed CpG set', F)
cpg.bot.bar <- plot.overview.bar(data.overview.reduced[91:100, -1], 'D', 'Seed CpG set', T)

pdf(file = paste0(PATHS$PLOT.DIR, 'hervS2_cpg_ggm_data_bar.pdf'), width = 6.3, height = 6)
grid.arrange(snp.top.bar, snp.bot.bar, cpg.top.bar, cpg.bot.bar, ncol = 2, widths = c(0.41, 0.59))
dev.off()

colnames(overview.top.flat) <- c('measure', 'cpg.set', 'count')
overview.top.bar <- ggplot(data=overview.top.flat, aes(x=cpg.set, y=count, fill=measure)) + geom_bar(stat="identity")

overview.bot.flat <-  melt(as.matrix(t(data.overview.reduced[134:143, -1])))

ggplot(data.frame()) + geom_point() + xlim(0, 100) + ylim(0, 100) + annotate("text", x=50, y = 50, 'C - Same graph for SNP-GGMs')

colnames(overview.bot.flat) <- c('measure', 'cpg.set', 'count')
ggplot(data=overview.bot.flat, aes(x=cpg.set, y=count, fill=measure)) + geom_bar(stat="identity")


