source('Scripts/R/paths.R')
source('Scripts/R/lib.R')

require(graph)
require(Rgraphviz)

##cpg seeded ggms

annotate.cpg.ggm <- function(g, d.meta, g.meta) {
  gn <- nodes(g)
  
  nodeDataDefaults(g,'seed.cpg') <- F
  nodeData(g, gn, 'seed.cpg') <- gn %in% d.meta$seed.meth.ids
  
  nodeDataDefaults(g,'extra.cpg') <- F
  nodeData(g, gn, 'extra.cpg') <- gn %in% d.meta$extra.meth.ids
  
  nodeDataDefaults(g,"snp") <- F
  nodeData(g, gn, 'snp') <- grepl("^rs", gn)
  
  nodeDataDefaults(g,"snp.gene") <- F
  nodeData(g, gn, "snp.gene") <- gn %in% c(d.meta$snp.genes, d.meta$snp.no.gene.probes)
  
  nodeDataDefaults(g,"cpg.gene") <- F
  nodeData(g, gn, "cpg.gene") <- gn %in% c(d.meta$seed.meth.genes, d.meta$seed.meth.no.gene.probes, d.meta$extra.meth.genes, d.meta$extra.meth.no.gene.probes)
  
  nodeDataDefaults(g, "tf") <- F
  nodeData(g, gn, "tf") <- gn %in% c(d.meta$seed.tfbs.genes, d.meta$extra.tfbs.genes)
  
  nodeDataDefaults(g, "sp.gene") <- F
  nodeData(g, gn, "sp.gene") <- gn %in% d.meta$path.genes
  
  edgeDataDefaults(g, "iseQTL") <- FALSE
  
  present.eqtls <- g.meta$eqtl.pairs[g.meta$eqtl.pairs$direct, ]
  if(nrow(present.eqtls) > 0) {
    for(i in 1:nrow(present.eqtls)) {
      eqtl.pair <- present.eqtls[i, ]
      expr <- ifelse(is.null(eqtl.pair$gene), eqtl.pair$expr.id, eqtl.pair$gene)
      graph::edgeData(g,eqtl.pair$cpg, expr,"iseQTL") <- T
    }
  }
  
  edgeDataDefaults(g, "iseQTM") <- FALSE 
  present.eqtms <- g.meta$eqtm.pairs[g.meta$eqtm.pairs$direct, ]
  if(nrow(present.eqtms) > 0) {
    for(i in 1:nrow(present.eqtms)) {
      eqtm.pair <- present.eqtms[i, ]
      expr <- ifelse(is.null(eqtm.pair$gene), eqtm.pair$expr.id, eqtm.pair$gene)
      graph::edgeData(g,eqtm.pair$cpg, expr,"iseQTM") <- T
    }
  }
  
  
  return(g)
}

plot.cpg.ggm <- function(g, id, dot.out=NULL){
  library(graph)
  library(Rgraphviz)
  library(GenomicRanges)
  
  # remove any unconnected nodes not within the biggest cc
  ccs <- connComp(g)
  cc.lengths <- unlist(lapply(ccs, length))
  if(length(ccs) > 1) {
    g <- removeNode(nodes(g)[!nodes(g) %in% ccs[[which.max(cc.lengths)]]], g)
  }
  
  n <- nodes(g)
  
  # get trans and cpg gene symbols
  snp.genes <- n[unlist(nodeData(g,n,"snp.gene"))]
  cpg.genes <- n[unlist(nodeData(g,n,"cpg.gene"))]
  tfs <- n[unlist(nodeData(g,n,"tf"))]
  
  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=14, 
                          style="filled", fontname="helvetica"), 
                graph=list(overlap="false", root=id, outputorder="edgesfirst"))
  
  shape = rep("ellipse", numNodes(g))
  names(shape) = n
  shape[unlist(nodeData(g,n,"seed.cpg"))] = "box"
  shape[unlist(nodeData(g,n,"extra.cpg"))] = "box"
  shape[unlist(nodeData(g,n,"snp"))] = "box"
  
  width = rep(0.8, numNodes(g))
  names(width) = n
  width[grep("cg", n)] = 0.4
  
  height = rep(0.3, numNodes(g))
  names(height) = n
  height[grep("cg", n)] = 0.4
  
  label = n
  names(label) = n
  label[grep("cg", n)] = ""
  
  col = rep("#ffffff", numNodes(g))
  names(col) = n
  col[grep("^rs", n)] = "#fab4ad";
  col[unlist(nodeData(g,n,"seed.cpg"))] = "#b5a380";
  col[unlist(nodeData(g,n,"extra.cpg"))] = "#e4d7bc";
  col[cpg.genes] = "#e0e810";
  col[snp.genes] = "#a32c2c"
  if(!is.null(tfs)){
    col[tfs] = "green"
  }
  
  penwidth = rep(1, numNodes(g))
  names(penwidth) = n
  penwidth[snp.genes] = 3
  penwidth[cpg.genes] = 3
  if(!is.null(tfs)){
    penwidth[tfs] = 3  
  }
  
  bordercol = rep("black", numNodes(g));
  names(bordercol) = n;
  #bordercol[cpg.genes] = "#e4d7bc";
  bordercol[unlist(nodeData(g,n,"seed.cpg"))] = "#fab4ad";
  
  nAttrs = list(shape=shape, label=label, width=width, 
                height=height, penwidth=penwidth, fillcolor=col, 
                color=bordercol)
  
  # default color for edges: black
  ecol = rep("black", numEdges(g))
  names(ecol) = edgeNames(g)
  for(i in snp.genes) {
    # color any edge from a SNP to one of our snp genes red
    ecol[grepl("^rs|~rs", names(ecol)) & grepl(i, names(ecol))] = "red"
  }
  
  # set also color for cpgs
  for(cg in cpg.genes){
    # color any edge from a cpg to one of its cpg genes blue (proximity edges)
    ecol[grepl("^cg|~cg", names(ecol)) & grepl(cg, names(ecol))] = "#b3cde2"
  }
  
  # check edgeData and add to colors
  for(edge in names(ecol)){
    n1 <- strsplit(edge,"~")[[1]][1]
    n2 <- strsplit(edge,"~")[[1]][2]
    
    if(unlist(graph::edgeData(g,n1,n2, "iseQTL"))){
      ecol[edge] <- "#decae3"
    }
    if(unlist(graph::edgeData(g,n1,n2, "iseQTM"))){
      ecol[edge] <- "#ccebc5"
    }
  }
  
  dir = rep("none", numEdges(g))
  names(dir) = edgeNames(g)
  
  eAttrs = list(color=ecol, dir=dir)
  
  if(numEdges(g)>500){
    warning("Skipping plotting on device due to large amount of edges")
  } else{
    plot(g, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  if(!is.null(dot.out)){
    # output the dot-information
    toDot(g, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  # return the list of created plotattributes and the possibly modified graph 
  # object
  return(list(graph=g, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs))
}


check.cpg.ggm <- function (set.name, cutoff) {
  GGM.DIR <- '/home/icb/julian.schmidt/data/ggm/hervS2.meqtl.meth.250kb.string/'
  load(paste0(GGM.DIR, 'data.meta.RData'))
  load(paste0(GGM.DIR, 'ggm.meta', cutoff, '.RData'))
  load(paste0(GGM.DIR, 'ggm.overview.', cutoff, '.RData'))
  load(paste0(GGM.DIR, 'ggm/', set.name, '.RData'))
  
  d.meta <- data.meta[[set.name]]
  g.meta <- ggm.meta[[set.name]]
  
  g <- graph.from.fit(ggm, cutoff)
  g <- annotate.cpg.ggm(g, d.meta, g.meta)
  
  DOT.OUT <- paste0(GGM.DIR, 'dot/', set.name, '_', cutoff, '.dot')
  
  out <- plot.cpg.ggm(g,  d.meta$seed.meth.ids[1], DOT.OUT)
  
  return(out)
}

cg06889971.graph <- check.cpg.ggm('cg06889971', 0.9)

##snp seeded ggms
snp <- 'rs11079148'
snp2 <- 'rs114676898'
snp3 <- 'rs1833977'

annotate.snp.ggm <- function(g, meta) {
  gn <- nodes(g)
  
  nodeDataDefaults(g,'cpg') <- F
  nodeData(g, gn, 'cpg') <- gn %in% meta$meth.ids
  
  nodeDataDefaults(g,"snp") <- F
  nodeData(g, gn, 'snp') <- grepl("^rs", gn)
  
  nodeDataDefaults(g,"snp.gene") <- F
  nodeData(g, gn, "snp.gene") <- gn %in% c(meta$snp.genes, meta$snp.no.gene.probes)
  
  nodeDataDefaults(g,"cpg.gene") <- F
  nodeData(g, gn, "cpg.gene") <- gn %in% c(meta$meth.genes, meta$meth.no.gene.probes)
  
  nodeDataDefaults(g, "tf") <- F
  nodeData(g, gn, "tf") <- gn %in% meta$tfbs.genes
  
  nodeDataDefaults(g, "sp.gene") <- F
  nodeData(g, gn, "sp.gene") <- gn %in% meta$path.genes
  
  edgeDataDefaults(g, "isChipSeq") <- FALSE
  edgeDataDefaults(g, "isPPI") <- FALSE 
  
  return(g)
}

plot.snp.ggm <- function(g, id, dot.out=NULL){
  library(graph)
  library(Rgraphviz)
  library(GenomicRanges)
  
  # remove any unconnected nodes not connected to the sentinel snp
  ccs <- connComp(g)
  if(length(ccs) > 1) {
    g <- removeNode(nodes(g)[!nodes(g) %in% ccs[[1]]], g)
  }
  
  
  n <- nodes(g)
  
  # get trans and cpg gene symbols
  snp.genes <- n[unlist(nodeData(g,n,"snp.gene"))]
  cpg.genes <- n[unlist(nodeData(g,n,"cpg.gene"))]
  tfs <- n[unlist(nodeData(g,n,"tf"))]
  
  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=14, 
                          style="filled", fontname="helvetica"), 
                graph=list(overlap="false", root=id, outputorder="edgesfirst"))
  
  
  shape = rep("ellipse", numNodes(g))
  names(shape) = n
  shape[unlist(nodeData(g,n,"cpg"))] = "box"
  shape[unlist(nodeData(g,n,"snp"))] = "box"
  
  width = rep(0.8, numNodes(g))
  names(width) = n
  width[grep("cg", n)] = 0.4
  
  height = rep(0.3, numNodes(g))
  names(height) = n
  height[grep("cg", n)] = 0.4
  
  label = n
  names(label) = n
  label[grep("cg", n)] = ""
  
  col = rep("#ffffff", numNodes(g))
  names(col) = n
  col[grep("^rs", n)] = "#fab4ad";
  col[grep("^cg", n)] = "#e4d7bc";
  col[cpg.genes] = "#e0e810";
  col[snp.genes] = "#a32c2c"
  if(!is.null(tfs)){
    col[tfs] = "green"
  }
  
  penwidth = rep(1, numNodes(g))
  names(penwidth) = n
  penwidth[snp.genes] = 3
  penwidth[cpg.genes] = 3
  if(!is.null(tfs)){
    penwidth[tfs] = 3  
  }
  
  bordercol = rep("black", numNodes(g));
  names(bordercol) = n;
  #bordercol[cpg.genes] = "#e4d7bc";
  bordercol[id] = "#fab4ad";
  
  nAttrs = list(shape=shape, label=label, width=width, 
                height=height, penwidth=penwidth, fillcolor=col, 
                color=bordercol)
  
  # default color for edges: black
  ecol = rep("black", numEdges(g))
  names(ecol) = edgeNames(g)
  for(i in snp.genes) {
    # color any edge from a SNP to one of our snp genes red
    ecol[grepl("^rs|~rs", names(ecol)) & grepl(i, names(ecol))] = "red"
  }
  
  # set also color for cpgs
  for(cg in cpg.genes){
    # color any edge from a cpg to one of its cpg genes blue (proximity edges)
    ecol[grepl("^cg|~cg", names(ecol)) & grepl(cg, names(ecol))] = "#b3cde2"
  }
  
  # check edgeData and add to colors
  for(edge in names(ecol)){
    n1 <- strsplit(edge,"~")[[1]][1]
    n2 <- strsplit(edge,"~")[[1]][2]
    
    if(unlist(graph::edgeData(g,n1,n2, "isPPI"))){
      ecol[edge] <- "#decae3"
    }
    if(unlist(graph::edgeData(g,n1,n2, "isChipSeq"))){
      ecol[edge] <- "#ccebc5"
    }
  }
  
  dir = rep("none", numEdges(g))
  names(dir) = edgeNames(g)
  
  eAttrs = list(color=ecol, dir=dir)
  
  if(numEdges(g)>500){
    warning("Skipping plotting on device due to large amount of edges")
  } else{
    plot(g, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  if(!is.null(dot.out)){
    # output the dot-information
    toDot(g, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  # return the list of created plotattributes and the possibly modified graph 
  # object
  return(list(graph=g, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs))
}


check.snp.ggm <- function (snp, cutoff = 0.9) {
  GGM.DIR <- '/home/icb/julian.schmidt/data/ggm/hervS2.meqtl.snp.250kb.string/'
  load(paste0(GGM.DIR, 'ggm/', snp, '.RData'))
  load(paste0(GGM.DIR, 'data.meta.RData'))
  load(paste0(GGM.DIR, 'ggm.meta', cutoff, '.RData'))
  
  d.meta <- data.meta[[snp]]
  g.meta <- ggm.meta[[snp]]
  
  g <- graph.from.fit(ggm, cutoff)
  g <- annotate.snp.ggm(g, data.meta[[snp]])
  
  DOT.OUT <- paste0(GGM.DIR, 'dot/', snp, '_2_', cutoff, '.dot')
  
  return(plot.snp.ggm(g, snp, DOT.OUT))
}

rs1833977.plot <- check.snp.ggm('rs1833977', 0.9)
rs7805005.plot <- check.snp.ggm('rs7805005', 0.9)
