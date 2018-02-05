#!/usr/bin/env Rscript
#' Preprocess the data  and generate data matrix as input for the network/ggm
#' analysis
#' 
#' @author Johann Hawe
#' 


library(optparse)
library(GenomicRanges)
library(GenomicFeatures)
library(FDb.InfiniumMethylation.hg19)
library(BDgraph)
library(data.table)
library(parallel)
library(illuminaHumanv3.db)
library(rtracklayer)
library(graph)
library(GeneNet)
library(Rgraphviz)
library(RBGL) # for shortest paths
library(Matrix)

source("Scripts/R/lib.R")

################################################################################
# Method definitions
################################################################################

#' Gets the shortest paths between two sets of genes
#'
#' Uses the validated string network and identified the genes on the shortest paths
#' between two given genesets.
#'
#' @param cis List of nodes
#' @param trans List of nodes
#' @param string.db Instance of the string.db to be used (graphNEL)
#' @return A vector of gene symbols being on the shortest path between the two lists
#' of genes as found in the validated string network
#'
#' @author Johann Hawe, Matthias Heinig
#' 
get.string.shortest.paths <- function(cis, trans, snp.genes, string.db) {
  
  g <- string.db
  g.nodes <- nodes(string.db)
  
  # ensure to have only nodes in our giant cluster
  cis <- cis[which(cis %in% g.nodes)]
  trans <- trans[which(trans %in% g.nodes)]
  snp.genes <- snp.genes[which(snp.genes %in% g.nodes)]
  if(length(cis) == 0 | length(trans) == 0 | length(snp.genes) == 0) {
    return(NULL)
  }
  
  # calculate weights for nodes
  prop = propagation(graph2sparseMatrix(g), n.eigs=500, 
                     from=cis, to=trans, sum="both")
  node.weights <- rowSums(prop)
  
  # get the best trans gene
  best.trans = snp.genes[which.max(prop[snp.genes,"from"])]
  cat("best trans genes: ", best.trans, "\n")
  
  ## find the shortest path with maximal weight
  ## we have an algorithm that finds minimum node weight paths so we need
  ## to turn the weighting around
  node.weights = max(node.weights) - node.weights + 1
  
  # extract shortest paths
  sp = min.node.weight.path(g, node.weights, from=cis, to=best.trans)
  nodes <- setdiff(unlist(lapply(sp, "[", "path_detail")), NA)
  nodes <- setdiff(nodes, c(trans, cis))
  return(nodes)
}

#' Gets ranges-data for a specific sentinel SNP
#'
#' Check for each sentinel SNP which genes are near it and its associated CpGs and gets their ranges.
#' Also saves the ranges of the sentinel and its other associated SNPs.
#' Filters by number of CpGs per sentinel.
#'
#' @param id The id (char) of sentinel to be processed.
#' @param cosmo
#' @param trans.meQTL
#' 
#' @value A list containing a relevant ranges for the sentinel snp:
#'  cpg ranges -> cpgs
#'  snp ranges -> snp.ranges
#'  sentinel ranges -> sentinel.range and sentinel.ext.range (extended by window)
#'  gene ranges -> cpg.genes and snp.genes
#' 
#' @author Johann Hawe
#' 
#' @date 02/07/2017
#'
collect.ranges <- function(id, cosmo, trans.meQTL) {
  
  pairs = which(trans.meQTL[,"sentinel.snp"] == id)
  
  # get sentinel idx
  sentinel.idx <- which(cosmo$snp==id)[1]
  
  # the large interval range for the sentinel
  chr <- paste0("chr", trans.meQTL[pairs,"chr.snp"][1])
  start <- trans.meQTL[pairs,"interval.start.snp"][1]
  end <- trans.meQTL[pairs,"interval.end.snp"][1]
  
  sentinel.extrange <- GRanges(chr, IRanges(start,end))
  sentinel.range <- get.snp.range(id)
  names(sentinel.range) <- names(sentinel.extrange) <- id
  
  if(!SENTINEL_ONLY) {
    # get ranges of SNP which the sentinel represents
    sentinel.subset <- trans.meQTL[pairs,]
    snp.ranges <- GRanges(paste0("chr", sentinel.subset[,"chr.snp"]), 
                          IRanges(sentinel.subset[,"pos.snp"],width=1))
    names(snp.ranges) <- sentinel.subset[,"snp"]
    # remove the sentinel from the list of additional SNPs
    snp.ranges <- unique(snp.ranges[!grepl(paste0(id,"$"), names(snp.ranges))])
  }
  
  # load string.db to only filter on genes which are in string
  # the loaded string graph is different to the one saved in 
  # results/current/networks/string_v9_wb.RData...
  load.string.db()
  
  # get cosmo subset
  idxs <- get.trans.cpgs(id, trans.meQTL, cosmo)
  
  # get related genes, i.e. genes near meQTL loci (snp+cpg)
  
  # extend cpgs
  cosmosub <- cosmo[idxs,]
  
  croi <- with(cosmosub, GRanges(paste0("chr", cpg.chr), 
                                 IRanges(cpg.pos,width=2)))
  names(croi) <- as.character(cosmosub[,"cpg"])
  croi <- unique(croi)
  
  # extended sentinel region
  sroi <- sentinel.extrange
  names(sroi) <- id
  
  # for now use the results from the enrichment analysis to get the
  # relevant SNP genes
  genes.sroi <- GENE.ANNOTATION[GENE.ANNOTATION %over% sroi]
  genes.sroi$ids <- probes.from.symbols(genes.sroi$SYMBOL, as.list=T)
  if(is.null(genes.sroi$ids)){
    warning(paste0("No expr probe ids available for snp genes for sentinel ",id))
    return(NULL)
  }
  
  # get genes near our cpg regions
  genes.croi <- get.nearby.ranges(croi, promoters(GENE.ANNOTATION))
  genes.croi$ids <- probes.from.symbols(genes.croi$SYMBOL, as.list=T)
  if(is.null(genes.croi$ids)){
    warning(paste0("No expr probe ids available for cpg genes for sentinel ",id))
    return(NULL)
  }
  
  tfs <- NULL
  sp <- NULL
  
  ## now get also the SYMBOLS of the genes being on the shortest paths
  ## between any snp<->cpg gene combinations
  cat("Collecting shortest paths.\n")
  
  cpgs <- names(croi)
  snp.genes <- unique(genes.sroi$SYMBOL)
  
  # modify string.db to contain our CpGs  
  
  # load the cpg-tf context
  tfbs.ann <- get.chipseq.context(names(croi))
  cpgs.with.tfbs <- cpgs[cpgs %in% rownames(tfbs.ann[rowSums(tfbs.ann)>0,])]
  snp.genes.in.string <- snp.genes[snp.genes %in% nodes(STRING.DB)]
  
  string.db <- STRING.DB
  string.db <- add.to.graphs(list(string.db), id, snp.genes, cpgs.with.tfbs, tfbs.ann)[[1]]
  
  # get tfs connected to cpgs
  tfs = unique(unlist(adj(string.db, cpgs.with.tfbs)))
  cat("Annotated TFs: ", tfs, "\n")
  if(length(tfs)<1){
    warning("No TFs, skipping shortest paths calculation.")
  } else {
    # the nodes we want to keep
    nodeset = c(nodes(STRING.DB), setdiff(tfs, "KAP1"), snp.genes.in.string, cpgs.with.tfbs)
    string.db = subGraph(intersect(nodes(string.db), nodeset), string.db)
    
    syms.sp <- get.string.shortest.paths(cis = cpgs.with.tfbs, 
                                         trans=unique(c(snp.genes.in.string, tfs)), 
                                         snp.genes=snp.genes.in.string,
                                         string.db)
    cat("Shortest paths genes: ", paste0(syms.sp, sep=";"))
    if(length(syms.sp) < 1){
      warning("No shortest path genes.")
    } else {
      sp <- GENE.ANNOTATION[GENE.ANNOTATION$SYMBOL %in% syms.sp]
      sp$ids <- probes.from.symbols(sp$SYMBOL, as.list=T)
    }
    tf.ranges <- GENE.ANNOTATION[GENE.ANNOTATION$SYMBOL %in% tfs]
    tf.ranges$ids <- probes.from.symbols(tf.ranges$SYMBOL, as.list=T)
  }
  
  # construct our result list
  result <-  list(cpgs=croi,sentinel.range=sentinel.range,
                  sentinel.ext.range=sentinel.extrange, 
                  snp.genes=genes.sroi, cpg.genes=genes.croi)
  
  if(!SENTINEL_ONLY) {
    result$snp.ranges <- snp.ranges
  }
  if(!is.null(sp)){
    result$spath <- sp
  }
  if(!is.null(tf.ranges)){
    result$tf.ranges <- tf.ranges
  }
  return(result)
}

#' Adjust the data collected for a specific sentinel SNP
#'
#' Adjusts all data for a given sentinel SNP, based on the previously
#' identified ranges, genes etc done in \code{collect.ranges}. 
#' Needs methylation and expression data available
#' ('meth' and 'expr' variables, respectively). 
#'
#' @param id the id of the sentinel SNP to be processed (needs to be in the list "sentinels")
#' @param ranges The collected ranges for the sentinel SNP
#' @param data A Matrix containing the data retrieved from the ranges object
#'
#' @value NULL
#' 
#' @author Johann Hawe
#'
adjust.data <- function(id, ranges, data, cohort) {
  
  # collect gene symbols to get summarized probe levels
  symbols <- unique(c(ranges$cpg.genes$SYMBOL, 
                      ranges$snp.genes$SYMBOL,
                      ranges$spath$SYMBOL,
                      ranges$tfs$SYMBOL))
  symbols <- symbols[!sapply(symbols, is.na)]
  
  # retrieve the genotype data
  s <- data[,grepl("^rs", colnames(data)),drop=F]
  
  # correct for covariate variation
  g.resid <- get.residuals(data[,!grepl("^rs|^cg", colnames(data))],"expr")
  c.resid <- get.residuals(data[,!grepl("^rs", colnames(data))],"meth")
  
  cat("Loading QTLs.\n")
  eqtls <- load.eqtls(colnames(g.resid))
  meqtls <- load.meqtls(colnames(c.resid))
  
  # get genotype data for eqtls and meqtls
  temp <- meqtls
  if(!is.null(temp)){
    mcols(temp) <- NULL
  }
  temp2 <- eqtls
  if(!is.null(temp2)){
    mcols(temp2) <- NULL
  }
  if(is.null(temp)){
    roi <- temp2 
  } else if(is.null(temp2)){
    roi <- temp
  } else {
    roi<- unique(c(temp,temp2)) 
  }
  
  cat("Number of QTL snps to load:",length(roi), ".\n")
  geno <- get.genotypes(roi, names(roi), cohort)
  geno <- geno[data$geno.ids,,drop=F]
  
  # get rid of cis-eqtl/meqtl effects
  cat("Adjusting for cis eqtls.\n")
  g.resid <- adjust.cis.eqtls(g.resid, eqtls, geno)
  cat("Adjusting for cis meqtls.\n")
  c.resid <- adjust.cis.meqtls(c.resid, meqtls, geno)
  
  cat("Summarizing probe levels.\n")
  # summarize to gene level estimates from expression probes
  g.resid <- summarize(g.resid,symbols = symbols)
  
  # create complete matrix, containing all the information
  data <- cbind.data.frame(g.resid,c.resid,s,stringsAsFactors=F)

  return(data)
}
