source('Scripts/R/paths.R')

if (!file.exists(PATHS$GSC.DATA)) {
  library(Homo.sapiens)
  
  get.gsc <- function(db, idtype) {
    ## idtype can by ENSEMBL, SYMBOL, ENTREZID, etc. see keytypes(db) for id types
    require(AnnotationDbi)
    require(GSEABase)
    
    frameData = select(db, keys=keys(db, idtype), keytype=idtype, columns=c("GO", "EVIDENCE", idtype))
    frameData = frameData[!is.na(frameData$EVIDENCE),]
    
    
    frame = GOFrame(frameData[,c("GO", "EVIDENCE", idtype)])
    allFrame = GOAllFrame(frame)
    
    gsc <- GeneSetCollection(allFrame, setType=GOCollection())
    return (gsc)
  }
  
  gsc = get.gsc(Homo.sapiens, "SYMBOL")
  dir.create(dirname(gsc.file), recursive=T)
  save(gsc, file=PATHS$GSC.DATA)
} else {
  load(PATHS$GSC.DATA)
}


go.enrichment <- function(genes, universe, gsc, ontologies=c("MF", "BP", "CC")) {
  require(GSEABase)
  require(GOstats)

  go.tab = NULL
  for (ontology in ontologies) {
    params = GSEAGOHyperGParams(name="GO", ontology=ontology, geneIds=genes, universeGeneIds=universe, pvalueCutoff=1, testDirection="over", geneSetCollection=gsc, conditional=FALSE)
    hgt = hyperGTest(params)
    res = data.frame(ontology, summary(hgt, pvalue=1))
    res = data.frame(res, q=p.adjust(res[,"Pvalue"]))
    colnames(res) = gsub("^GO.*ID$", "GOID", colnames(res))
    go.tab = rbind(go.tab, res)
  }
  return(go.tab)
}

extract.significant <- function(herv.set.name, condition, enrichment, cutoff = 0.05) {
  if(enrichment[1, 'q'] > cutoff) {
    significant <- enrichment[1,]
  } else {
    significant <- enrichment[enrichment$q <= cutoff,]
  }
  significant$set <- c(herv.set.name)
  significant$condition <- c(condition)
  significant <- significant[, c(10, 11, 1:9)]
  return(significant)
}


