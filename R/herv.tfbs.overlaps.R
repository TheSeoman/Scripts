source('Scripts/R/paths.R')

library(rtracklayer)
library(data.table)

if (!file.exists(PATHS$TFBS.DATA)) {
  tfbs = import(paste0(PATHS$TFBS.DIR, "/tfbs/filPeaks_public.bed"))
  ann = t(matrix(unlist(strsplit(
    values(tfbs)[, "name"], ".", fixed = T
  )), nrow = 3))
  colnames(ann) = c("geo_id", "TF", "condition")
  values(tfbs) = DataFrame(name = values(tfbs)[, "name"], data.frame(ann, stringsAsFactors =
                                                                       F))
  
  ## we write out a table with all conditions and select the blood related ones
  conditions = t(matrix(unlist(strsplit(
    unique(values(tfbs)[, "name"]), ".", fixed = T
  )), nrow = 3))
  colnames(conditions) = c("geo_id", "TF", "condition")
  conditions = conditions[order(conditions[, "condition"]), ]
  conditions = conditions[, c(1, 3)]
  conditions = conditions[!duplicated(paste(conditions[, 1], conditions[, 2])), ]
  
  conditions = data.frame(conditions, blood.related = F)
  
  for (term in c(
    "amlpz12_leukemic",
    "aplpz74_leukemia",
    "bcell",
    "bjab",
    "bl41",
    "blood",
    "lcl",
    "erythroid",
    "gm",
    "hbp",
    "k562",
    "kasumi",
    "lymphoblastoid",
    "mm1s",
    "p493",
    "plasma",
    "sem",
    "thp1",
    "u937"
  )) {
    conditions[grep(term, conditions[, 2]), "blood.related"] = TRUE
  }
  
  selected = tfbs[values(tfbs)[, "condition"] %in% conditions[conditions[, "blood.related"], "condition"]]
  
  ## load the encode tfs separately
  encode = as.data.frame(fread(
    paste0(
      PATHS$TFBS.DIR,
      "tfbs/wgEncodeRegTfbsClusteredWithCellsV3.bed"
    ),
    header = F
  ))
  encode = GRanges(
    seqnames = encode[, 1],
    ranges = IRanges(encode[, 2] + 1, encode[, 3]),
    name = paste("ENCODE", encode[, 4], tolower(encode[, 6]), sep = "."),
    geo_id = "ENCODE",
    TF = encode[, 4],
    condition = tolower(encode[, 6])
  )
  
  encode.lcl = encode[grep("gm", values(encode)[, "condition"])]
  values(encode.lcl)[, "condition"] = "lcl"
  encode.k562 = encode[grep("k562", values(encode)[, "condition"])]
  values(encode.k562)[, "condition"] = "k562"
  
  blood.tfbs = c(selected, encode.lcl, encode.k562)
  
  save(blood.tfbs, file = PATHS$TFBS.DATA)
} else {
  load(PATHS$TFBS.DATA)
}

calc.herv.tfbs.overlap <- function(herv.ranges, tfbs.ranges) {
  out <- list()
  out$hits <- findOverlaps(herv.ranges, tfbs.ranges) 
  out$tfbs.ranges <- tfbs.ranges[unique(subjectHits(out$hits))]
  out$herv.ranges <- herv.ranges[unique(queryHits(out$hits))]
  return(out)
}

load(PATHS$HERV.DATA)

for (set in c('S1', 'S2', 'S3')) {
  for (flanking in c('', '.1kb', '.2kb')) {
    cat(paste0('Processing: herv', set, flanking), fill = TRUE)
    overlap.name <- paste0('herv', set, flanking, '.tfbs.overlap')
    assign(overlap.name, calc.herv.tfbs.overlap(get(paste0('herv', set, flanking, '.ranges')), blood.tfbs))
  }
}

