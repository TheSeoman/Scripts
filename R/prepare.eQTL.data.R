source('Scripts/R/paths.R')

get.expression.positions <- function () {
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
  table <- cbind(ids, chrs, start, end)
  colnames(table) <- c('geneid', 'chr', 's1', 's2')
  return(table)
}

expr.pos <- get.expression.positions();
write.table(expr.pos, PATHS$F.EXPR.POS, sep = '\t', quote = FALSE, row.names = FALSE)

covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)

id.map <- covariates.all[!is.na(covariates.all$expr_s4f4ogtt) & !is.na(covariates.all$axio_s4f4), c('axio_s4f4', 'expr_s4f4ogtt')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]

covariates.expr <- covariates.all[!is.na(covariates.all$expr_s4f4ogtt), c('ucsex', 'utalteru', 'utbmi', 'ul_wbc', 'expr_s4f4ogtt')]
rownames(covariates.expr) <- covariates.expr$expr_s4f4ogtt
covariates.expr$expr_s4f4ogtt <- NULL

load(PATHS$EXPR.DATA)
load(PATHS$EXPR.CORRECTED.DATA)
expr.corrected <- t(resMatrix)
colnames(expr.corrected) <- colnames(f4.norm)
ordered.samples <- sort(intersect(colnames(expr.corrected), id.map$expr_s4f4ogtt))

id.map <- id.map[id.map$expr_s4f4ogtt %in% ordered.samples,]

covariates.expr <- covariates.expr[as.character(ordered.samples),]
write.table(t(covariates.expr), F.COVARIATES.FILTERED, sep = '\t', quote = FALSE)

samples.snp <- scan(PATHS$F.SNP.SAMPLES)

sample.index.map <- cbind(samples.snp, c(6:(length(samples.snp)+5)))
rownames(sample.index.map) <- sample.index.map[,1]

indices.snp <- as.vector(sample.index.map[as.character(id.map$axio_s4f4), 2])
write(paste(indices.snp, collapse = ' "\\t" $'), PATHS$F.SNP.SAMPLE.INDICES)

write(paste(samples.snp, collapse = "\t"), PATHS$F.SNP.SAMPLES)

expr.ordered <- expr.corrected[,as.character(ordered.samples)]
write.table(expr.ordered, PATHS$F.EXPR.FILTERED, sep = '\t', quote = FALSE)

# filter for probes overlapping with 2kb flanking region for S2
#load(PATHS$EXPR.OVERLAP.DATA)

#probes <-  rownames(expr.S2.2kb.overlap$essay.data)
# expr.ordered.filtered <- expr.ordered[probes,]
# write.table(expr.ordered.filtered, PATHS$F.EXPR.S2.2KB, sep = '\t', quote = FALSE)

# get ids of snps directly in S2
# load(PATHS$S2.SNP.DATA)
# S2.snp.ids <- rownames(hervS2.SNPs$snpInfo)
# load list of all snp ids and get indices of snps in S2
# all.snp.ids <- scan(PATHS$F.SNP.IDS, sep = '\n', what = 'character')
# S2.snp.indices <- match(S2.snp.ids, all.snp.ids)
# write(paste0('1p;', paste(S2.snp.indices, collapse = 'p;')), PATHS$F.S2.SNP.IDS) 

