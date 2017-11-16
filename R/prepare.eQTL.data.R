DATA.DIR = '/media/data/Masterarbeit/data/'
COVARIATES = paste0(DATA.DIR, 'individuals_all_covariates.csv')
COVARIATES.OUT = paste0(DATA.DIR, 'eQTL/covariates.tsv')
EXPR.DATA = paste0(DATA.DIR, 'Expression/kora_f4_normalized.Rdata')
EXPR.OUT = paste0(DATA.DIR, 'eQTL/expression.tsv')
EXPR.FILT.OUT = paste0(DATA.DIR, 'eQTL/expression_S2_2kb.tsv')
EXPR.POS.OUT = paste0(DATA.DIR, 'eQTL/expression.pos.tsv')

SAMP.SNPS = paste0(DATA.DIR, 'SNPs/individuals.txt')
IND.SAMP.SNP.OUT = paste0(DATA.DIR, 'SNPs/snp.sample.indices.txt')
SAMPLES.SNPS.OUT= paste0(DATA.DIR, 'SNPs/snp.samples.txt')
IND.SNP.ID.OUT = paste0(DATA.DIR, 'SNPs/snp.indices.txt')

EXPR.OVERLAP.DATA <- paste0(DATA.DIR, 'overlaps/expression.RData')
S2.SNP.DATA <- paste0(DATA.DIR, 'SNPs/hervS2.SNP.RData')
SNP.IDS <- paste0(DATA.DIR, 'SNPs/snp.ids.tsv')

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
write.table(expr.pos, EXPR.POS.OUT, sep = '\t', quote = FALSE, row.names = FALSE)

covariates.all <- read.table(COVARIATES, sep = ";", header = TRUE)

id.map <- covariates.all[!is.na(covariates.all$expr_s4f4ogtt) & !is.na(covariates.all$axio_s4f4), c('axio_s4f4', 'expr_s4f4ogtt')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]

covariates.expr <- covariates.all[!is.na(covariates.all$expr_s4f4ogtt), c('ucsex', 'utalteru', 'utbmi', 'ul_wbc', 'expr_s4f4ogtt')]
rownames(covariates.expr) <- covariates.expr$expr_s4f4ogtt
covariates.expr$expr_s4f4ogtt <- NULL

load(EXPR.DATA)
ordered.samples <- sort(intersect(colnames(f4.norm), id.map$expr_s4f4ogtt))

id.map <- id.map[id.map$expr_s4f4ogtt %in% ordered.samples,]

covariates.expr <- covariates.expr[as.character(ordered.samples),]
write.table(t(covariates.expr), COVARIATES.OUT, sep = '\t', quote = FALSE)

samples.snp <- scan(SAMP.SNPS)

sample.index.map <- cbind(samples.snp, c(6:(length(samples.snp)+5)))
rownames(sample.index.map) <- sample.index.map[,1]

indices.snp <- as.vector(sample.index.map[as.character(id.map$axio_s4f4), 2])
write(paste(indices.snp, collapse = ' "\\t" $'), IND.SAMP.SNP.OUT)

write(paste(samples.snp, collapse = "\t"), SAMPLES.SNPS.OUT)

expr.ordered <- f4.norm[,as.character(ordered.samples)]
write.table(expr.ordered, EXPR.OUT, sep = '\t', quote = FALSE)

# filter for probes overlapping with 2kb flanking region for S2
load(EXPR.OVERLAP.DATA)

probes <-  rownames(expr.S2.2kb.overlap$essay.data)
expr.ordered.filtered <- expr.ordered[probes,]
write.table(expr.ordered.filtered, EXPR.FILT.OUT, sep = '\t', quote = FALSE)

# get ids of snps directly in S2
load(S2.SNP.DATA)
S2.snp.ids <- rownames(hervS2.SNPs$snpInfo)
# load list of all snp ids and get indices of snps in S2
all.snp.ids <- scan(SNP.IDS, sep = '\n', what = 'character')
S2.snp.indices <- match(S2.snp.ids, all.snp.ids)
write(paste0('1p;', paste(S2.snp.indices, collapse = 'p;')), IND.SNP.ID.OUT) 

