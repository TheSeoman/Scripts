DATA.DIR = '/media/data/Masterarbeit/data/'
COVARIATES = paste0(DATA.DIR, 'individuals_all_covariates.csv')
COVARIATES.OUT = paste0(DATA.DIR, 'eQTL/covariates.tsv')
EXPR.DATA = paste0(DATA.DIR, 'F4/Expression/kora_f4_normalized.Rdata')
EXPR.DATA.OUT = paste0(DATA.DIR, 'eQTL/expression.tsv')

SAMP.SNPS = paste0(DATA.DIR, 'SNPs/individuals.txt')
IND.SNPS.OUT = paste0(DATA.DIR, 'SNPs/snp_indices.txt')
SAMPLES.SNPS.OUT= paste0(DATA.DIR, 'SNPs/snp_samples.txt')

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
write(paste(indices.snp, collapse = ' "\\t" $'), IND.SNPS.OUT)

write(paste(samples.snp, collapse = "\t"), SAMPLES.SNPS.OUT)

expr.ordered <- f4.norm[,as.character(ordered.samples)]
write.table(expr.ordered, EXPR.DATA.OUT, sep = '\t', quote = FALSE)
