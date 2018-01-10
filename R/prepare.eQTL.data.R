source('Scripts/R/paths.R')

require(GenomicRanges)

load(PATHS$EXPR.RANGES.DATA)
expr.pos <- data.frame(cbind(names(expr.ranges), as.character(seqnames(expr.ranges)), start(expr.ranges), end(expr.ranges)),
                       row.names = names(expr.ranges))
colnames(expr.pos) <- c('geneid', 'chr', 's1', 's2')

load(PATHS$EXPR.DATA)
samples.snp <- scan(PATHS$F.SNP.SAMPLES)
covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)
#extract ids for used samples
id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(f4.norm) 
                         & covariates.all$axio_s4f4 %in% samples.snp, c('expr_s4f4ogtt', 'axio_s4f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]
expr.pos <- expr.pos[rownames(expr.pos) %in% rownames(f4.norm),]

#filter expression data for probes with known pos and used samples
expr.filtered <- f4.norm[rownames(expr.pos), as.character(id.map$expr_s4f4ogtt)]

#extract covariates for used samples
covariates.filtered <- covariates.all[covariates.all$expr_s4f4ogtt %in% colnames(f4.norm) 
                                  & covariates.all$axio_s4f4 %in% samples.snp, 
                                  c('expr_s4f4ogtt', 'ucsex', 'utalteru', 'utbmi', 'ul_wbc')]
covariates.filtered <- covariates.filtered[order(covariates.filtered$expr_s4f4ogtt),]
colnames(covariates.filtered) <- c('id', 'sex', 'age', 'bmi', 'wbc')

#calculate column indices of used samples in snp-file
sample.index.map <- c(6:(length(samples.snp) + 5))
names(sample.index.map) <- samples.snp

indices.snp <- as.vector(sample.index.map[as.character(id.map$axio_s4f4)])

# generate snp.pos file from snp.ranges, no reordering needed, as snp.ranges are in same order as in snp-file
load(PATHS$SNP.RANGES.DATA)
snp.pos <- data.frame(cbind(names(snp.ranges), as.character(seqnames(snp.ranges)), start(snp.ranges)),
                      row.names = names(snp.ranges))
colnames(snp.pos) <- c('snp', 'chr', 'pos')

write.table(
  expr.pos,
  PATHS$F.EXPR.POS,
  sep = '\t',
  quote = FALSE,
  row.names = FALSE)

write.table(expr.filtered,
            PATHS$F.EXPR.FILTERED,
            sep = '\t',
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)

write.table(t(covariates.filtered),
            PATHS$F.COVARIATES.FILTERED,
            sep = '\t',
            quote = FALSE,
            col.names = FALSE)

write(paste(indices.snp, collapse = ' "\\t" $'),
      PATHS$F.SNP.SAMPLE.INDICES)

