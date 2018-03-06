source('Scripts/R/paths.R')

require(gplots)
library(VennDiagram)

load(PATHS$METH.RESIDUALS.DATA)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$SNP.SAMPLES.DATA)

covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)
covariates.all$expr_s4f4ogtt <- as.character(covariates.all$expr_s4f4ogtt)
covariates.all$axio_s4f4 <- as.character(covariates.all$axio_s4f4)
covariates.all$meth_f4 <- as.character(covariates.all$meth_f4)


id.map <- covariates.all[as.character(covariates.all$expr_s4f4ogtt) %in% rownames(expr.residuals) 
                         | as.character(covariates.all$axio_s4f4) %in% snp.samples 
                         | covariates.all$meth_f4 %in% rownames(meth.residuals), ]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]


get.covariate.overview <- function(covariates.map) {
  overview <- data.frame(matrix(ncol = 4, nrow = 3))
  rownames(overview) <- c('Age', 'BMI', 'WBC')
  colnames(overview) <- c('Covariate', 'Minimum', 'Maximum', 'Mean')
  overview[1, ] <- c('Age', min(covariates.map$utalteru), max(covariates.map$utalteru), format(mean(covariates.map$utalteru), digits = 4))
  overview[2, ] <- c('BMI' ,min(na.omit(covariates.map$utbmi)), max(na.omit(covariates.map$utbmi)), format(mean(na.omit(covariates.map$utbmi)), digits = 4))
  overview[3, ] <- c('WBC', min(na.omit(covariates.map$ul_wbc)), max(na.omit(covariates.map$ul_wbc)), format(mean(na.omit(covariates.map$ul_wbc)), digits = 3))
  return(overview)
}

full.kora.overview <- get.covariate.overview(id.map)

select.id.map <- id.map[id.map$expr_s4f4ogtt %in% rownames(expr.residuals) 
                        & id.map$axio_s4f4 %in% snp.samples 
                        & id.map$meth_f4 %in% rownames(meth.residuals),]

select.kora.overview <- get.covariate.overview(select.id.map) 

write.table(full.kora.overview, file = paste0(PATHS$TABLE.DIR, 'full.covariate.overview.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)

write.table(select.kora.overview, file = paste0(PATHS$TABLE.DIR, 'select.covariate.overview.tsv'), 
            quote = F, sep = '\t', row.names = F, col.names = T)

expr.count <-  sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals))
snp.count <- sum(id.map$axio_s4f4 %in% snp.samples)
meth.count <- sum(id.map$meth_f4 %in% rownames(meth.residuals))
pdf(file = paste0(PATHS$PLOT.DIR, 'samples_venn.pdf'), width = 4.5, height = 2)
grid.newpage()
venn.plot <- draw.triple.venn(area1 = expr.count,
                              area2 = snp.count,
                              area3 = meth.count,
                              n12 = sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals) & id.map$axio_s4f4 %in% snp.samples),
                              n13 = sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals) & id.map$meth_f4 %in% rownames(meth.residuals)),
                              n23 = sum(id.map$axio_s4f4 %in% snp.samples & id.map$meth_f4 %in% rownames(meth.residuals)),
                              n123 =  sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals) & id.map$axio_s4f4 %in% snp.samples & id.map$meth_f4 %in% rownames(meth.residuals)),
                              category = c(paste0('Expression (', expr.count, ')'),
                                           paste0('Genotype (', snp.count, ')'), 
                                           paste0('Methylation (', meth.count, ')')),
                              fill = c('red', 'blue', 'yellow'),
                              cex = c(1),
                              cat.cex = c(1),
                              euler.d = T,
                              scaled = T,
                              cat.dist = c(0.07, 0.07, 0.04),
                              cat.pos = c(330, 30, 180))
dev.off()

grid.draw(venn.plot)
                              