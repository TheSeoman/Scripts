source('Scripts/R/paths.R')

require(gplots)
library(VennDiagram)

load(PATHS$METH.RESIDUALS.DATA)
load(PATHS$EXPR.RESIDUALS.DATA)
load(PATHS$SNP.SAMPLES.DATA)

covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)

id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                         | covariates.all$axio_s4f4 %in% snp.samples 
                         | covariates.all$meth_f4 %in% rownames(meth.residuals), c('expr_s4f4ogtt', 'axio_s4f4', 'meth_f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]
id.map$expr_s4f4ogtt <- as.character(id.map$expr_s4f4ogtt)
id.map$axio_s4f4 <- as.character(id.map$axio_s4f4)
id.map$meth_f4 <- as.character(id.map$meth_f4)

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
                              