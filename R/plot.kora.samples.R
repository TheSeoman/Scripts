source('Scripts/R/paths.R')

require(gplots)
library(VennDiagram)


covariates.all <- read.table(PATHS$F.COVARIATES, sep = ";", header = TRUE)

id.map <- covariates.all[covariates.all$expr_s4f4ogtt %in% rownames(expr.residuals) 
                         & covariates.all$axio_s4f4 %in% snp.samples 
                         & covariates.all$meth_f4 %in% rownames(meth.residuals), c('expr_s4f4ogtt', 'axio_s4f4', 'meth_f4')]
id.map <- id.map[order(id.map$expr_s4f4ogtt),]
id.map$expr_s4f4ogtt <- as.character(id.map$expr_s4f4ogtt)
id.map$axio_s4f4 <- as.character(id.map$axio_s4f4)
id.map$meth_f4 <- as.character(id.map$meth_f4)

indices <- list(Expression = which(id.map$expr_s4f4ogtt %in% rownames(expr.residuals)), 
                Genotype = which(id.map$axio_s4f4 %in% snp.samples),
                Methylation = which(id.map$meth_f4 %in% rownames(meth.residuals))) 

venn(indices)

pdf(file = paste0(PATHS$PLOT.DIR, 'samples_venn.pdf'), width = 6, height = 3)
grid.newpage()
venn.plot <- draw.triple.venn(area1 = sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals)),
                              area2 = sum(id.map$axio_s4f4 %in% snp.samples),
                              area3 = sum(id.map$meth_f4 %in% rownames(meth.residuals)),
                              n12 = sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals) & id.map$axio_s4f4 %in% snp.samples),
                              n13 = sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals) & id.map$meth_f4 %in% rownames(meth.residuals)),
                              n23 = sum(id.map$axio_s4f4 %in% snp.samples & id.map$meth_f4 %in% rownames(meth.residuals)),
                              n123 =  sum(id.map$expr_s4f4ogtt %in% rownames(expr.residuals) & id.map$axio_s4f4 %in% snp.samples & id.map$meth_f4 %in% rownames(meth.residuals)),
                              category = c('Expression', 'Genotype', 'Methylation'),
                              fill = c('red', 'blue', 'yellow'),
                              cex = c(1.3),
                              cat.cex = c(1.3),
                              euler.d = T,
                              scaled = T)
dev.off()

grid.draw(venn.plot)
                              