source('Scripts/R/paths.R')
source('Scripts/R/residuals.R')

cat('Loading methylation + covariates matrix', fill = T)
load(PATHS$METH.COV.MATRIX.DATA)

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 1000

meth.matrix.batch <- meth.matrix[, c(1:25, (26+(batch-1)*batch.size):min(25+batch*batch.size, ncol(meth.matrix)))]
cat(paste0('Calculating residuals for cpgs ', (26+(batch-1)*batch.size), ' to ', min(26+batch*batch.size, ncol(meth.matrix))), fill = T)
meth.residuals.batch <- get.residuals(meth.matrix.batch, 'meth')

save(meth.residuals.batch, file = paste0(PATHS$DATA.DIR, 'Methylation/residual.batches/', batch))
