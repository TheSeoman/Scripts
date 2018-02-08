source('Scripts/R/paths.R')
source('Scripts/R/residuals.R')

cat('Loading methylation + covariates matrix', fill = T)
load(PATHS$METH.COV.MATRIX.DATA)

args <- commandArgs(TRUE)
batch <- as.integer(args[1])
batch.size <- 10000

data.cols <- (26+(batch-1)*batch.size):min(25+batch*batch.size, ncol(meth.matrix))
meth.matrix.batch <- meth.matrix[, c(1:25, data.cols)]
cat(paste0('Calculating residuals for cpgs ', data.cols[1], data.cols[batch.size]), fill = T)
meth.residuals.batch <- get.residuals(meth.matrix.batch, 'meth', colnames(meth.matrix)[data.cols])

save(meth.residuals.batch, file = paste0(PATHS$DATA.DIR, 'Methylation/residual.batches/', batch))
