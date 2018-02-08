source('Scripts/R/paths.R')

batch.dir <- paste0(PATHS$DATA.DIR, 'Methylation/residual.batches/')

batch.files <- list.files(batch.dir)

load(paste0(batch.dir, '1'))
meth.residuals <- meth.residuals.batch
for(i in 2:49) {
  cat(i, fill = T)
  load(paste0(batch.dir, i))
  meth.residuals <- cbind.data.frame(meth.residuals, meth.residuals.batch, stringsAsFactors = F)
}

save(meth.residuals, file = PATHS$METH.RESIDUALS.DATA)
