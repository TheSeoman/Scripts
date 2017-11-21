source('Scripts/R/paths.R')

load(PATHS$F.EXPR)
load(PATHS$F.COVA)

covData <- merge(t(f4.norm), covars.f4, by.x = 'row.names', by.y = 'ZZ.NR')

source('./residuals.R')

resMatrix <- get.Residuals(covData, 'expr')

save(resMatrix, file = 'resMatrix.rdata')