source('Scripts/R/paths.R')

load(PATHS$EXPR.DATA)
load(PATHS$COVARIABLES.DATA)


expr.cov <- merge(t(f4.norm), covars.f4, by.x = 'row.names', by.y = 'ZZ.NR')

message('Data loaded... Calculating residuals.')
source('Scripts/R/residuals.R')

expr.corrected <- get.residuals(expr.cov, 'expr')

save(expr.corrected, file = PATHS$EXPR.CORRECTED.DATA)