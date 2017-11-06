KORA.DIR <- "/storage/groups/groups_epigenereg/analyses/PV_K14115g_Heinig/"
F.EXPR <- paste0(KORA.DIR, "data/20160204/Expression/kora_f4_normalized.Rdata");
F.COVA <- paste0(KORA.DIR, "data/20160204/Expression/technical_covariables_kora_f4.Rdata");

load(F.EXPR)
load(F.COVA)

covData <- merge(t(f4.norm), covars.f4, by.x = 'row.names', by.y = 'ZZ.NR')

source('./residuals.R')

resMatrix <- get.Residuals(covData, 'expr')

save(resMatrix, file = 'resMatrix.rdata')