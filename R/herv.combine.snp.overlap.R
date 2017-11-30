source('Scripts/R/paths.R')


batch.dir = paste0(PATHS$DATA.DIR, 'SNPs/S2.2kb/')
files = list.files(batch.dir, pattern = 'batch\\d+\\.rds', all.files=FALSE)

snpInfo.list <- list()
snpData.list <- list()

snpInfo.list <- lapply(files, function(file) { return(readRDS(paste0(batch.dir, file))$snpInfo)  })
snpData.list <- lapply(files, function(file) { return(readRDS(paste0(batch.dir, file))$snps)})

unique.snp.ids = unique(unlist(lapply(snpInfo.list, row.names)))

n = length(snpInfo.list)
snpInfo = snpInfo.list[[1]]
for (i in c(2:n)) {
  snpInfo = rbind(snpInfo, snpInfo.list[[i]])
}


snpData = snpData.list[[1]]
for (i in c(2:n)) {
  snpData = rbind(snpData, snpData.list[[i]])
}

hervS2.2kb.snp.info = snpInfo[rownames(snpInfo) %in% unique.snp.ids,]
save(hervS2.2kb.snp.info, file = PATHS$HERV.S2.2KB.SNP.INFO.DATA)
