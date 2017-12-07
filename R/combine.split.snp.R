source('Scripts/R/paths.R')

batch.dir = paste0(PATHS$DATA.DIR, 'SNPs/S2.2kb/')

hervS2.2kb.snp.info.list = lapply(list.files(batch.dir), function(file) {
  return(readRDS(paste0(batch.dir, file))$snpInfo)
})

hervS2.2kb.snp.info <- hervS2.2kb.snp.info.list[[1]]

for (i in c(2:length(hervS2.2kb.snp.info.list))) {
  hervS2.2kb.snp.info <- rbind(hervS2.2kb.snp.info, hervS2.2kb.snp.info.list[[i]])
}

hervS2.2kb.snp.info <- hervS2.2kb.snp.info[order(hervS2.2kb.snp.info$chr, hervS2.2kb.snp.info$pos),]

save(hervS2.2kb.snp.info, file = PATHS$HERVS2.2KB.SNP.INFO.DATA)
