#source('Scripts/R/paths.R')

PATHS = c()
PATHS$DATA.DIR = '/media/data/Masterarbeit/data/'

PATHS$F.CIS.EQTL = paste0(PATHS$DATA.DIR, 'eQTL/journal.pone.0093844.s005.XLS')
PATHS$F.TRANS.EQTL = paste0(PATHS$DATA.DIR, 'eQTL/journal.pone.0093844.s004.XLS')

PATHS$S2.SNP.DATA <- paste0(PATHS$DATA.DIR, 'SNPs/hervS2.SNP.RData')
PATHS$S2.SNP.INFO.DATA <- paste0(PATHS$DATA.DIR, 'SNPs/hervS2.snp.info.RData')
PATHS$EXPR.OVERLAP.DATA <- paste0(PATHS$DATA.DIR, 'overlaps/expression.RData')

PATHS$HERV.SMALL.ME <- paste0(PATHS$DATA.DIR, 'eQTL/old_snps/me.RData')
PATHS$HERV.MAF001.ME <- paste0(PATHS$DATA.DIR, 'eQTL/me.RData')

PATHS$F.SNP.FILTERED <- paste0(PATHS$DATA.DIR, 'SNPs/withexpr_sorted')
PATHS$F.EXPR.FILTERED <- paste0(PATHS$DATA.DIR, 'eQTL/expression.tsv')
PATHS$F.COVARIATES.FILTERED <- paste0(PATHS$DATA.DIR, 'eQTL/covariates.tsv')
PATHS$F.CIS.EQTL.OUT <- paste0(PATHS$DATA.DIR, 'eQTL/cis.tsv')
PATHS$F.TRANS.EQTL.OUT <- paste0(PATHS$DATA.DIR, 'eQTL/trans.tsv')
PATHS$F.SNP.POS = paste0(PATHS$DATA.DIR, 'SNPs/snp.pos.tsv')
PATHS$F.EXPR.POS = paste0(PATHS$DATA.DIR, 'eQTL/expression.pos.tsv')

PATHS$F.EXPR <- paste0(PATHS$DATA.DIR, 'Expression/kora_f4_normalized.Rdata')
PATHS$F.COVA <- paste0(PATHS$DATA.DIR, 'Expression/technical_covariables_kora_f4.Rdata')

PATHS$HERV.DATA <- paste0(PATHS$DATA.DIR, 'herv/ranges.RData')
PATHS$HERV.S2.DATA <- paste0(PATHS$DATA.DIR, 'herv/S2.ranges.RData')
PATHS$EXPR.DATA = paste0(PATHS$DATA.DIR, 'F4/Expression/kora_f4_normalized.Rdata')
PATHS$METH.DATA = paste0(PATHS$DATA.DIR, '/F4/KORAF4_illuminamethylation450k_qn_bmiq_n1727/KF4_beta_qn_bmiq.RData')

PATHS$EXPR.OVERLAP.DATA <- paste0(PATHS$DATA.DIR, 'overlaps/expression.RData')
PATHS$METH.OVERLAP.DATA <- paste0(PATHS$DATA.DIR, 'overlaps/methylation.RData')

PATHS$ROADMAP.DIR = paste0(PATHS$DATA.DIR, 'roadmap/') #"/storage/groups/groups_epigenereg/datasets/roadmap/"
PATHS$CHROMHMM.OUT.DIR = paste0(PATHS$DATA.DIR, 'chromHMM/')

PATHS$HERVS1.ANNOT = paste0(PATHS$DATA.DIR, 'herv/hervS1.bed')
PATHS$HERVS2.ANNOT = paste0(PATHS$DATA.DIR, 'herv/hervS2.bed')
PATHS$HERVS3.ANNOT = paste0(PATHS$DATA.DIR, 'herv/hervS3.bed')

PATHS$F.SNP = paste0(PATHS$DATA.DIR, 'SNPs/full_sorted.bgz')
PATHS$F.SAMPLES = paste0(PATHS$DATA.DIR,"F4/individuals.txt")

PATHS$HERVS1.SNP.DATA = paste0(PATHS$DATA.DIR, 'SNPs/hervS1.SNP.RData')
PATHS$HERVS2.SNP.DATA = paste0(PATHS$DATA.DIR, 'SNPs/hervS2.SNP.RData')
PATHS$HERVS3.SNP.DATA = paste0(PATHS$DATA.DIR, 'SNPs/hervS3.SNP.RData')

PATHS$MEQTL.DATA = paste0(PATHS$DATA.DIR, 'meQTL/cosmopairs_combined_151216.RData')
PATHS$MEQTL.HERV.DATA <- paste0(PATHS$DATA.DIR, 'meQTL/meqtl.herv.RData')
PATHS$MEQTL.HERV.1KB.DATA <- paste0(PATHS$DATA.DIR, 'meQTL/meqtl.herv.1kb.RData')
PATHS$MEQTL.HERV.2KB.DATA <- paste0(PATHS$DATA.DIR, 'meQTL/meqtl.herv.2kb.RData')

PATHS$F.COVARIATES = paste0(PATHS$DATA.DIR, 'individuals_all_covariates.csv')
PATHS$F.EXPR.S2.2KB = paste0(PATHS$DATA.DIR, 'eQTL/expression_S2_2kb.tsv')
PATHS$F.SNP.SAMPLES = paste0(PATHS$DATA.DIR, 'SNPs/individuals.txt')
PATHS$F.SNP.SAMPLE.INDICES = paste0(PATHS$DATA.DIR, 'SNPs/snp.sample.indices.txt')
PATHS$F.S2.SNP.INDICES = paste0(PATHS$DATA.DIR, 'SNPs/snp.indices.txt')

PATHS$S2.SNP.DATA = paste0(PATHS$DATA.DIR, 'SNPs/hervS2.SNP.RData')
PATHS$F.SNP.IDS = paste0(PATHS$DATA.DIR, 'SNPs/snp.ids.tsv')
PATHS$F.S2.SNP.IDS = paste0(PATHS$DATA.DIR, 'SNPs/snp.S2.ids.tsv')
