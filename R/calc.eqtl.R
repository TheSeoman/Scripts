source("Scripts/R/paths.R")

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

library(MatrixEQTL)

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR
# modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = PATHS$F.SNP.FILTERED


# Gene expression file name
expression_file_name = PATHS$F.EXPR.FILTERED

# Covariates file name
# Set to character() for no covariates
covariates_file_name = PATHS$F.COVARIATES.FILTERED


# Output file name
output_file_name.cis = PATHS$F.CIS.EQTL.OUT
output_file_name.trans = PATHS$F.TRANS.EQTL.OUT

# Only associations significant at this level will be saved
pvOutputThreshold.cis = 1e-6

pvOutputThreshold.trans = 1e-8


cisDist = 5e5

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()



## Load genotype data

snps = SlicedData$new()

snps$fileDelimiter = "\t"
# the TAB character
snps$fileOmitCharacters = "NA"
# denote missing values;
snps$fileSkipRows = 1
# one row of column labels
snps$fileSkipColumns = 1
# one column of row labels
snps$fileSliceSize = 2000
# read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)


## Load gene expression data

gene = SlicedData$new()

gene$fileDelimiter = "\t"
# the TAB character
gene$fileOmitCharacters = "NA"
# denote missing values;
gene$fileSkipRows = 1
# one row of column labels
gene$fileSkipColumns = 1
# one column of row labels
gene$fileSliceSize = 2000
# read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)


## Load covariates

cvrt = SlicedData$new()

cvrt$fileDelimiter = "\t"
# the TAB character
cvrt$fileOmitCharacters = "NA"
# denote missing values;
cvrt$fileSkipRows = 1
# one row of column labels
cvrt$fileSkipColumns = 1
# one column of row labels
if (length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
}

snpspos = read.table(PATHS$F.SNP.POS,
                     header = TRUE,
                     stringsAsFactors = FALSE)
genepos = read.table(PATHS$F.EXPR.POS,
                     header = TRUE,
                     stringsAsFactors = FALSE)

## Run the analysis

eqtl.me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name.trans,
  pvOutputThreshold = pvOutputThreshold.trans,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name.cis,
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)


save(eqtl.me, file = PATHS$MAF001.ME.DATA)
