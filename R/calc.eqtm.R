source('Scripts/R/paths.R')

library(MatrixEQTL)

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR
# modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Only associations significant at this level will be saved
pvOutputThreshold.cis = 1e-6

pvOutputThreshold.trans = 1e-8

cisDist = 5e5

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()

## Load gene expression data

expr = SlicedData$new()

expr$fileDelimiter = "\t"
# the TAB character
expr$fileOmitCharacters = "NA"
# denote missing values;
expr$fileSkipRows = 1
# one row of column labels
expr$fileSkipColumns = 1
# one column of row labels
expr$fileSliceSize = 2000
# read file in slices of 2,000 rows
expr$LoadFile(PATHS$F.EQTM.EXPRESSION.FILTERED)


## Load methylation data
meth = SlicedData$new()

meth$fileDelimiter = "\t"
# the TAB character
meth$fileOmitCharacters = "NA"
# denote missing values;
meth$fileSkipRows = 1
# one row of column labels
meth$fileSkipColumns = 1
# one column of row labels
meth$fileSliceSize = 2000
# read file in slices of 2,000 rows
meth$LoadFile(PATHS$F.EQTM.METHYLATION.FILTERED)

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
if (length(PATHS$F.EQTM.COVARIATES.FILTERED) > 0) {
  cvrt$LoadFile(PATHS$F.EQTM.COVARIATES.FILTERED)
}

snpspos = read.table(PATHS$F.EQTM.METHYLATION.POS,
                     header = TRUE,
                     stringsAsFactors = FALSE)
genepos = read.table(PATHS$F.EQTM.EXPRESSION.POS,
                     header = TRUE,
                     stringsAsFactors = FALSE)

## Run the analysis

eqtm.me = Matrix_eQTL_main(
  snps = meth,
  gene = expr,
  cvrt = cvrt,
  output_file_name = PATHS$F.EQTM.TRANS,
  pvOutputThreshold = pvOutputThreshold.trans,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = PATHS$F.EQTM.CIS,
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

save(eqtm.me, file = PATHS$EQTM.ME.DATA)

