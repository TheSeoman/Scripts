source('Scripts/R/paths.R')

library(MatrixEQTL)

use.residuals = T
## Settings

useModel = modelLINEAR
pvOutputThreshold.cis = 1e-6
pvOutputThreshold.trans = 1e-8
cisDist = 5e5
errorCovariance = numeric()

expr = SlicedData$new()
expr$fileDelimiter = "\t"
expr$fileOmitCharacters = "NA"
expr$fileSkipRows = 1
expr$fileSkipColumns = 1
expr$fileSliceSize = 2000
expr$LoadFile(PATHS$F.EQTM.EXPRESSION.FILTERED)

meth = SlicedData$new()
meth$fileDelimiter = "\t"
meth$fileOmitCharacters = "NA"
meth$fileSkipRows = 1
meth$fileSkipColumns = 1
meth$fileSliceSize = 2000
meth$LoadFile(PATHS$F.EQTM.METHYLATION.FILTERED)

## Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
if (!use.residuals) {
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

