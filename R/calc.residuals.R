source('Scripts/R/paths.R')

require(GenomicRanges)

#' Gets residuals based on linear models from a data matrix (individuals in the rows)
#'
#' Calculates the residual matrix for a given expression or methylation matrix 
#' considering available covariates. Uses linear model for methylation data and linear mixed
#' model with plate as random effect for the expression data. Expects expression probe ids to start
#' with "ILMN_" and methylation probe ids to start with "cg". Supplied covariates must not start with 
#' either of the prefixes
#'
#' @param data the matrix for which to calculate the covariates. Needs to contain covariates themselfes plus
#' either the methylation or expression probes for which to get the residuals
#'
#' @param data.type the type of data in the matrix: either "meth" or "expr". Depending on the type
#' different formulas are used for calculating the linear model.
##
#' @param col.names Optional. The col.names over which to iterate in the dataframe to calculate
#' the residuals on (e.g. probe.ides, gene.names,..)
#' 
#' @return A matrix  where in the colums are the measured entities (e.g.probes) 
#' and in the rows are the samples, containing the calculated residuals. Covariates supplied in the
#' input matrix are discared
#'
get.residuals <- function(data, data.type, col.names=NULL) {
  
  cols <- col.names
  
  if(data.type=="meth") {
    if(is.null(cols)) {
      cols <- colnames(data)[grepl("^cg", colnames(data))];
    }
    res <- lapply(cols, function(n) {
      fm <- as.formula(paste0(n,"~",
                              paste0("1+CD4T+CD8T+NK+Bcell+Mono+",
                                     paste0("PC", 
                                            paste(1:20, "cp", sep="_"), 
                                            collapse="+"))))
      return(lm(fm, data=data, na.action=na.exclude));
    });
  } else if(data.type == "expr") {
    if(is.null(cols)){
      cols <- colnames(data)[grepl("^ILMN_", colnames(data))];
    }
    res <- lapply(cols, function(n) {      
      fm <- as.formula(paste0(n, "~",
                              "1+age+sex+RIN+plate+storage.time"))     
      return(lm(fm, data=data, na.action=na.exclude));
    });
  } else {
    stop("Data type not supported for residual calculation.")
  }
  residual.mat <- matrix(data=unlist(lapply(res, resid)), nrow=nrow(data))
  colnames(residual.mat) <- cols;
  rownames(residual.mat) <- rownames(data);
  
  return(residual.mat)
}

if(!file.exists(PATHS$EXPR.COV.MATRIX)) {
  cat('Loading expression data and covariates...', fill = TRUE)
  load(PATHS$EXPR.DATA)
  load(PATHS$EXPR.RANGES)
  load(PATHS$COVARIABLES.DATA)
  covars.f4$sex <- as.factor(covars.f4$sex);
  covars.f4$plate <- as.factor(covars.f4$plate)
  covars.f4 <- covars.f4[order(covars.f4$ZZ.NR),]
  
  expr.ranges <- expr.ranges[order(expr.ranges)]
  
  expr.data <- t(f4.norm[names(expr.ranges), colnames(f4.norm)])
  expr.matrix <- cbind.data.frame(covars.f4[, 2:6], expr.data, stringsAsFactors = FALSE)
  rownames(expr.matrix) <- rownames(expr.data)
  save(expr.matrix, file = PATHS$EXPR.COV.MATRIX)
} else {
  load(PATHS$EXPR.COV.MATRIX)
}

cat('Calculation expression residuals...', fill = TRUE)                                
expr.residuals <- get.residuals(expr.matrix, 'expr')
cat('Saving expression residuals...', fill = TRUE)                                
save(expr.residuals, file = PATHS$EXPR.RESIDUALS.DATA)

if(!file.exists(PATHS$METH.COV.MATRIX.DATA)) {
  cat('Loading methylation data and covariates...', fill = TRUE)
  houseman <- read.table(PATHS$F.HOUSEMAN ,sep=";", header=T,row.names=1);
  houseman <- houseman[order(rownames(houseman)), ]
  load(PATHS$METH.PCS.DATA)
  pcs <- pcs[order(rownames(pcs)),]
  colnames(pcs) <- paste(colnames(pcs), "cp", sep="_")
  load(PATHS$METH.DATA)
  load(PATHS$METH.RANGES.DATA)
  meth.data <- t(beta[names(meth.ranges), order(colnames(beta))])
  
  meth.matrix <- cbind.data.frame(pcs[rownames(meth.data), 1:20], houseman[rownames(meth.data), 1:5], meth.data ,stringsAsFactors = FALSE)
  rownames(meth.matrix) <- rownames(meth.data)
  save(meth.matrix, file = PATHS$METH.COV.MATRIX.DATA)
} else {
  load(PATHS$METH.COV.MATRIX.DATA)
}
cat('Calculating methylation residuals...', fill = TRUE)

meth.residuals40 <- get.residuals(meth.matrix, 'meth')
cat('Saving methylation residuals...', fill = TRUE)                                
save(meth.residuals, file = PATHS$METH.RESIDUALS.DATA)
