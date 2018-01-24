
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
