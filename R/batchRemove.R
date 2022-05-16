#' Elimination of covariates batch effect or bias
#'
#' It eliminates the effects of batch or bias of the covariates
#'
#' @param expressionMatrix A matrix or data frame with genes in rows and samples
#' in columns. An ExpressionSet object can be used too
#'
#' @param pheno A dataframe with samples in rows and covariates in colums.
#'
#' @param formula Formula of the covariates that are wanted to be corrected
#'
#' @param mainCov Name of the main covariate to be corrected
#'
#' @param nameGroup Name of the column of the Phenodata object in which the 
#' reference groups (cases and controls) are
#'
#' @param ... Other arguments are passed to lmFit of limma package
#'
#' @return The Expression Matrix with the bias or batch effect corrected.
#' Moreover a plot of the visualization of the association between principal 
#' components and covariates is shown.
#'
#'
#' @author Juan Antonio Villatoro Garcia,
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#'
#' @references Martin Lauss (2019). swamp: Visualization, Analysis and 
#' Adjustment of High-Dimensional Data in Respect to Sample Annotations. 
#' R package version 1.5.1. \url{https://CRAN.R-project.org/package=swamp}
#'
#' @examples
#' data(DExMAExampleData)
#' batchRemove(listMatrixEX$Study2, listPhenodatas$Study2, formula=~gender+race,
#' nameGroup="condition")
#'
#' @export
#'


batchRemove <- function(expressionMatrix, pheno, formula, mainCov=NULL,
    nameGroup, ...){
    # Remove invariable
    pheno<-pheno[,apply(pheno,2,function(x) length(table(x)))>1]
    Diagnosis <- pheno[,nameGroup]
    design.disease <- model.matrix(~Diagnosis, data=pheno)
    design.covariates <- model.matrix(formula, data=pheno)
    
    if(is.null(mainCov)){
        ebatch <- removeBatchEffect(expressionMatrix,
            covariates=design.covariates,
            design = design.disease, ... = ...)
    }else{
        ebatch <- removeBatchEffect(expressionMatrix, batch=pheno[,mainCov],
            covariates=design.covariates,
            design=design.disease, ... = ...)
    }
    # Character variables must be converted in numeric
    pheno.factors <- pheno
    for(i in seq_len(ncol(pheno.factors))){
        pheno.factors[,i] <- factor(pheno.factors[,i])
    }
    # Number of PC to plot
    if(ncol(expressionMatrix<25)){
        top=ncol(expressionMatrix)
    }else{
        top=25
    }
    # Visualization
    prince.plot(prince(ebatch,pheno.factors,top=top),note = TRUE ,notecex = 0.5)
    return(ebatch)
}