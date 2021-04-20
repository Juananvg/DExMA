#' Visualization the effect of covariates in data variability
#'
#' It uses the prince and prince.plot function from the swamp package 
#' to visualiza the effect of covariates in data variability
#'
#' @param expressionMatrix A matrix or data frame with genes in rows and samples
#' in columns. An ExpressionSet object can be used too
#' @param pheno A dataframe with samples in rows and covariates in colums. 
#' It should contain only the most important covariates
#'
#' @return A visualization (heatmap) in which it can be seen how the data
#' variability is affected by the covariates. The plot represents the p-values 
#' of each principal component associated with the covariates.
#' 
#' @note Requires the package swamp
#'
#' @author Juan Antonio Villatoro Garcia, 
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#' @seealso \code{\link{batchRemove}}
#' 
#' @references Martin Lauss (2019). swamp: Visualization, Analysis and 
#' Adjustment of High-Dimensional Data in Respect to Sample Annotations. 
#' R package version 1.5.1. \url{https://CRAN.R-project.org/package=swamp}

#' 
#' @examples
#' data(DExMAExampleData)
#' seeCov(listMatrixEX$Study2, listPhenodatas$Study2)
#'
#' @export



seeCov <- function(expressionMatrix, pheno){
    # Remove invariable
    pheno <- pheno[,apply(pheno,2,function(x) length(table(x)))>1]
    # Character variables must be converted in numeric
    pheno <- data.frame(apply(pheno, 2, factor), stringsAsFactors = TRUE)
    print(pheno$race)
    # Number of PC to plot
    if(ncol(expressionMatrix<25)){
        top=ncol(expressionMatrix)
    }else{
        top=25
    }
    res_prince <- prince(expressionMatrix, pheno,top=top)
    prince.plot(res_prince,note=TRUE , notecex=0.5)
}