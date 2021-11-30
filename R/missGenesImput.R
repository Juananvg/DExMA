#' Function to use sampleKNN method for missing genes imputation
#' 
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA should be used.
#' 
#' @param k Number of neighbors to be used in the imputation (default=7)
#' 
#' @details
#' 
#' missGenesImput uses k-nearest neighbors in the space of samples
#' to impute the unmeasured genes of the different datasets. 
#' 
#' @return The same object with missing genes imputed
#'
#' @author Juan Antonio Villatoro Garcia,
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#' @seealso \code{\link{createObjectMA}} \code{\link{metaAnalysisDE}}
#' 
#' @references
#' 
#' Christopher A Mancuso, Jacob L Canfield, Deepak Singla, Arjun Krishnan, 
#' A flexible, interpretable, and accurate approach for imputing the expression 
#' of unmeasured genes, Nucleic Acids Research, Volume 48, Issue 21, 
#' 2 December 2020, Page e125, \url{https://doi.org/10.1093/nar/gkaa881}
#' 
#' Alberto Franzin, Francesco Sambo, Barbara di Camillo. bnstruct: an R 
#' package for Bayesian Network structure learning in the presence of 
#' missing data. Bioinformatics, 2017; 33 (8): 1250-1252, 
#' Oxford University Press, \url{https://doi.org/10.1093/bioinformatics/btw807}
#' 
#' @examples
#' 
#' data(DExMAExampleData)
#'
#' missGenesImput(maObject)
#'
#' @export


missGenesImput <- function(objectMA, k = 7){
    fusionMatrix <- .matrixmergeExpression(objectMA)
    #Scale for imputation
    fusionMatrix.mean <- apply(fusionMatrix, 2, calc.mean)
    fusionMatrix.sd <- apply(fusionMatrix, 2, calc.sd)
    fsM.scale <- .scaleS(fusionMatrix, fusionMatrix.mean, fusionMatrix.sd)
    #Imputation
    fsM.scale <- knn.impute(t(fsM.scale), k = k, cat.var = NA)
    fsM.scale <- t(fsM.scale)
    #Get the objectMA back
    fusionMatrix <- .unscaleS(fsM.scale, fusionMatrix.mean, fusionMatrix.sd)
    fusionMatrix <- cbind(fusionMatrix, rep(NA, nrow(fusionMatrix)))
    for(est in seq(length(objectMA))){
        objectMA[[est]][[1]] <- fusionMatrix[,seq(1,ncol(objectMA[[est]][[1]]))]
        fusionMatrix <- fusionMatrix[,-seq(1,ncol(objectMA[[est]][[1]]))]
    }
    return(objectMA)
}


#Function to calculate the mean per sample
calc.mean <- function(x){
    media.x <- mean(x, na.rm = TRUE)
    return(media.x)
}

#Function to calculate the standard deviation per sample
calc.sd <- function(x){
    sd.x <- sd(x, na.rm = TRUE)
    return(sd.x)
}

#Function to scale by sample
.scaleS <- function(x, means, sds){
    for(i in seq(ncol(x))){
        x[,i] <- (x[,i] - means[i]) / sds[i]
    }
    return(x)
}

#Function to unscale by sample
.unscaleS <- function(x, means, sds){
    for(i in seq(ncol(x))){
        x[,i] <- (x[,i] * sds[i]) + means[i]
    }
    return(x)
}


# Merge Expression Matrix
.matrixmergeExpression <- function(lista){
    lista <- lapply(lista, function(x){x<-x[[1]]})
    t.lista <- lapply(lista, t)
    fused <- plyr::rbind.fill.matrix(t.lista)
    sampleNames <- lapply(t.lista,function(x){x <- rownames(x)})
    sampleNames<- plyr::rbind.fill.matrix(sampleNames)
    fused <- t(fused)
    colnames(fused) <- sampleNames
    return(fused)
}

# Merge class vector
.mergeClass <- function(objectMA){
    classes <- lapply(objectMA, function(x){x<-x[[2]]})
    classes <- c(plyr::rbind.fill.matrix(classes))
    return(classes)
} 