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
#' @return A list formed by two elements:
#' \itemize{
#' \item{First element (objectMA) the same objectMA with missign genes imputed}
#' \item {Second element (imputIndicators) a list with 4 different objects:}
#'     \itemize{
#'         \item{imputValuesSample: Number of missing values imputed 
#'         per sample}
#'         \item{imputPercentageSample: Percentage of missing values 
#'         imputed per gene}
#'         \item{imputValuesGene: Number of missing values imputed per sample}
#'         \item{imputPercentageGene: Percentage of missing values imputed 
#'         per gene}
#'     }
#'}
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
    preObjectMA <- objectMA
    fusionMatrix <- .matrixmergeExpression(objectMA)
    preFus <- fusionMatrix
    #Scale for imputation
    fusionMatrix.mean <- apply(fusionMatrix, 2, calc.mean)
    fusionMatrix.sd <- apply(fusionMatrix, 2, calc.sd)
    fsM.scale <- .scaleS(fusionMatrix, fusionMatrix.mean, fusionMatrix.sd)
    #Imputation
    fsM.scale <- knn.impute(t(fsM.scale), k = k, cat.var = NA)
    fsM.scale <- t(fsM.scale)
    #Get the objectMA back
    fusionMatrix <- .unscaleS(fsM.scale, fusionMatrix.mean, fusionMatrix.sd)
    posFus <- fusionMatrix
    fusionMatrix <- cbind(fusionMatrix, rep(NA, nrow(fusionMatrix)))
    for(est in seq(length(objectMA))){
        objectMA[[est]][[1]] <- fusionMatrix[,seq(1,ncol(objectMA[[est]][[1]]))]
        fusionMatrix <- fusionMatrix[,-seq(1,ncol(objectMA[[est]][[1]]))]
    }
    #Checking imputed genes
    #Number of values imputed
    print(paste0("Number of values imputed ", sum(is.na(preFus)), " (",
        round(sum(is.na(preFus)) / ((dim(posFus)[1] * dim(posFus)[2])) ,3), 
        " %)"))
    #Number of imputed genes in each datasets
    percenta <- 0
    for(i in seq_len(length(objectMA))){
        percen <- (nrow(objectMA[[i]][[1]]) - nrow(preObjectMA[[i]][[1]])) / 
            nrow(objectMA[[i]][[1]]) * 100
        percenta[i] <- percen
        numGen <- nrow(objectMA[[i]][[1]]) - nrow(preObjectMA[[i]][[1]])
        print(paste0("Number of genes imputed in ", names(preObjectMA)[i], ": ", 
            numGen, " of ", nrow(objectMA[[i]][[1]]),  
            " (", round(percen,1), "%)"))}
    #Number and percentage values imputed per sample:
    imputSample <- apply(preFus, 2, function(x){sum(is.na(x))})
    imputPerSample <- apply(preFus, 2, function(x){
        sum(is.na(x)) / length(x) * 100})
    #Number and percentage values imputed per gen:
    imputGene <- apply(preFus, 1, function(x){sum(is.na(x))})
    imputPerGene <- apply(preFus, 1, function(x){
        sum(is.na(x)) / length(x) * 100})
    #Final object
    imputIndicators <- list(imputValuesSample = imputSample,
        imputPercentageSample = imputPerSample, imputValuesGene = imputGene, 
        imputPercentageGene = imputPerGene)
    ImputResults <- list(objectMA = objectMA, imputIndicators = imputIndicators)
    return(ImputResults)
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