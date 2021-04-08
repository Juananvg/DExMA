#'Calculation of Effects Sizes and their variance
#'
#'This function uses the Hedges'g estimator to calulate the different Effects
#'size and their variances for each genes and for each dataset.
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the diffenrent samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA can be used too.
#'
#' @param missAllow a number that indicates the maximun proportion of missing
#' values allowed in a sample. If the sample has more proportion of missing
#' values the sample will be eliminated. In the other case the missing values
#' will be imputed using the K-NN algorithm.
#'
#'
#' @return A list formed by three elements:
#' First element (ES) is a dataframe were columns are each of the studies
#' (datasets) and rows are the genes. Each element of the dataframe represents
#' the Effect Size.
#' Second element (Var) is a dataframe were columns are each of the studies
#' (datasets) and rows are the genes. Each element of the dataframe represents
#' the variance of the Effect size.
#'
#' @author Juan Antonio Villatoro Garcia, \email{juan.villatoro@@genyo.es}
#'
#' @seealso \code{\link{createObjectMA}}, \code{\link{metaAnalysisDE}}
#'
#' @examples
#'
#' data(DExMAExampleData)
#'
#' resultsEffects <- calculateES(objectMA = maObject, missAllow = 0.3)
#' resultsEffects
#'
#' @export


#Function that combines the others for calculating effects size
calculateES <- function(objectMA, missAllow = 0.3){
    objectMA <- .metaImpute(objectMA, missAllow=missAllow)
    K <- length(objectMA)
    if(is.null(names(objectMA))){
        names(objectMA)<- paste("Study",seq_len(K),sep="")
    }
    Effect <- list(0)
    Variance <- list(0)
    for (i in seq_len(K)) {
        res <- .getES(objectMA[i])
        colnames(res$ES) <- colnames(res$Var) <- names(objectMA[i])
        Effect[[i]] <- res$ES
        Variance[[i]] <- res$Var
    }
    names(Effect) <- names(objectMA)
    names(Variance) <- names(objectMA)
    Total.Effect <- .matrixmerge(Effect)
    Total.Var <- .matrixmerge(Variance)
    Prop.dataset <- as.matrix(1-rowMeans(is.na(Total.Effect)))
    result <- list(ES = Total.Effect, Var = Total.Var)
    return(result)
}





#FUNCTION FOR MERGING MATRIX TAKING INTO ACCOUNT MISSING ROWS
.matrixmerge <- function(lista){
    t.lista <- lapply(lista, t)
    fused <- plyr::rbind.fill.matrix(t.lista)
    fused <- t(fused)
    colnames(fused) <- names(lista)
    return(fused)
}


#Function for delete samples with missing values
.deleteNa <- function(df) {
    no.miss <- colSums(is.na(df[[1]])) <= 0
    df[[1]] = df[[1]][,no.miss]
    df[[2]] = df[[2]][no.miss]
    return(df)
}
#FUNCTION FOR FILTERING SAMPLES WITH MORE THAN % MISSING VALUES
.metaImpute <- function(objectMA,missAllow){
    index.miss <- which(vapply(objectMA,
        FUN = function(y)any(is.na(y[[1]])), 
        FUN.VALUE = TRUE))
    if(length(index.miss)>0){
        for(j in index.miss){
            k<-nrow(objectMA[[j]][[1]])
            rnum<-which(apply(objectMA[[j]][[1]],2,
                function(y) sum(is.na(y))/k)<missAllow)
            print(length(rnum))
            if(length(rnum)>1){
                objectMA[[j]][[1]][,rnum]<-impute.knn(objectMA[[j]][[1]][,rnum],
                    k=10)$data}
            objectMA[[j]]<-.deleteNa(objectMA[[j]])
        }
    }
    return(objectMA)
}


#Functions for calculating effects size
.indCalES<-function(y,l){
    l <- unclass(factor(l))
    n <- table(factor(l))
    ind <- diag(rep(1,length(n)))[l,]
    ym<-y%*%ind%*%diag(1/n)
    ntilde <- 1/sum(1/n)
    m <-sum(n)-2
    cm <- cm<-1-(3/(4*sum(n)-9))
    s <- sqrt((1/(sum(n)-2)*((y^2%*%ind)%*%diag(1/(n-1))-
            ym^2%*%diag(n/(n-1)))%*%(n-1)))
    d <- (ym[,2]-ym[,1])/s
    dprime <- d*cm
    terme1 <- 1/ntilde
    vard <- (sum(l==1)^(-1)+sum(l==2)^(-1))+(d^2)/(2*(sum(l==1)+sum(l==2)))
    vardprime <- sum(1/n)+dprime^2/(2*sum(n))
    result <- cbind( dprime, vardprime)
    colnames(result) <- c( "dprime", "vardprime")
    rownames(result) <- rownames(y)
    return(result)
}


.getES <- function(x){
    K <- length(x)
    ES.m <- Var.m <- N <- n <- NULL
    y <- x[[1]][[1]]
    l <- x[[1]][[2]]
    temp <- .indCalES(y,l)
    ES.m <- as.matrix(temp[,"dprime"])
    Var.m <- as.matrix(temp[,"vardprime"])
    N <- c(N,length(l))
    n <- c(n,table(l))
    rownames(ES.m) <- rownames(y)
    rownames(Var.m) <- rownames(y)
    res <- list(ES = ES.m,Var = Var.m)
    return(res)
}