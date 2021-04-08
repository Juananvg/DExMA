#' Auxiliary function to check if data are log transfromed and 
#' transformed if it are not log-transformed
#' 
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA should be used.
#'
#' @return The same object with log-transformed  expression matrix 
#'
#' @author Juan Antonio Villatoro Garcia, \email{juan.villatoro@@genyo.es}
#'
#' @seealso \code{\link{createObjectMA}}
#'
#' @examples
#' 
#' data(DExMAExampleData)
#'
#' dataLog(maObject)
#'
#' @export


dataLog <- function(objectMA){
    resultLog <- lapply(objectMA, function(l){
        l <- l[[1]] 
        return(l)})
    resultLog <- lapply(resultLog, function(l){
        l <- .logTransform(l)
        return(l)
    })
    for(i in seq_len(length(objectMA))){
        objectMA[[i]][[1]] <- resultLog[[i]]
    }
    return(objectMA)
}



.logTransform <- function(expressionMatrix){
    if (!any(is.data.frame(expressionMatrix) |
            is.matrix(expressionMatrix) |
            class(expressionMatrix)[1] == "ExpressionSet")){
        stop("ExpressionMatrix must be a dataframe,
            a matrix or a ExpressionSet object")
    }
    if(class(expressionMatrix)[1] == "ExpressionSet"){
        ex <- exprs(expressionMatrix)}
    else{ex <- expressionMatrix}
    ## Check if data are log-tranfromed (adapted from GEO2R)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), 
        na.rm=TRUE))
    LogC <- (qx[5] > 100) ||(qx[6] - qx[1] > 50 &&
            qx[2] > 0) || (qx[2] > 0 &&
                    qx[2] < 1 &&
                    qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)}
    ## Return the matrix changed
    if(class(expressionMatrix)[1] == "ExpressionSet"){
        exprs(expressionMatrix) <- ex}
    else{expressionMatrix <- ex}
    return(expressionMatrix)
}