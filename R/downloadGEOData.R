#' Download datasets from GEO database
#'
#' Download different ExpressionSets objects from the GEO database

#' @param GEOobject a vector of character where each element represents a
#' GEO object for downloading.
#' 
#' @param directory The directory where the different downloaded GSE Series 
#' Matrix files from GEO will be stored. By default they are downloaded to 
#' the working directory.
#'
#' @return a list of the different ExpressionSets
#'
#' @author Juan Antonio Villatoro Garcia,
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#' @details This function internally uses getGEO function of GEOquery package. 
#' However, downloadGEO allows you to download multiple files at the same time
#' 
#' @references Davis, S. and Meltzer, P. S. GEOquery: a bridge between the 
#' Gene Expression Omnibus (GEO) and BioConductor. 
#' Bioinformatics, 2007, 14, 1846-1847
#' 
#' @examples
#' \dontrun{
#' GEOobjects<- c("GSE4588", "GSE10325")
#' dataGEO<-downloadGEOData(GEOobjects)
#' dataGEO
#' }
#' @export

downloadGEOData <- function(GEOobject, directory=getwd()){
    gset <- list(0)
    i <- 1
    j <- 1
    allnames <- character(0)
    while(i <= length(GEOobject)){
        IDGEO <- GEOobject[i]
        downloadGEO <- getGEO(IDGEO, GSEMatrix=TRUE, 
            destdir=directory, getGPL=FALSE)
        k <- 1
        while(k <= length(downloadGEO)){
            gset[[j]] <- downloadGEO[[k]]
            ##Check if data are log-tranfromed (adapted from GEO2R)
            ex <- exprs(gset[[j]])
            qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5,
                0.75, 0.99, 1.0), na.rm=TRUE))
            LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 &&
                    qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 &&
                            qx[4] > 1 && qx[4] < 2)
            if(LogC){ex[which(ex <= 0)] <- NaN
            exprs(gset[[j]]) <- log2(ex)}
            allnames[j] <- names(downloadGEO[k])
            k <- k+1
            j <- j+1
        }
        i <- i+1
    }
    if(length(allnames) == length(GEOobject)){names(gset) = GEOobject}
    else{names(gset) <- allnames}
    return(gset)
}