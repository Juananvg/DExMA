#' Set all datasets in the same ID (Official Gene Symbol, Entrez or Ensembl)
#'
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA can be used too.
#' @param initialIDs A vector in which each element indicates the ID of the 
#' equivalent element of listExpMatrix.
#' To know the avalible ids, you can write avaliableIDs.
#' @param finalID A character that indicates the final ID all the different
#' studies will have. To know the available ids, you can write avaliableIDs.
#' @param organism A character that indicates the organism of the data.
#' To know the avaliable organisms write avaliableOrganism
#'
#' @return The same list with all the datasets in the same selected annotation
#'
#' @author Juan Antonio Villatoro Garcia, \email{juan.villatoro@@genyo.es}
#'
#' @seealso \code{\link{createObjectMA}}
#'
#' @examples
#' data(DExMAExampleData)
#' sameData <- allSameID(objectMA = maObjectDif, initialIDs = c("Entrez",
#' "Entrez", "GeneSymbol", "GeneSymbol"), finalID = "GeneSymbol", 
#' organism = "Homo sapiens")
#' sameData
#'
#' @export

allSameID<- function(objectMA, initialIDs, finalID = "GeneSymbol",
    organism="Homo sapiens"){
    listExpMatrix <- lapply(objectMA, function(l) l[[1]])
    if(length(listExpMatrix) != length(initialIDs)){
        stop("objectMA and initialIDs must have the same length")}
    if(!is.list(objectMA)) {stop("objectMA must be a list")}
    if(length(organism) != 1) {stop("Only one organism can be included")}
    if(!is.character(initialIDs))
    {stop("initialIDs must an object of class character")}
    if(!is.character(finalID))
    {stop("finalID must an object of class character")}
    if(!(finalID %in% DExMAdata::avaliableIDs)){
        stop("finalID not available")}
    for (j in 1:length(initialIDs)){
        if(!(initialIDs[j] %in% DExMAdata::avaliableIDs)){
                stop(initialIDs[j], " not available")}
    }
    k<-length(listExpMatrix)
    newlist <- list(0)
    message("Changing IDs")
    #Progress bar
    pb <- txtProgressBar(min = 0, max = k, style = 3)
    for (i in seq_len(k)) {
        newlist[[i]]<-.toFinalID(listExpMatrix[[i]], initialIDs[[i]], 
            finalID, organism)
        #Progress bar upgrade
        setTxtProgressBar(pb, i)
    }
    close(pb)
    names(newlist) <- names(listExpMatrix)
    for(i in seq_len(length(objectMA))){
        objectMA[[i]][[1]] <- newlist[[i]]
    }
    return(objectMA)
}


#Function to review GeneSymbol
.reviewGeneSymbol <- function(expressionMatrix, tableGeneSymbol,
    Organism = c("Bos taurus", 
        "Caenorhabditis elegans",
        "Canis familiaris", "Danio rerio",
        "Drosophila melanogaster",
        "Gallus gallus",
        "Homo sapiens", "Mus musculus",
        "Rattus norvegicus",
        "Arabidopsis thaliana", 
        "Saccharomyces cerevisiae",
        "Escherichia coli"))
{
    Organism <- match.arg(Organism)
    tableGeneSymbol <- tableGeneSymbol[tableGeneSymbol$Organism == Organism,]
    if(sum(rownames(expressionMatrix) %in% tableGeneSymbol[,"Name"]) >0){
        tableGeneSymbol <- tableGeneSymbol[tableGeneSymbol[,"Name"] %in%
                rownames(expressionMatrix),]
        geneSymbols <- unique(tableGeneSymbol[,"GeneSymbol"])
        annotationMatrix <- matrix(NA,
            ncol=ncol(expressionMatrix),
            nrow=(length(geneSymbols)))
        rownames(annotationMatrix) <- geneSymbols
        colnames(annotationMatrix) <- colnames(expressionMatrix)
        gene <- geneSymbols[1]
        for (gene in geneSymbols){
            id <- as.character(tableGeneSymbol[tableGeneSymbol[,2] == gene,
                'Name'])
            id <- as.character(id)
            if (length(id) > 1){
                unification <- apply(expressionMatrix[id,],2, median)
            } else{
                unification <- expressionMatrix[id,]
            }
            annotationMatrix[gene,] <- unification
        }
        expressionMatrix <- annotationMatrix
    }
    return(expressionMatrix)
}

#Function for change the ids in rownames
.convertID <- function(expressionMatrix, initID, finalID, tableAnnot){
    tableAnnot <- tableAnnot[tableAnnot[,
        initID] %in% rownames(expressionMatrix),]
    finalAnnot <- unique(tableAnnot[,finalID])
    annotationMatrix <- matrix(NA,
        ncol = ncol(expressionMatrix),
        nrow = (length(finalAnnot)))
    rownames(annotationMatrix) <- finalAnnot
    colnames(annotationMatrix) <- colnames(expressionMatrix)
    gene <- finalAnnot[1]
    for (gene in finalAnnot){
        codes <- as.character(
            unique(tableAnnot[tableAnnot[,finalID] == gene,initID]))
        if (length(codes)>1){
            unification<-apply(expressionMatrix[codes,],2, median)
        } else{
            unification<-expressionMatrix[codes,]
        }
        annotationMatrix[gene,] <- unification
    }
    expressionMatrix <- annotationMatrix
    return(expressionMatrix)
}

#Function that includes both function
.toFinalID <- function(expressionMatrix, initialID, finalID,
    organism = c("Bos taurus", "Caenorhabditis elegans",
        "Canis familiaris", "Danio rerio",
        "Drosophila melanogaster", 
        "Gallus gallus",
        "Homo sapiens", "Mus musculus",
        "Rattus norvegicus",
        "Arabidopsis thaliana",
        "Saccharomyces cerevisiae",
        "Escherichia coli")){
    organism <- match.arg(organism)
    #annot <- DExMAdata::IDsDExMA[[id]]
    annot <- DExMAdata::IDsDExMA
    if(initialID == "GeneSymbol"){
        expressionMatrix <- .reviewGeneSymbol(expressionMatrix, 
            tableGeneSymbol=DExMAdata::SynonymsDExMA, organism)
    }
    #End of geneSymbol
    if(initialID==finalID){
        newExpMatrix <- expressionMatrix
    }
    else{
        newExpMatrix <- .convertID(expressionMatrix, initID=initialID, 
            finalID=finalID, tableAnnot=annot)
    }
    #Return values
    return(newExpMatrix)
}
