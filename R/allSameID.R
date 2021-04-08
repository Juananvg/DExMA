#' Set all datasets in the same ID (Official Gene Symbol)
#'
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA can be used too.
#' @param ids A vector in which each element indicates the ID of the equivalent
#' element of listExpMatrix.
#' To know the avalible ids, you can write avaliableIDs
#' @param organism A character that indicates the organism of the data.
#' To know the avaliable organisms write avaliableOrganism
#'
#' @return The same list with all the datasets in Official Gene Symbol
#'
#' @author Juan Antonio Villatoro Garcia, \email{juan.villatoro@@genyo.es}
#'
#' @seealso \code{\link{createObjectMA}}
#'
#' @examples
#' data(DExMAExampleData)
#' sameData <- allSameID(objectMA = maObject, ids = c("Entrez",
#' "Entrez", "Genesymbol", "Genesymbol"), organism = "Homo sapiens")
#' sameData
#'
#' @export



allSameID<- function(objectMA, ids, organism="Homo sapiens"){
    listExpMatrix <- lapply(objectMA, function(l) l[[1]])
    if(length(listExpMatrix) != length(ids)){
        stop("objectMA and ids must have the same length")}
    if(!is.list(objectMA)) {stop("objectMA must be a list")}
    if(length(organism) != 1) {stop("Only one organism can be included")}
    if(!is.character(ids))
    {stop("ids must an object of class character")}
    k<-length(listExpMatrix)
    newlist <- list(0)
    message("Changing IDs")
    #Progress bar
    pb <- txtProgressBar(min = 0, max = k, style = 3)
    for (i in seq_len(k)) {
        newlist[[i]]<-.toGeneSymbol(listExpMatrix[[i]], ids[[i]], organism)
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

#Function to review gene symbol
.reviewGeneSymbol <- function(expressionMatrix, tableGeneSymbol, tableSynonyms,
    Organism = c("Bos taurus", 
        "Caenorhabditis elegans",
        "Canis familiaris", "Danio rerio",
        "Drosophila melanogaster",
        "Gallus gallus",
        "Homo sapiens", "Mus musculus",
        "Rattus norvegicus", "Sus scrofa",
        "Arabidopsis thaliana", 
        "Oryza sativa",
        "Saccharomyces cerevisiae",
        "Aspergillus nidulans",
        "Escherichia coli"))
{
    Organism <- match.arg(Organism)
    tableGeneSymbol <- tableGeneSymbol[tableGeneSymbol$Organism == Organism,]
    tableSynonyms <- tableSynonyms[tableSynonyms$Organism == Organism,]
    if(sum(rownames(expressionMatrix) %in% tableSynonyms[,1]) >0){
        tableGeneSymbol <- tableGeneSymbol[tableGeneSymbol[,1] %in%
                rownames(expressionMatrix),]
        geneSymbols <- unique(tableGeneSymbol[,2])
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


#Function for change to geneSymbol
.convertID <- function(expressionMatrix,tableAnnot){
    tableAnnot <- tableAnnot[tableAnnot[,1] %in% rownames(expressionMatrix),]
    geneSymbols <- unique(tableAnnot[,2])
    annotationMatrix <- matrix(NA,
        ncol = ncol(expressionMatrix),
        nrow = (length(geneSymbols)))
    rownames(annotationMatrix) <- geneSymbols
    colnames(annotationMatrix) <- colnames(expressionMatrix)
    gene <- geneSymbols[1]
    for (gene in geneSymbols){
        probes <- as.character(tableAnnot[tableAnnot[,2] == gene,1])
        if (length(probes)>1){
            unification<-apply(expressionMatrix[probes,],2, median)
        } else{
            unification<-expressionMatrix[probes,]
        }
        annotationMatrix[gene,] <- unification
    }
    expressionMatrix <- annotationMatrix
    return(expressionMatrix)
}



.toGeneSymbol <- function(ExpMatrix, id = "Genesymbol",
    organism = c("Bos taurus", "Caenorhabditis elegans",
        "Canis familiaris", "Danio rerio",
        "Drosophila melanogaster", 
        "Gallus gallus",
        "Homo sapiens", "Mus musculus",
        "Rattus norvegicus", "Sus scrofa",
        "Arabidopsis thaliana", "Oryza sativa",
        "Saccharomyces cerevisiae",
        "Aspergillus nidulans",
        "Escherichia coli")){
    if (!any(is.data.frame(ExpMatrix) | is.matrix(ExpMatrix) |
            class(ExpMatrix)[1] == "ExpressionSet")){
        
        stop("Elements of listExpMatrix must be dataframe,
            matrix or ExpressionSet objects")
    }
    if(class(ExpMatrix)[1] == "ExpressionSet"){
        expressionMatrix <- exprs(ExpMatrix)
    }
    else{
        expressionMatrix <- ExpMatrix
    }
    organism <- match.arg(organism)
    annot <- DExMAdata::IDsDExMA[[id]]
    if(id == "Genesymbol"){
        newExpMatrix <- .reviewGeneSymbol(expressionMatrix, annot,
            DExMAdata::SynonymsDExMA, organism)
    }
    #End of geneSymbol
    else{
        newExpMatrix <- .convertID(expressionMatrix, annot)
    }
    #Return values
    if(class(ExpMatrix)[1] == "ExpressionSet"){
        newData <- ExpressionSet(newExpMatrix)
        pData(newData) <- pData(ExpMatrix)
        annotation(newData) <- "Genesymbol"
    }
    else{newData <- newExpMatrix}
    return(newData)
}
