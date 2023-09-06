#' Set all datasets in the same ID (Official Gene Symbol, Entrez or Ensembl)
#'
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA can be used too.
#' @param finalID A character that indicates the final ID all the different
#' studies will have. To know the available ids, you can write avaliableIDs.
#' @param organism A character that indicates the organism of the data.
#' To know the avaliable organisms write avaliableOrganism
#'
#' @return The same list with all the datasets in the same selected annotation
#'
#' @author Juan Antonio Villatoro Garcia,
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#' @seealso \code{\link{createObjectMA}}
#'
#' @examples
#' data(DExMAExampleData)
#' sameData <- allSameID(objectMA = maObjectDif, 
#' finalID = "GeneSymbol", 
#' organism = "Homo sapiens")
#' sameData
#'
#' @export

allSameID<- function(objectMA, finalID = "GeneSymbol",
    organism="Homo sapiens"){
    listExpMatrix <- lapply(objectMA, function(l) l[[1]])
    if(!is.list(objectMA)) {stop("objectMA must be a list")}
    if(length(organism) != 1) {stop("Only one organism can be included")}
    if(!is.character(finalID))
    {stop("finalID must an object of class character")}
    if(!(finalID %in% DExMAdata::avaliableIDs)){
        stop("finalID not available")}
    initialIDs <- .identifyIDS(objectMA, organism)
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

#Fucntion to identify genes IDs
.identifyIDS <- function(objectMA,organism){
    annot <- DExMAdata::IDsDExMA
    annot <- annot[annot$Organism == organism,]
    genesName <- lapply(objectMA, function(x){return(rownames(x[[1]]))})
    ids <- NA
    for(i in seq_len(length(objectMA))){
        cont <- c(FALSE, FALSE, FALSE)
        res <- genesName[[i]] %in% annot$GeneSymbol
        if (TRUE %in% res){cont[1] <- TRUE}
        res <- genesName[[i]] %in% annot$Entrez
        if (TRUE %in% res){cont[2] <- TRUE}
        res <- genesName[[i]] %in% annot$Ensembl
        if (TRUE %in% res){cont[3] <- TRUE}
        if(sum(cont)<1){stop(
            "One of the datsets has genes IDs that are not avalible")}
        if(sum(cont)>1){stop("One of the datasets has more than a gene ID")}
        if(cont[1]){ids[i] <- "GeneSymbol"}
        if(cont[2]){ids[i] <- "Entrez"}
        if(cont[3]){ids[i] <- "Ensembl"}
    }
    return(ids)
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
.convertID <- function(expressionMatrix, initID, finalID, tableAnnot, 
    organism){
    tableAnnot <- tableAnnot[tableAnnot[,
        initID] %in% rownames(expressionMatrix),]
    tableAnnot <- subset(tableAnnot, tableAnnot$Organism == organism)
    tableAnnot <- tableAnnot[,c(initID, finalID)]
    genes <- rownames(expressionMatrix)
    mEx <- cbind(genes, data.frame(expressionMatrix))
    annotMatrix <- merge(tableAnnot, mEx, by.x = initID, by.y = 1)
    annotMatrix <- annotMatrix[!duplicated(annotMatrix[,c(1,2)]),]
    annotMatrix <-annotMatrix[,-which(colnames(annotMatrix)==initID)]
    agregMatrix <- aggregate(annotMatrix[,-1], by = list(annotMatrix[,finalID]), 
        FUN = median, simplify = TRUE, drop = TRUE)
    rownames(agregMatrix) <- agregMatrix[,1]
    expressionMatrix <- as.matrix(agregMatrix[,-1])
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
            finalID=finalID, tableAnnot=annot, organism=organism)
    }
    #Return values
    return(newExpMatrix)
}