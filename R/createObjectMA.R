#' Creation of the object to use in meta-analysis
#'
#' It allows the creation of an object to perform meta-analysis.
#'
#' @param listEX A list of dataframes or matrix (genes in rows and sample
#' in columns). A list of ExpressionSets can be used too
#'
#' @param listPheno A list of phenodatas (dataframes or matrix). If the object
#' listEX is a list of ExpressionSets this element can be null.
#'
#' @param namePheno A list or vector of the different colunm names or 
#' positions from the phenodatas where the experimental and reference groups 
#' are identified. Each element of namePheno correspont to its equivalent 
#' element in  the listPheno (default a vector of 1, all the first columns of 
#' each elements of listPheno are selected).
#'
#' @param  expGroups A list of vectors or a vector containing the names or the 
#' positions with which we identify the elements of the experiment groups 
#' (cases) of the namePheno element (default a vector  of 1, all the first 
#' groups are selected)
#'
#' @param refGroups A list of vectors or a vector containing the names or the 
#' positions with which we identify the elements of the reference groups 
#' (control) of the namePheno elements (default a vector  of 1, all the first 
#' groups are selected)
#'
#'
#' @return The object needed to perform meta-analysis. Each list contains 
#' two elements: The first element is the expression matrix (genes in rows 
#' and sample in columns) The second element is a vector of zeros and ones 
#' that represents the state of the diffenrent samples of the expression 
#' matrix. 0 represents reference group (controls) and 1 represents 
#' experimental group (cases).
#'
#' @author Juan Antonio Villatoro Garcia,
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#' @seealso \code{\link{elementObjectMA}}
#'
#' @examples
#' data(DExMAExampleData)
#'
#' phenoGroups = c("condition", "condition", "state", "state")
#' phenoCases = list(Study1 = "Diseased", Study2 = c("Diseased", "ill"),
#'                   Study3 = "Diseased", Study4 = "ill")
#' phenoControls = list(Study1 = "Healthy", Study2 = c("Healthy", "control"),
#'                      Study3 = "Healthy", Study4 = "control")
#'
#' newObjectMA <- createObjectMA(listEX=listMatrixEX,
#'                     listPheno = listPhenodatas, namePheno=phenoGroups, 
#'                     expGroups=phenoCases, refGroups = phenoControls)
#' newObjectMA
#'
#' @export

createObjectMA <- function(listEX, listPheno = NULL,
    namePheno =  c(rep(1, length(listEX))),
    expGroups= c(rep(1, length(listEX))),
    refGroups = c(rep(2, length(listEX)))){
    #check input objects
    if(!is.list(listEX)){
        stop("listEX  must be a list")
    }
    if(length(listEX)<2){stop("There must be at least 2 studies")}
    if (!any(is.list(listPheno) | is.null(listPheno))){
        stop("listPheno must be a List or NULL")
    }
    if(!any(is.list(namePheno) | is.numeric(namePheno) |
            is.character(namePheno))){
        stop("namePheno  must be a list or numeric or character")
    }
    if(!any(is.list(expGroups) | is.numeric(expGroups) |
            is.character(expGroups))){
        stop("expGroups  must be a list or numeric or character")
    }
    if(!any(is.list(refGroups) | is.numeric(refGroups) |
            is.character(refGroups))){
        stop("refGroups  must be a list or numeric or character")
    }
    #Obtaining the object
    container <- list(0)
    for (i in seq_len(length(listEX))) {
        container[[i]] <- elementObjectMA(expressionMatrix=listEX[[i]],
            pheno=listPheno[[i]],
            groupPheno=namePheno[[i]],
            expGroup=expGroups[[i]],
            refGroup=refGroups[[i]])
    }
    #Names of the results
    if(is.null(names(listEX))){
        names(container) <- paste("Study",seq_len(length(listEX)),sep="")}
    else{
        if(length(names(listEX)) >
                length(names(listEX)[is.na(names(listEX))==FALSE])){
            names(container) <- paste("Study",seq_len(length(listEX)),sep="")}
        else{names(container) <- names(listEX)}}
    return(container)
}