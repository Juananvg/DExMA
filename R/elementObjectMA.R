#' Creation of the object to use in meta-analysis
#'
#' It allows the creation of a element of the object needed to perform
#' meta-analysis
#'
#' @param expressionMatrix A dataframe or matrix that contanining
#' genes in rows and samples if columns. An ExpressionSet object 
#' can be used too.
#'
#' @param pheno  A data frame or a matrix containing samples in rows and
#' covariates in columns. If NULL (default), pheno is extracted from the
#' ExpressionSet object
#'
#' @param groupPheno The column name or position from pheno where experimental 
#' group (cases) and reference group (control) are identified
#'
#' @param  expGroup The group name or position from groupPheno variable used as
#' experimental group (cases). By default the first group (character) is taken
#'
#' @param refGroup The group name or position from groupPheno variable used as
#' reference group (control). By default the second group (character) is taken
#'
#'
#' @return An element that can be included in meta-analysis object.
#'
#' @author Juan Antonio Villatoro Garcia, \email{juan.villatoro@@genyo.es}
#'
#' @seealso \code{\link{createObjectMA}}
#'
#' @examples
#' data(DExMAExampleData)
#'
#' ExpressionSetStudy5
#' newElem <-elementObjectMA(expressionMatrix = ExpressionSetStudy5,
#'                           groupPheno = "condition",
#'                           expGroup = c("Diseased", "ill"),
#'                           refGroup = c("Healthy", "control"))
#'                           
#' @export
#'                      

elementObjectMA <- function(expressionMatrix, pheno=NULL, groupPheno,
    expGroup=1, refGroup=2){
    if (!any(is.data.frame(expressionMatrix) | is.matrix(expressionMatrix) |
            is(expressionMatrix,"ExpressionSet"))){
        stop("Elements of listEX must be a  dataframe, matrix or
            ExpressionSet object")
    }
    if (!any(is.data.frame(pheno) | is.matrix(pheno) | is.null(pheno))){
        stop("Elements of listPheno must be a dataframe, a matrix or NULL")}
    if (!any(is.character(groupPheno) |is.numeric(groupPheno))){
        stop("groupPheno must be a character or numeric object")}
    if (!any(is.character(expGroup) | is.numeric(expGroup))){
        stop("expGroup must be a character or numeric object")}
    if (!any(is.character(refGroup) | is.numeric(refGroup))){
        stop("refGroup must be a character or numeric object")}
    ## Get data from expression set
    if(is(expressionMatrix, "ExpressionSet")){
        if (is.null(pheno)){
            pheno <- pData(expressionMatrix)}
        expressionMatrix <- exprs(expressionMatrix)
    }
    else{
        if(is.null(pheno)) {
            stop("If DATA is not a ExpressionSet, you must provide
                pheno parameter")}
    }
    colnames(expressionMatrix) = as.character(colnames(expressionMatrix))
    rownames(pheno) = as.character(rownames(pheno))
    pheno = pheno[colnames(expressionMatrix),,drop=FALSE]
    if (is.numeric(groupPheno)){
        groupPheno <- colnames(pheno)[groupPheno]}
    if (is.numeric(expGroup)){
        pheno[,groupPheno] <- factor(pheno[,groupPheno])
        expGroup <- levels(pheno[,groupPheno])[expGroup]}
    if (is.numeric(refGroup)){
        pheno[,groupPheno] <- factor(pheno[,groupPheno])
        refGroup <- levels(pheno[,groupPheno])[refGroup]}
    Inclusion <- rownames(pheno)[pheno[,groupPheno] %in%
            c(expGroup, refGroup)]
    expressionMatrix <- expressionMatrix[,Inclusion]
    pheno <- data.frame(pheno[Inclusion, groupPheno, drop=FALSE])
    colnames(pheno) <- groupPheno
    pheno <- droplevels(pheno)
    ## Convert pheno to 0 and 1 for use the rest of functions
    store<-rep(0, length(pheno))
    for (k in seq_len(length(pheno[,groupPheno]))) {
        if(pheno[,groupPheno][k] %in% expGroup){store[k]=1}
        if(pheno[,groupPheno][k] %in% refGroup){store[k]=0}}
    metalObject <- list(mExpres=expressionMatrix, condition=store)
    return(metalObject)
}