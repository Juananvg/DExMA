#'Calculation p-value for each gene and study
#'
#'This function uses t-test based on limma package in other to obtain the
#'individual p-values for each study and gene
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the diffenrent samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateObjectMA can be used too.
#'
#' @param missAllow a number that indicates the maximun proportion of missing
#' values allowed in a sample. If the sample has more proportion of missing
#' values the sample will be eliminated. In the other case the missing values
#' will be imputed using the K-NN algorithm.
#'
#'
#' @return A list formed by two elements:
#' - First element (p) is a dataframe were columns are each of the studies
#' (datasets) and rows are the genes. Each element of the dataframe represents
#' the p-value.
#' - Second element (FC) is a dataframe were columns are each of the studies
#' (datasets) and rows are the genes. Each element of the dataframe the logFC.
#'
#'
#' @author Juan Antonio Villatoro Garcia,
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#' @seealso \code{\link{createObjectMA}}, \code{\link{metaAnalysisDE}}
#'
#' @examples
#'
#' data(DExMAExampleData)
#'
#' pvalues <- pvalueIndAnalysis(objectMA=maObject, missAllow=0.3)
#' pvalues
#'
#' @export


pvalueIndAnalysis <- function(objectMA, missAllow=0.3){
    objectMA <- .metaImpute(objectMA, missAllow=missAllow)
    K <- length(objectMA)
    if(is.null(names(objectMA))){
        names(objectMA) <- paste("Study",seq_len(K),sep="")
    }
    storeP <- list(0)
    storeF <- list(0)
    #Prepaing data for indiviual analysis
    for(k in seq_len(K)){
        M <- objectMA[[k]][[1]]
        Y <- objectMA[[k]][[2]]
        design <- matrix(data=0, ncol=2, nrow=ncol(M))
        rownames(design) <- colnames(M)
        colnames(design) <- c("CONTROL","CASE")
        for(i in seq_len(length(Y))){
            if(Y[[i]] == 1){
                design[i,1] = 0
                design[i,2] = 1
            }
            else{
                design[i,1] = 1
                design[i,2] = 0
            }
        }
        fit <- lmFit(M, design)
        contrast.matrix <- makeContrasts("CASE-CONTROL", levels = design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        limma<- topTable(fit2, coef=1, adjust.method="BH",number=nrow(M),
            sort.by="none")
        pvalores <- limma$P.Value
        pvalores <- as.matrix(pvalores)
        rownames(pvalores) <- rownames(objectMA[[k]][[1]])
        colnames(pvalores) <- names(objectMA)[k]
        storeP[[k]] <- pvalores
        #selecting Fold-Change
        foldChange <- limma$logFC
        foldChange <- as.matrix(foldChange)
        rownames(foldChange) <- rownames(objectMA[[k]][[1]])
        colnames(foldChange) <- names(objectMA)[k]
        storeF[[k]] <- foldChange
    }
    names(storeF) <- names(objectMA)
    names(storeP) <- names(objectMA)
    fc <- .matrixmerge(storeF)
    p <- .matrixmerge(storeP)
    resultP <- list(p = p, FC = fc)
    return(resultP)
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