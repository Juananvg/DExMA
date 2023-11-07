#'Performing Meta-analysis
#'
#'It performs meta-analysis using eight different methods.
#'
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA can be used too.
#' 
#' @param effectS A list of three elements. The first element is a dataframe
#' with genes in rows and studies in columns. Each component of the dataframe
#' is the effect of a gene in a study. 
#' The second element of the list is also a dataframe
#' with the same structure, but in this case each component of the dataframe
#' represent the variance of the effect of a gene in a study. This argument
#' should be only used in the case that objectMA argument is null.
#' The third element of the list is also a dataframe
#' with the same structure, but in this case each component of the dataframe
#' represent the log fold change of a gene in a study. 
#' This argument should be only used in the case that objectMA argument is null.
#' 
#' @param pvalues A list of two elements. The first element is a dataframe
#' with genes in rows and studies in columns. Each component of the dataframe
#' is the p-value of a gene in a study. 
#' The second element of the list is also a dataframe
#' with the same structure, but in this case each component of the dataframe
#' represent the log fold change of a gene in a study. This argument
#' should be only used in the case that objectMA argument is null.
#' 
#' @weight A vector of the weights of each dataset. This argument should only 
#' be included in case objectMA is null and you want to 
#' use Stouffer or ACAT method.
#'
#' @param typeMethod A character that indicates the method to be performed.
#' See "Details"for more information
#'
#' @param missAllow a number that indicates the maximum proportion of missing
#' values allowed in a sample. If the sample has more proportion of missing
#' values the sample will be eliminated. In the other case the missing values
#' will be imputed using the K-NN algorithm.
#' 
#' 
#' @param proportionData The minimum proportion of datasets in which a gene
#' must be contained to be included. By default, the gene must be contained 
#' in at least half of the datasets
#'
#' 
#' @details The different meta-analysis methods that can be applied are:
#'\enumerate{
#'    \item \bold{Effects sizes methods:}
#'    \itemize{
#'         \item "FEM": Fixed Effects model
#'         \item "REM": Random Effects model
#'     }
#'     \item \bold{P-value combination mehods}
#'     \itemize{
#'         \item "Fisher": Fisher's methods
#'         \item "Stouffer": Stouffer's method
#'         \item "maxP": maximum p-value method (Wilkinson's method)
#'         \item "minP": minimum p-value method (Tippett's method)
#'         \item "ACAT": Aggregated Cauchy Association Test method
#'     }
#' }
#' 
#'
#' @return A dataframe with the meta-analysis results. Depending on the 
#' applied method, a different dataframe is obtained. For more information 
#' see the package vignette.
#' 
#' @references 
#' 
#' Daniel Toro-Domínguez, Juan Antonio Villatoro-García, 
#' Jordi Martorell-Marugán, Yolanda Román-Montoya, Marta E Alarcón-Riquelme, 
#' Pedro Carmona-Sáez, 
#' A survey of gene expression meta-analysis: methods and applications, 
#' Briefings in Bioinformatics, 2020;, bbaa019, 
#' \url{https://doi.org/10.1093/bib/bbaa019}
#' 
#' Michael Dewey (2020). metap: meta-analysis of significance values. 
#' R package version 1.4
#' 
#' Liu, Y., Chen, S., Li, Z., Morrison, A. C., Boerwinkle, E., & Lin, X. (2019).
#'  ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant 
#'  Analysis in Sequencing Studies. The American Journal of Human Genetics, 
#'  104(3), 410-421. https://doi.org/10.1016/j.ajhg.2019.01.002
#'
#' @author Juan Antonio Villatoro Garcia, 
#' \email{juanantoniovillatorogarcia@@gmail.com}
#'
#'
#' @examples
#'
#' data(DExMAExampleData)
#'
#' ResultsMA <- metaAnalysisDE(objectMA=maObject, typeMethod="REM",
#'                             missAllow=0.3, proportionData=0.5)
#' ResultsMA
#'
#' @export


metaAnalysisDE<-function(objectMA = NULL, effectS = NULL, 
    pvalues = NULL, weight = NULL,
    typeMethod=c("FEM", "REM", "maxP","minP","Fisher","Stouffer", "ACAT"),
    missAllow=0.3, proportionData = 0.5){
    typeMethod <- match.arg(typeMethod)
    if(typeMethod == "FEM" | typeMethod == "REM"){
        if(!is.null(objectMA)){
            effectS <- calculateES(objectMA, missAllow = missAllow)}
        else{
            if(is.null(effectS)){
                stop("You have to add an effectS argument is objectMA is null")}
            names(effectS) <- c("ES", "Var")    
        }
        results <- .metaES(effectS, metaMethod = typeMethod,
            proportionData = proportionData)
    }
    if(typeMethod == "Fisher" | typeMethod == "Stouffer" | 
            typeMethod == "maxP" | typeMethod == "minP" | typeMethod == "ACAT"){
        if(!is.null(objectMA)){
            pvalues <- pvalueIndAnalysis(objectMA, missAllow = missAllow)}
        else{
            if(is.null(pvalues)){
                stop("You have to add a pvalues argument is objectMA is null")}
            names(pvalues) <- c("p", "logFC")
            if(typeMethod == "Stouffer" | typeMethod == "ACAT"){
                if(is.null(weight)){
                    stop("weight argument is required for this method")
                }
                pvalues[[3]] <- weight
                for(i in seq_len(nrow(pvalues[[1]]))[-1]){
                    pvalues[[3]] <- rbind(pvalues[[3]], weight)
                }
                names(pvalues) <- c("p", "logFC", "weights_z")
            }
        }
        results <- .metaPvalue(pvalues, metaMethod = typeMethod,
            proportionData = proportionData)
    }
    return(results)
}

##EFFECTS SIZE FUNCTIONS
.metaES <- function(calESResults, metaMethod=c("FEM","REM"),
    proportionData = 0.5){
    metaMethod <- match.arg(metaMethod)
    K<-ncol(calESResults$ES)
    if(metaMethod == "REM"){
        print("Performing Random Effects Model")
        res <- .getREM(calESResults$ES, calESResults$Var, calESResults$logFC)
        tempFDR <- matrix(res$FDR, ncol=1)
        rownames(tempFDR) <- rownames(calESResults$ES)
        colnames(tempFDR) <- "FDR"
        meta.res <- data.frame(matrix(0, ncol=10,
            nrow = length(rownames(calESResults$ES))))
        rownames(meta.res) <- rownames(calESResults$ES)
        meta.res[,1] <- res$mu.hat
        meta.res[,2] <- res$mu.var
        meta.res[,3] <- res$Qval
        meta.res[,4] <- res$Qpval
        meta.res[,5] <- res$tau2
        meta.res[,6] <- res$zval
        meta.res[,7] <- res$pval
        meta.res[,8] <- tempFDR
        meta.res[,9] <- res$AveFC
        meta.res[,10] <- as.matrix(1-rowMeans(is.na(calESResults$ES)))
        colnames(meta.res) <- c("Com.ES", "ES.var", "Qval", "Qpval", "tau2",
            "Zval", "Pval", "FDR", "AveFC" ,"propDataset")
    }else{
        print("Performing Fixed Effects Model")
        res <- .getFEM(calESResults$ES,calESResults$Var, calESResults$logFC)
        tempFDR <- matrix(res$FDR,ncol=1)
        rownames(tempFDR) <- rownames(calESResults$ES)
        colnames(tempFDR) <- "FDR"
        meta.res <- data.frame(matrix(0, ncol = 7,
            nrow = nrow(calESResults$ES)))
        rownames(meta.res) <- rownames(calESResults$ES)
        meta.res[,1] <- res$mu.hat
        meta.res[,2] <- res$mu.var
        meta.res[,3] <- res$zval
        meta.res[,4] <- res$pval
        meta.res[,5] <- tempFDR[,1]
        meta.res[,6] <- res$AveFC
        meta.res[,7] <- as.matrix(1-rowMeans(is.na(calESResults$ES)))
        colnames(meta.res) <- c("Com.ES", "ES.var", "Zval", "Pval",
            "FDR", "AveFC", "propDataset")
    }
    meta.res<- subset(meta.res, 
        subset = meta.res[,"propDataset"] > 1/
            ncol(calESResults$ES))
    meta.res<- subset(meta.res, 
        subset = meta.res[,"propDataset"] >= proportionData)
    meta.res <- as.data.frame(as.matrix(meta.res))
    attr(meta.res,"metaMethod") <- metaMethod
    return(meta.res)
}

## Fixed Effects Model (FEM)
.getFEM <- function(em,vm,logFC){
    wt <- 1/vm
    mu.hat <- rowSums(wt*em, na.rm=TRUE)/rowSums(wt, na.rm=TRUE)
    mu.var <- 1/rowSums(wt, na.rm=TRUE)
    z.score <- mu.hat/sqrt(mu.var)
    z.p <- 2*(1-pnorm(abs(z.score)))
    qval <- p.adjust(z.p,method = "BH")
    AveFC <- rowSums(wt*logFC, na.rm=TRUE)/rowSums(wt, na.rm=TRUE)
    res <- list(mu.hat = mu.hat,mu.var = mu.var,
        zval = z.score,pval = z.p, FDR = qval, AveFC = AveFC)
    return(res)
}


## RANDOM Effects Model (REM)
## Obtaining the Q that represents the total variance
.getQ <- function(em,vm){
    wt <- 1/vm
    temp1 <- wt * em
    mu.hat <- rowSums(temp1, na.rm=TRUE)/rowSums(wt, na.rm=TRUE)
    Q <- rowSums(wt * (em - mu.hat)^2, na.rm = TRUE)
    return(Q)
}

## Obtaining variance between studies
.getTau2 <- function(Q,vm,k){
    wt <- 1/vm
    s1 <- rowSums(wt, na.rm=TRUE)
    s2 <- rowSums(wt^2, na.rm=TRUE)
    temp <- (Q - (k - 1))/(s1 - (s2/s1))
    tau2 <- pmax(temp,0)
    return(tau2)
}
## Calculating the model
.getREM <- function(em,vm,logFC){
    k <- ncol(em)
    Q.val <- .getQ(em,vm)
    tau2 <- .getTau2(Q.val,vm,k)
    temp.wt <- 1/(vm+tau2)
    mu.hat <- rowSums(temp.wt*em, na.rm=TRUE)/rowSums(temp.wt, na.rm = TRUE)
    mu.var <- 1/rowSums(temp.wt, na.rm=TRUE)
    Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
    z.score <- mu.hat/sqrt(mu.var)
    z.p <- 2*(1-pnorm(abs(z.score)))
    qval <- p.adjust(z.p,method="BH")
    AveFC <- rowSums(temp.wt*logFC, na.rm=TRUE)/rowSums(temp.wt, na.rm = TRUE)
    res <- list(mu.hat = mu.hat, mu.var = mu.var, Qval = Q.val, Qpval = Qpval,
        tau2 = tau2, zval = z.score, pval = z.p, FDR = qval, AveFC = AveFC)
    return(res)
}


## P-VALUES COMBINATION FUNCTIONS
.metaPvalue <-function(resultP, metaMethod=c("maxP", "minP", "Fisher",
    "Stouffer", "ACAT"), proportionData = 0.5){
    p <- resultP$p
    K <- ncol(p)
    nm <- length(metaMethod)
    meta.res <- list(stat=NA,pval=NA,FDR=NA, propDataset=NA)
    meta.res$stat <- meta.res$pval <- meta.res$FDR <- matrix(NA, nrow(p), nm)
    for( i in seq_len(nm)){
        temp <- switch(metaMethod[i], maxP = {.getMaxP(p)},
            minP = {.getMinP(p)},
            Fisher = {.getFisher(p)},
            Stouffer = {.getStouffer(resultP)},
            ACAT = {.getACAT(resultP)})
        meta.res$stat[,i] <- temp$stat
        meta.res$pval[,i] <- temp$pval
        meta.res$FDR[,i] <- temp$FDR
    }
    meta.res$propDataset <- as.matrix(1-rowMeans(is.na(p)))
    colnames(meta.res$stat) <- "Stat"
    colnames(meta.res$pval) <- "Pval"
    colnames(meta.res$FDR) <- "FDR"
    rownames(meta.res$stat) <- rownames(meta.res$pval) <-
        rownames(meta.res$FDR) <- rownames(p)
    attr(meta.res, "nstudy") <- K
    attr(meta.res, "metaMethod") <- metaMethod
    ind.res <- as.data.frame(matrix(rowMeans(resultP$logFC, na.rm=TRUE), ncol=1))
    colnames(ind.res) <- "logFC"
    rownames(ind.res) <- rownames(resultP$logFC)
    res <- data.frame(matrix(0, ncol=5, nrow= nrow(p)))
    colnames(res) <- c("Stat", "Pval", "FDR", "AveFC", "propDataset")
    rownames(res) <-rownames(p)
    res[,1] <- meta.res$stat
    res[,2] <- meta.res$pval
    res[,3] <- meta.res$FDR
    res[,4] <- ind.res
    res[,5] <- meta.res$propDataset
    res<- subset(res, subset = res[,"propDataset"] > 1/ncol(p))
    res<- subset(res, subset = res[,"propDataset"] >= proportionData)
    res <- as.data.frame(as.matrix(res))
    return(res)
}

## Function to obtain the statistic
.sumlogp<-function (p){
    lnp <- log(p)
    chisq <- (-2) * sum(lnp, na.rm = TRUE)
    df <- 2 * length(lnp)
    res <- list(chisq = chisq, df = df, p = pchisq(chisq, df,
        lower.tail = FALSE),
        validp = p)
    return(res)
}

## Fisher's Method
.getFisher <- function(p){
    print("Performing Fisher's method")
    stat <- res <- pval <- fdr <- 0
    todos <- apply(p, 1, .sumlogp)
    ##Extraction of p-values of each one
    for(i in seq_len(nrow(p))){
        stat[i] <- todos[[i]]$chisq
        pval[i] <- todos[[i]]$p
    }
    fdr <- p.adjust(pval, method = "BH")
    res <- list(stat = stat, pval = pval, FDR = fdr)
    names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}

## STOUFFER
## Function to obtain the statistic
.sumzp<-function(pw){
    size <- pw[length(pw)]
    p <- pw[seq_len(size)]
    weights_z <- pw[(size+1):(length(pw)-1)]
    if (length(p) != length(weights_z))
        warning("Length of p and weights differ")
    keep <- (p > 0) & (p < 1) & (is.na(p)==FALSE)
    zp <- (qnorm(p[keep], lower.tail=FALSE) %*%
            weights_z[keep])/sqrt(sum(weights_z[keep]^2))
    res <- list(z = zp, p = pnorm(zp, lower.tail=FALSE),
        weights = weights_z[keep])
    return(res)
}

## Stouffer's method
.getStouffer <- function(resultP){
    p <- resultP$p
    weights_z <- resultP$weights_z
    print("Performing Stouffer's method")
    stat <- res <- pval <- fdr <- 0
    size <- rep(ncol(p), nrow(p))
    pw <- cbind(p, weights_z)
    pw <- cbind(pw,size)
    todos <- apply(pw, 1, .sumzp)
    ## Extraction of p-values of each one
    for(i in seq_len(nrow(p))){
        stat[i] <- todos[[i]]$z
        pval[i] <- todos[[i]]$p
    }
    fdr <- p.adjust(pval, method = "BH")
    res<-list(stat = stat, pval = pval, FDR = fdr)
    names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}

## MINIMUM
## Functions to obtain the statistic
.statisticp<-function (p, r = 1, alpha = 0.05){
    alpha <- ifelse(alpha > 1, alpha/100, alpha)
    stopifnot(alpha > 0, alpha < 1)
    alpha <- ifelse(alpha > 0.5, 1 - alpha, alpha)
    keep <- (is.na(p)==FALSE)
    pi <- p[keep]
    k <- length(pi)
    pi <- sort(pi)
    if(r != 1) {r<- length(pi)}
    pr <- pi[r]
    res <- list(p = pbeta(pr, r, k + 1 - r), pr = pr, r = r,
        critp = qbeta(alpha, r, k + 1 - r), alpha = alpha, validp = pi)
    return(res)
}

## Minimum Tippet's method
.minimumTp<-function (p, alpha = 0.05){
    res <- .statisticp(p, r = 1, alpha)
    return(res)
}

## Minimun P-value method
.getMinP <- function(p){
    print("Performing MinP's method")
    stat <- res <- pval <- fdr <- 0
    pval.all <- apply(p, 1, .minimumTp)
    for(i in seq_len(nrow(p))){
        pval[i] <- pval.all[[i]]$p
        stat[i] <- pval.all[[i]]$pr
    }
    fdr <- p.adjust(pval, method = "BH")
    res <- list(stat = stat, pval = pval, FDR = fdr)
    names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}

## MAXIMUM Wilkinson method
.maximumWp <- function(p, alpha = 0.05){
    k <- length(p)
    res <- .statisticp(p, r = k, alpha)
    res
}

## Maximun P value
.getMaxP<- function(p){
    print("Performing MaxP's method")
    stat <- res <- pval <- fdr <- 0
    pval.all <- apply(p, 1, .maximumWp)
    for(i in seq_len(nrow(p))){
        pval[i] <- pval.all[[i]]$p
        stat[i] <- pval.all[[i]]$pr
    }
    fdr <- p.adjust(pval, method = "BH")
    res <- list(stat = stat, pval = pval, FDR = fdr)
    names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}

## ACAT (Aggregated Cauchy Association Test) method
# Functions to obtain the statistic
.sumCauchy <- function(pw){
    size <- pw[length(pw)]
    p <- pw[seq_len(size)]
    weights_z <- pw[(size+1):(length(pw)-1)]
    if (length(p) != length(weights_z)){
        warning("Length of p and weights differ")}
    keep <- (p > 0) & (p < 1) & (is.na(p)==FALSE)
    weights_c <- weights_z[keep] / sum(weights_z[keep])
    tangent <- tan((0.5 - p[keep])*pi)
    tan_wi <- sum(weights_c * tangent)
    res <- list(q = tan_wi, p = pcauchy(tan_wi, lower.tail = FALSE), 
        weights_c = weights_c)
    return(res)
}

.getACAT <- function(resultP){
    p <- resultP$p
    weights_z <- resultP$weights_z
    print("Performing ACAT's method")
    stat <- res <- pval <- fdr <- 0
    size <- rep(ncol(p), nrow(p))
    pw <- cbind(p, weights_z)
    pw <- cbind(pw,size)
    todos <- apply(pw, 1, .sumCauchy)
    ## Extraction of p-values of each one
    for(i in seq_len(nrow(p))){
        stat[i] <- todos[[i]]$q
        pval[i] <- todos[[i]]$p
    }
    fdr <- p.adjust(pval, method = "BH")
    res<-list(stat = stat, pval = pval, FDR = fdr)
    names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
    return(res)
}
