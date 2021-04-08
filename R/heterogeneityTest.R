#' It allows checking the heterogeneity of the different studies
#'
#' Shows a qq-plot of the Cochran's test 
#' and the quantiles of I^2 statistic values to measure heterogeneity
#'
#'
#' @param objectMA A list of list. Each list contains two elements. The first
#' element is the expression matrix (genes in rows and sample in columns) and
#' the second element is a vector of zeros and ones that represents the state
#' of the different samples of the expression matrix. 0 represents one group
#' (controls) and 1 represents the other group (cases).
#' The result of the CreateobjectMA can be used too.
#'
#' @details If in the QQ-plot of the Cochran’s test most of the values are 
#' close to the central line (most of the Cochran’s test values
#' are close to the expected distribution, 
#' it can be said that there is homogeneity. In the case that these 
#' values deviate greatly from the expected distribution, it
#' must be assumed that there is heterogeneity.
#' I^2 measures the percentage of variation across studies due to 
#' heterogeneity. To assume homogeneity in the gene expression meta-analysis, 
#' almost all I^2 values (quantiles) must be 0 or at least less than 0.25. 
#'
#' @return Deciles of the I^2 values and a QQ-plot of the Cochran's test
#'
#' @author Juan Antonio Villatoro Garcia, \email{juan.villatoro@@genyo.es}
#' 
#' @references Higgins JPT, Thompson SG. Quantifying heterogeneity in a 
#' meta-analysis. Stat Med 2002;21:1539–58.
#' 
#' Higgins JPT, Thompson SG, Deeks JJ, et al. Measuring inconsistency in 
#' meta-analyses. BMJ 2003;327:557–60.
#'
#' @seealso \code{\link{createObjectMA}}
#'
#' @examples
#'
#' data(DExMAExampleData)
#'
#' heterogeneityTest(maObject)
#'
#' @export

#Heterogeneity test
heterogeneityTest<-function(objectMA){
    allES <- calculateES(objectMA)
    wt <- 1/allES$Var
    temp1 <- wt * allES$ES
    mu.hat <- rowSums(temp1, na.rm = TRUE)/rowSums(wt, na.rm = TRUE)
    Qs <- rowSums(wt * (allES$ES - mu.hat)^2, na.rm = TRUE)
    Qs <- as.data.frame(Qs)
    I2 <- cbind((Qs-length(objectMA)+1)/Qs,0)
    colnames(I2) <- c("Estatistic", "other")
    I2 <- as.data.frame(apply(I2, 1, max))
    colnames(I2) <- "I2"
    quantileI2 <- quantile(I2$I2, probs=seq(from=0, to=1, by=0.05))
    I2 <- as.vector(I2$I2)
    print("I^2 Quantiles")
    colnames(Qs) <- "Qs"
    Qs <- Qs$Qs
    qq.chisq(Qs, df=length(objectMA)-1, slope.one=TRUE,
        main = "QQ plot heterogeneity", col.shade="white")
    return(quantileI2)
}