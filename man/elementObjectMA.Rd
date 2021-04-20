\name{elementObjectMA}
\alias{elementObjectMA}
\title{Creation of the object to use in meta-analysis}
\usage{
elementObjectMA(
    expressionMatrix,
    pheno = NULL,
    groupPheno,
    expGroup = 1,
    refGroup = 2
)
}
\arguments{
\item{expressionMatrix}{A dataframe or matrix that contanining
genes in rows and samples if columns. An ExpressionSet object 
can be used too.}

\item{pheno}{A data frame or a matrix containing samples in rows and
covariates in columns. If NULL (default), pheno is extracted from the
ExpressionSet object}

\item{groupPheno}{The column name or position from pheno where experimental 
group (cases) and reference group (control) are identified}

\item{expGroup}{The group name or position from groupPheno variable used as
experimental group (cases). By default the first group (character) is taken}

\item{refGroup}{The group name or position from groupPheno variable used as
reference group (control). By default the second group (character) is taken}
}
\value{
An element that can be included in meta-analysis object.
}
\description{
It allows the creation of a element of the object needed to perform
meta-analysis
}
\examples{
data(DExMAExampleData)

ExpressionSetStudy5
newElem <-elementObjectMA(expressionMatrix = ExpressionSetStudy5,
                            groupPheno = "condition",
                            expGroup = c("Diseased", "ill"),
                            refGroup = c("Healthy", "control"))
}
\seealso{
\code{\link{createObjectMA}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
