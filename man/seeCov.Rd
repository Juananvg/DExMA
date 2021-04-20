\name{seeCov}
\alias{seeCov}
\title{Visualization the effect of covariates in data variability}
\usage{
seeCov(expressionMatrix, pheno)
}
\arguments{
\item{expressionMatrix}{A matrix or data frame with genes in rows and samples
in columns. An ExpressionSet object can be used too}

\item{pheno}{A dataframe with samples in rows and covariates in colums. 
It should contain only the most important covariates}
}
\value{
A visualization (heatmap) in which it can be seen how the data
variability is affected by the covariates. The plot represents the p-values 
of each principal component associated with the covariates.
}
\description{
It uses the prince and prince.plot function from the swamp package 
to visualiza the effect of covariates in data variability
}
\note{
Requires the package swamp
}
\examples{
data(DExMAExampleData)
seeCov(listMatrixEX$Study2, listPhenodatas$Study2)

}
\references{
Martin Lauss (2019). swamp: Visualization, Analysis and 
Adjustment of High-Dimensional Data in Respect to Sample Annotations. 
R package version 1.5.1. \url{https://CRAN.R-project.org/package=swamp}
}
\seealso{
\code{\link{batchRemove}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
