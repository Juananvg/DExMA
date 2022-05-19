\name{batchRemove}
\alias{batchRemove}
\title{Elimination of covariates batch effect or bias}
\usage{
batchRemove(expressionMatrix, pheno, formula, mainCov = NULL, nameGroup, ...)
}
\arguments{
\item{expressionMatrix}{A matrix or data frame with genes in rows and samples
in columns. An ExpressionSet object can be used too}

\item{pheno}{A dataframe with samples in rows and covariates in colums.}

\item{formula}{Formula of the covariates that are wanted to be corrected}

\item{mainCov}{Name of the main covariate to be corrected}

\item{nameGroup}{Name of the column of the Phenodata object in which the 
reference groups (cases and controls) are}

\item{...}{other arguments are passed to lmFit fucntion of limma package}

}
\value{
The Expression Matrix with the bias or batch effect corrected.
Moreover a plot of the visualization of the association between principal 
components and covariates is shown.
}
\description{
It eliminates the effects of batch or bias of the covariates
}

\examples{
data(DExMAExampleData)
batchRemove(listMatrixEX$Study2, listPhenodatas$Study2, formula=~gender+race,
nameGroup="condition")

}
\references{
Martin Lauss (2019). swamp: Visualization, Analysis and 
Adjustment of High-Dimensional Data in Respect to Sample Annotations. 
R package version 1.5.1. \url{https://CRAN.R-project.org/package=swamp}
}

\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
