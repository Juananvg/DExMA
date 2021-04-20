\name{downloadGEOData}
\alias{downloadGEOData}
\title{Download datasets from GEO database}
\usage{
downloadGEOData(GEOobject, directory = getwd())
}
\arguments{
\item{GEOobject}{a vector of character where each element represents a
GEO object for downloading.}

\item{directory}{The directory where the different downloaded GSE Series 
Matrix files from GEO will be stored. By default they are downloaded to 
the working directory}
}
\value{
A list of the different ExpressionSets
}
\description{
Download different ExpressionSets objects from the GEO database
}
\details{
This function internally uses getGEO function of GEOquery package. 
However, downloadGEO allows you to download multiple files at the same time.
}
\examples{
\dontrun{
GEOobjects<- c("GSE4588", "GSE10325")
dataGEO<-downloadGEOData(GEOobjects)
dataGEO
}
}
\references{
Davis, S. and Meltzer, P. S. GEOquery: a bridge between the 
Gene Expression Omnibus (GEO) and BioConductor. 
Bioinformatics, 2007, 14, 1846-1847
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
