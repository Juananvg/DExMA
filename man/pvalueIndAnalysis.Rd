\name{pvalueIndAnalysis}
\alias{pvalueIndAnalysis}
\title{Calculation p-value for each gene and study}
\usage{
pvalueIndAnalysis(objectMA, missAllow = 0.3)
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the diffenrent samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateObjectMA can be used too.}

\item{missAllow}{a number that indicates the maximun proportion of missing
values allowed in a sample. If the sample has more proportion of missing
values the sample will be eliminated. In the other case the missing values
will be imputed using the K-NN algorithm.}
}
\value{
A list formed by two elements:
\itemize{
\item First element (p) is a dataframe were columns are each of the studies
(datasets) and rows are the genes. Each element of the dataframe represents
the p-value.
\item Second element (FC) is a dataframe were columns are each of the studies
(datasets) and rows are the genes. Each element of the dataframe is the logFC.
}
}
\description{
This function uses t-test based on limma package in other to obtain the
individual p-values for each study and gene
}
\examples{

data(DExMAExampleData)

pvalues <- pvalueIndAnalysis(objectMA=maObject, missAllow=0.3)
pvalues

}
\seealso{
\code{\link{createObjectMA}}, \code{\link{metaAnalysisDE}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
