\name{calculateES}
\alias{calculateES}
\title{Calculation of Effects Sizes and their variance}
\usage{
calculateES(objectMA, missAllow = 0.3)
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the diffenrent samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateobjectMA can be used too.}

\item{missAllow}{a number that indicates the maximun proportion of missing
values allowed in a sample. If the sample has more proportion of missing
values the sample will be eliminated. In the other case the missing values
will be imputed using the K-NN algorithm.}
}
\value{
A list formed by three elements:
First element (ES) is a dataframe were columns are each of the studies
(datasets) and rows are the genes. Each element of the dataframe represents
the Effect Size.
Second element (Var) is a dataframe were columns are each of the studies
(datasets) and rows are the genes. Each element of the dataframe represents
the variance of the Effect size.
}
\description{
This function uses the Hedges'g estimator to calulate the different Effects
size and their variances for each genes and for each dataset.
}
\examples{

data(DExMAExampleData)

resultsEffects <- calculateES(objectMA = maObject, missAllow = 0.3)
resultsEffects

}
\seealso{
\code{\link{createObjectMA}}, \code{\link{metaAnalysisDE}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
