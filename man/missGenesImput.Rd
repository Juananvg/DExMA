\name{missGenesImput}
\alias{missGenesImput}
\title{Imputation of unmeasured genes}
\usage{
missGenesImput(objectMA, k = 7)
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the different samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateobjectMA can be used too.}

\item{k}{Number of neighbors to be used in the imputation (default=7).}
}
\value{
A list formed by two elements:
\itemize{
\item First element (objectMA) the same objectMA with missign genes imputed
\item Second element (imputIndicators) a list with 4 different objects:
    \itemize{
    \item imputValuesSample: Number of missing values imputed per sample
    \item imputPercentageSample: Percentage of missing values imputed per sample
    \item imputValuesGene: Number of missing values imputed per gene
    \item imputPercentageGene: Percentage of missing values imputed per gene
    }
}
}
\description{
missGenesImput uses k-nearest neighbors in the space of samples
to impute the unmeasured genes of the different datasets.
}
\examples{
data(DExMAExampleData)
missGenesImput(maObject)
}

\references{
Christopher A Mancuso, Jacob L Canfield, Deepak Singla, Arjun Krishnan, 
A flexible, interpretable, and accurate approach for imputing the expression 
of unmeasured genes, Nucleic Acids Research, Volume 48, Issue 21, 
2 December 2020, Page e125, \url{https://doi.org/10.1093/nar/gkaa881}

Alberto Franzin, Francesco Sambo, Barbara di Camillo. bnstruct: an R 
package for Bayesian Network structure learning in the presence of 
missing data. Bioinformatics, 2017; 33 (8): 1250-1252, 
Oxford University Press, \url{https://doi.org/10.1093/bioinformatics/btw807}
}

\seealso{
\code{\link{createObjectMA}}, \code{\link{metaAnalysisDE}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
