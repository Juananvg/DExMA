\name{metaAnalysisDE}
\alias{metaAnalysisDE}
\title{Performing Meta-analysis}
\usage{
metaAnalysisDE(
    objectMA,
    typeMethod = c("FEM", "REM", "maxP", "minP", "Fisher", 
                    "Stouffer"),
    missAllow = 0.3,
    proportionData = 0.5
)
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the diffenrent samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateobjectMA can be used too.}

\item{typeMethod}{A character that indicates the method to be peformed.
See "Details"for more information}

\item{missAllow}{a number that indicates the maximun proportion of missing
values allowed in a sample. If the sample has more proportion of missing
values the sample will be eliminated. In the other case the missing values
will be imputed using the K-NN algorithm.}

\item{proportionData}{The minimum proportion of datasets in which a gene
must be contained to be included. By default, the gene must be contained 
in at least half of the datasets}

}
\value{
A dataframe with the meta-analysis results. Depending on the 
applied method, a different dataframe is obtained. For more information 
see the package vignette.
}
\description{
It performs meta-analysis using eight different methods.
}
\details{
The different meta-analysis method that can be applied are:
\enumerate{
\item \bold{Effects sizes methods:}
    \itemize{
        \item "FEM": Fixed Effects model
        \item "REM": Random Effects model
        }
\item \bold{P-value combination mehods}
    \itemize{
        \item "Fisher": Fisher's methods
        \item "Stouffer": Stouffer's method
        \item "maxP": maximum p-value method (Wilkinson's method)
        \item "minP": minimum p-value method (Tippet's method)
        }
}
}
\examples{
data(DExMAExampleData)
ResultsMA <- metaAnalysisDE(objectMA=maObject, typeMethod="REM",
                            missAllow=0.3, proportionData=0.5)
ResultsMA
}
\references{
Daniel Toro-Domínguez, Juan Antonio Villatoro-García, 
Jordi Martorell-Marugán, Yolanda Román-Montoya, Marta E Alarcón-Riquelme, 
Pedro Carmona-Sáez, 
A survey of gene expression meta-analysis: methods and applications, 
Briefings in Bioinformatics, 2020;, bbaa019, 
\url{https://doi.org/10.1093/bib/bbaa019}

Michael Dewey (2020). metap: meta-analysis of significance values.
}
\author{
Juan Antonio Villatoro Garcia, \email{juan.villatoro@genyo.es}
}
