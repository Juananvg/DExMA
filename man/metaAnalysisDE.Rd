\name{metaAnalysisDE}
\alias{metaAnalysisDE}
\title{Performing Meta-analysis}
\usage{
metaAnalysisDE(
    objectMA = NULL,
    effectS = NULL,
    pvalues = NULL,
    weight = NULL,
    typeMethod = c("FEM", "REM", "maxP", "minP", "Fisher", 
                    "Stouffer", "ACAT"),
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

\item{effectS}{A list of two elements. The first element is a dataframe
with genes in rows and studies in columns. Each component of the dataframe
is the effect of a gene in a study. 
The second element of the list is also a dataframe
with the same structure, but in this case each component of the dataframe
represent the variance of the effect of a gene in a study. The third element 
of the list is also a dataframe with the same structure, but in this case each 
component of the dataframe represent the log fold change of a gene in a study. 
This argument should be only used in the case that objectMA argument is null.}

\item{pvalues}{A list of two elements. The first element is a dataframe
with genes in rows and studies in columns. Each component of the dataframe
is the p-value of a gene in a study. 
The second element of the list is also a dataframe
with the same structure, but in this case each component of the dataframe
represent the log fold change of a gene in a study. This argument
should be only used in the case that objectMA argument is null.}

\item{weight}{A vector of the weights of each dataset. This argument should 
only be included in case objectMA is null and you want to use "Stouffer" or 
"ACAT" method.}

\item{typeMethod}{A character that indicates the method to be peformed.
See "Details"for more information}

\item{missAllow}{A number that indicates the maximun proportion of missing
values allowed in a sample. If the sample has more proportion of missing
values the sample will be eliminated. In the other case the missing values
will be imputed using the K-NN algorithm. In case the objectMA has been 
previously imputed, this element is not necessary.}

\item{proportionData}{The minimum proportion of datasets in which a gene
must be contained to be included. By default, the gene must be contained 
in at least half of the datasets. In case the objectMA has been previously 
imputed, this element is not necessary.}

}
\value{
A dataframe with the meta-analysis results. Depending on the 
applied method, a different dataframe is obtained. For more information 
see the package vignette.
}
\description{
It performs meta-analysis using seven different methods.
}
\details{
The different meta-analysis methods that can be applied are:
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
        \item "minP": minimum p-value method (Tippett's method)
        \item "ACAT": Aggregated Cauchy Association Test method
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

Liu, Y., Chen, S., Li, Z., Morrison, A. C., Boerwinkle, E., & Lin, X. (2019). 
ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis 
in Sequencing Studies. The American Journal of Human Genetics, 104(3), 410-421.
\url{https://doi.org/10.1016/j.ajhg.2019.01.002}
}

\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
