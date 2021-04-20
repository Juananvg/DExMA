\name{createObjectMA}
\alias{createObjectMA}
\title{Creation of the object to use in meta-analysis}
\usage{
createObjectMA(
    listEX,
    listPheno = NULL,
    namePheno = c(rep(1, length(listEX))),
    expGroups = c(rep(1, length(listEX))),
    refGroups = c(rep(2, length(listEX)))
)
}
\arguments{
\item{listEX}{A list of dataframes or matrix (genes in rows and sample
in columns). A list of ExpressionSets can be used too}

\item{listPheno}{A list of phenodatas (dataframes or matrix). If the object
listEX is a list of ExpressionSets this element can be null.}

\item{namePheno}{A list or vector of the different colunm names or 
positions from the phenodatas where the experimental and reference groups 
are identified. Each element of namePheno correspont to its equivalent 
element in  the listPheno (default a vector of 1, all the first columns of 
each elements of listPheno are selected).}

\item{expGroups}{A list of vectors or a vector containing the names or the 
positions with which we identify the elements of the experiment groups 
(cases) of the namePheno element (default a vector  of 1, all the first 
groups are selected)}

\item{refGroups}{A list of vectors or a vector containing the names or the 
positions with which we identify the elements of the reference groups 
(control) of the namePheno elements (default a vector  of 1, all the first 
groups are selected)}
}
\value{
The object needed to perform meta-analysis. This object is list of nested lists.
Each list contains two elements:
\itemize{
    \item The first element is the expression matrix (genes in rows and 
    sample in columns)
    \item The second element is a vector of zeros and ones that represents 
    the state of the diffenrent samples of the expression matrix. 
    0 represents reference group (controls) and 1 represents 
    experimental group (cases).
}
}
\description{
It allows the creation of an object to perform meta-analysis.
}
\examples{
data(DExMAExampleData)

phenoGroups = c("condition", "condition", "state", "state")
phenoCases = list(Study1 = "Diseased", Study2 = c("Diseased", "ill"),
                    Study3 = "Diseased", Study4 = "ill")
phenoControls = list(Study1 = "Healthy", Study2 = c("Healthy", "control"),
                        Study3 = "Healthy", Study4 = "control")

newObjectMA <- createObjectMA(listEX=listMatrixEX, listPheno = listPhenodatas, 
                                namePheno=phenoGroups, expGroups=phenoCases, 
                                refGroups = phenoControls)
newObjectMA

}
\seealso{
\code{\link{elementObjectMA}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
