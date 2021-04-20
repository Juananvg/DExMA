\name{allSameID}
\alias{allSameID}
\title{Set all datasets in the same ID}
\usage{
allSameID(objectMA, initialIDs, finalID = "GeneSymbol",
    organism="Homo sapiens")
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the different samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateobjectMA can be used too.}

\item{initialIDs}{A vector in which each element indicates the ID of the 
equivalent element of listExpMatrix. To know the avalible ids, you can write
avaliableIDs}

\item{finalID}{A character that indicates the final ID all the different
studies will have. To know the available ids, you can write avaliableIDs.}

\item{organism}{A character that indicates the organism of the data.
To know the avaliable organisms write avaliableOrganism}
}
\value{
The same list with all the datasets in the same selected gene ID.
}
\description{
Set all datasets in the 
same ID (Official Gene Symbol, Entrez or Ensembl)
}
\examples{
data(DExMAExampleData)
sameData <- allSameID(objectMA = maObjectDif, initialIDs = c("Entrez",
"Entrez", "GeneSymbol", "GeneSymbol"), finalID = "GeneSymbol", 
organism = "Homo sapiens")
sameData
}
\seealso{
\code{\link{createObjectMA}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juan.villatoro@genyo.es}
}
