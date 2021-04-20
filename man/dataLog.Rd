\name{dataLog}
\alias{dataLog}
\title{Auxiliary function to check if data are log transfromed and 
transformed if it are not log-transformed}
\usage{
dataLog(objectMA)
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the different samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateobjectMA should be used.}
}
\value{
The same object with log-transformed  expression matrix
}
\description{
Auxiliary function to check if data are log transfromed and 
transformed if it are not log-transformed
}
\examples{

data(DExMAExampleData)

dataLog(maObject)

}
\seealso{
\code{\link{createObjectMA}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
