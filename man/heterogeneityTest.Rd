\name{heterogeneityTest}
\alias{heterogeneityTest}
\title{It allows checking the heterogeneity of the different studies}
\usage{
heterogeneityTest(objectMA)
}
\arguments{
\item{objectMA}{A list of list. Each list contains two elements. The first
element is the expression matrix (genes in rows and sample in columns) and
the second element is a vector of zeros and ones that represents the state
of the different samples of the expression matrix. 0 represents one group
(controls) and 1 represents the other group (cases).
The result of the CreateobjectMA can be used too.}
}
\value{
Deciles of the I^2 values
}
\description{
Shows a QQ-plot of the Cochran's test and the quantiles of I^2 statistic values
to mesuare heterogeneity
}
\details{
If in the QQ-plot of the Cochran’s test most of the values are close to the 
central line (most of the Cochran’s test values
are close to the expected distribution),
it can be said that there is homogeneity. In the case that these 
values deviate greatly from the expected distribution, it
must be assumed that there is heterogeneity.
I^2 measures the percentage of variation across studies due to heterogeneity. 
To assume homogeneity in the gene expression meta-analysis, 
almost all I^2 values (quantiles) must be 0 or at least less than 0.25. 
}
\examples{

data(DExMAExampleData)

heterogeneityTest(maObject)

}
\references{
Higgins JPT, Thompson SG. Quantifying heterogeneity in a 
meta-analysis. Stat Med 2002;21:1539–58.

Higgins JPT, Thompson SG, Deeks JJ, et al. Measuring inconsistency in 
meta-analyses. BMJ 2003;327:557–60.
}
\seealso{
\code{\link{createObjectMA}}
}
\author{
Juan Antonio Villatoro Garcia, \email{juanantoniovillatorogarcia@gmail.com}
}
