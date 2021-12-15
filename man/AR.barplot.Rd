\name{AR.barplot}
\alias{AR.barplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots a barplot with standardized allelic richness
}
\description{
Plots a barplot with standardized allelic richness for n equal to the smallest population size in the data. Uses the list object produced by \code{\link{rgenotypes.arich}} as input values. 95% Confidence intervals are optional.
}
\usage{
AR.barplot(Rlist,lwd=10,col="black",b.bar.x=5,popsize=0.7,error.bar=F)
}
\arguments{
  \item{Rlist}{a list like the one produced by \code{\link{rgenotypes.arich}} with elements R, A and A.sd}
  \item{lwd}{the width of each bar. The plot used type="h" to produce the bars so this is the lwd used for R line width}
  \item{col}{a vector of colors for each population. Defaults to black}
  \item{b.bar.x}{the space in x plotting units between consecutive bars. Defaults to 5.}
  \item{popsize}{The size for pop names printed in the x axis, controled by par(cex=popsize)}
  \item{error.bar}{A logical to print 95% confidence intervals. Defaults to FALSE}
}
\details{
Once  \code{\link{rgenotypes.arich} produces the randomized data set with permutation of n, this function can be used to plot a barplot of allelic richness for n equal to the smallest population size across all sites. A vector of colors can be used to plot different populations' bars with different color.
}
\value{
a barplot of Allelic richness
}
\references{
Alberto F, Arnaud-Haond S, Duarte CM, Serrao EA (in press) Genetic diversity of a clonal angiosperm near its range limit: the case of \emph{Cymodocea nodosa} in the Canary Islands. Marine Ecology Progress Series}
\author{
Filipe Alberto,
Department of Biological Sciences,
University of Wisconsin-Milwaukee}
\seealso{
\code{\link{rgenotypes.arich}},\code{\link{allele.genotype.plot}
\examples{
data(Exdata)
Exresults<-rgenotypes.arich(Exdata,3)
AR.barplot(Exresults, error.bar=T,lwd=50,b.bar.x = 1)
}
\keyword{ dplot }
