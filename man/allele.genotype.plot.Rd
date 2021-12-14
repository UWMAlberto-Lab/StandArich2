\name{allele.genotype.plot}
\alias{allele.genotype.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a single plot with population curves of allelic richness vs sample size  }
\description{Uses the Allelic richness table returned by the $A slot from rgenotypes.arich to print population curves of allelic richness as a function of sample size permutations
\usage{
allele.genotype.plot(results, g=0,xmin=0,xmax=50,xmark=10,main="",lwd=1,lty=1,lcol="black",print.pop=FALSE,pop.text.size=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{}
  \item{g}{set the n value to plot a vertical line over the plot to intersect at n. This helps visualize the standardized allelic richness values for a given n}
  \item{xmin}{sets the lower limit of the x axis}
  \item{xmax}{sets the higher limit of the x axis}
  \item{xmark}{sets the interval of tick marks for the x axis. Following the code for seq(xmin,xmax,xmark)}
  \item{main}{String to feed the main argument of plot}
  \item{lwd}{Line width for the lines in the plot. Can take a vector}
  \item{lty}{Line type for the lines in the plot. Can take a vector}
  \item{lcol}{Line color.Can take a vector}
  \item{print.pop}{Logical to plot populations names in front of their lines. Defaults to FALSE}
  \item{pop.text.size}{Controls the size of population names text.}
}
\value{
 Produces a plot with lines showing allelic richness as a function of genet number for each population.
  A vertical line at the g intersection in the x axis is drawn
}
\references{
Alberto F, Arnaud-Haond S, Duarte CM, Serrao EA (in press) Genetic diversity of a clonal angiosperm near its range limit: the case of \emph{Cymodocea nodosa} in the Canary Islands. Marine Ecology Progress Series
}
\author{
Filipe Alberto,
Department of Biological Sciences,
University of Wisconsin-Milwaukee
}

\seealso{
\code{\link{plot.allele.freq}},  \code{\link{allele.genotype.plot}}
}
\examples{
data(Exdata)
Exresults<-rgenotypes.arich(Exdata,3)
allele.genotype.plot(Exresults$R, g=10)

\keyword{ dplot }
