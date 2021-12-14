\name{plot.allele.freq}
\alias{plot.allele.freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots population allele frequencies per locus represented by dots of varying size}
\description{For each locus a plot is constructed in a new window. Allele frequencies are represented by
  dots of varying diameter. Allele codes are indicated on the x axis and populations on the y axis
}
\usage{
plot.allele.freq(data, LocusNames,Plogbase=20,textlegend=0.8,freqlegend=c(0.1,0.15,0.2,0.5,0.75,1),
                psize=7,poptext=0.5,allelesize=1,xaxispos=0,yplotlim=-3,alleleangle=45,
                pdfout=FALSE,pdfname="out.pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data.frame with rows corresponding to individuals and populations on top of each other.
   The first column has the population codes, the 2nd the individual codes. The following columns
  have the allele codes, in integers, using three or two digits  (but be consistent!).
  Missing data is represented by the value 999. There is no   header in this file}
  \item{LocusNames}{a character vector with the names of loci}
  \item{Plogbase}{used to control the points cex argument to plot allele frequencies as points, as in cex=log(allelefreq+1,base=Plogbase)*psize)}
  \item{textlegend}{controls the cex argument for the legend text size}
  \item{freqlegend}{a numeric vector with the allele frequencies used in the legend}
  \item{poptext}{controls the cex argument for the font size of population names in the y-axis}
  \item{allelesize}{controls the cex argument for the font size of allele names in the x-axis}
  \item{xaxispos}{the pos argument for horizontal axis in the plot, see help(axis)}
  \item{yplotlim}{The ylim argument of the plot, see help(plot.default)}
  \item{alleleangle}{text rotation for allele names, in degrees, defaults to 45}
  \item{pdfout}{logical to control if a pdf with the plot should be written, defaults to FALSE. Currently the function to export pdf is not working in RStudio}
  \item{psize}{Used to control the points cex argument to plot allele frequencies as points, as in cex=log(allelefreq+1,base=Plogbase)*psize}
  \item{pdfname}{Name for the output pdf file to write the plots if pdfout is TRUE, defaults to out.pdf}

}
\details{
Depending on the number of populations and alleles in the data set, different combinations of psize and Plogbase arguments should be tried}
\value{
As many different plot windows as there are loci in the data file, allowing a visualization  of allele frequencies at each locus across populations
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
allele.freq.plot(Exdata,LocusNames=paste("L",1:8,sep="."))
% Use only one concept per line.
\keyword{ dplot }
