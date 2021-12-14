\name{rgenotypes.arich}
\alias{rgenotypes.arich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Estimates allelic richness from  population data sets with varying number of resampled genotypes (g)

}
\description{
Permutes the data from n=1 to n=N (N sample size for each pop) for each population and for each n a number of times set by number of replicates. For each replicate the number of alleles for each locus in the permuted sample is returned as well as the mean and st dev across all loci. Then a tables with mean allelic richness across replicates for each n for each population is returned, as well as the st dev across replicates for each n for each population  summarizes the}
\usage{
rgenotypes.arich((Popdata, n.replicates,loc.labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Popdata}{A data table with all populations on top of each other. Lines correspond
to individual genotypes. The first column has the population codes, the 2nd column the individual codes. The following
columns have the allele codes, in integers, using three or two digits (but be consistent!). Missing data is
represented by the value 999.There is no header in this file}
  \item{n.replicates}{The number of replicates used for each n resampling step}
  \item{loc.labels}{A character vector with the names of loci used in the order found in Popdata}

}
\details{
In order to estimate population allelic richness that it comparable among populations,
the mean number of alleles per locus is estimated for  population data sets with varying number of
resampled genotypes (g), for all possible n values from n=1 to n=N, N being the total number of genotypes
observed in each population.

Standardization of N to an equal n for all populations  is necessary to obtain comparable estimates across populations.
The multiple random reduction (Leberg 2002b) employed here gives similar estimates to the rarefraction method (Petit et al. 1998)

With StandArich.2 this functions removes the need to use stand.Arich function after using rgenotypes.arich. All results are now returned in a list
}
\value{
The function returns a list with three elements R, A and stdA. R is a table containing information for all permutations, showing number of alleles observed at each n per loci, the mean number of alleles observed across loci (i.e., allelic richness per permutation), and the standard deviation of the number of alleles observed across loci. A is a table with the mean allelic richness across all permutations for each level of n for each population. Population data is organized by rows and different n values by columns. Finally, stdA is a table with the standard deviation across all permutations for each level of n for each population. The table format is identical to A, population by rows and n by columns. To get the allelic richness for a particular value of n for all populations, often the smallest sample size across all populations, just extract the column with corresponding n index from A.
}
\references{
Alberto F, Arnaud-Haond S, Duarte CM, Serrao EA (in press) Genetic diversity of a clonal angiosperm near
its range limit: the case of \emph{Cymodocea nodosa}in the Canary Islands. Marine Ecology Progress Series

Leberg PL (2002) Estimating allelic richness: Effects of sample size and bottlenecks. Mol Ecol 11: 2445-2449

Petit RJ, El Mousadik A, Pons O (1998) Identifying populations for conservation on the basis of genetic markers.
Conservation Biology, 12: 844-855
}
\author{
Filipe Alberto,
Department of Biological Sciences,
University of Wisconsin-Milwaukee
\email{albertof@uwm.edu}
}
\seealso{
\code{\link{plot.allele.freq}},  \code{\link{allele.genotype.plot}}
}
\examples{
data(Exdata)
rgenotypes.arich(Exdata,g=10,loc.labels=paste("L",1:8,sep="."))
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
\keyword{manip}
