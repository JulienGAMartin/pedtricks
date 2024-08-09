 \name{genotype.list}
\alias{genotype.list}
\title{Genotype Objects for all Loci}
\description{
  Creates a \code{list} of \code{genotype} objects from a \code{matrix} or \code{data.frame} of multilocus genotypes.
}
\usage{
genotype.list(G, marker.type="MSW")
}
\arguments{
  \item{G}{matrix or data.frame of multilocus genotypes with individuals down the rows and loci across columns. Adjacent columns are taken to be the same locus}
  \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).}
}
\value{
  list of \code{genotype} objects for all loci
}
\references{
Marshall, T. C. \emph{et al} (1998) Molecular Ecology 7 5 639-655
Kalinowski S.T. \emph{et al} (2007) Molecular Ecology 16 5 1099-1106
Hadfield J. D.  \emph{et al} (2009) \emph{in prep}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{genotype}}
\examples{
\dontrun{
  data(WarblerG)

  G<-genotype.list(WarblerG[,-1])
  summary(G[[1]])
}
}
\keyword{manip}