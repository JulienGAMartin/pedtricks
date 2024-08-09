\name{extractA}
\alias{extractA}
\title{Allele Frequencies}
\description{extracts allele frequencies from genotype data}
\usage{extractA(G, marker.type="MSW")}
\arguments{
  \item{G}{data frame or list of \code{genotype} objects}
  \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Marshall, 1998) or \code{"AFLP"} for dominant markers.}
}
\value{list of allele frequnecies at each loci}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{genotype.list}}, \code{genotype}}
\examples{
\dontrun{
  data(WarblerG)

  A<-extractA(WarblerG)
  A[[1]]
}
}
\keyword{manip}