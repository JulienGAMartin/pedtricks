% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{extractA}
\alias{extractA}
\title{Allele Frequencies}
\usage{
extractA(G, marker.type = "MSW")
}
\arguments{
\item{G}{data frame or list of \code{genotype} objects}

\item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Marshall, 1998) or \code{"AFLP"} for dominant markers.}
}
\value{
list of allele frequencies at each loci
}
\description{
Extracts allele frequencies from genotype data
}
\examples{
\donttest{
data(WarblerG)

A <- extractA(WarblerG)
A[[1]]
}

}
\author{
Jarrod Hadfield \email{j.hadfield@ed.ac.uk}
}
\keyword{internal}
