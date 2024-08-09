\name{pedStatSummary}
\alias{pedStatSummary}
\title{Post-processes output from pedigreeStats}
\description{
  Generates a manageable summary of pedigree-wide statistics
reported by pedigreeStats, either for a single pedigree or
for a comparison between pedigrees
}
\usage{
pedStatSummary(pedStats,pedStats2=NULL)
}
\arguments{
  \item{pedStats}{An output data list from pedgireeStats}
  \item{pedStats2}{An optional output data list from pedgireeStats}
}
\value{
  Returns a table of numbers of records, maternities,
paternities, pairwise sibship relationships, numbers
of different classes of grand-parental relationships,
pedigree depth, number of founders, mean subship sizes,
simple statistics of numbers of inbred and non-inbred
indivdiuals, and proportions of pairwise relationship 
coefficients equal to or greater than several
thresholds.
}
\references{Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@st-andrews.ac.uk}}
\keyword{manip}
