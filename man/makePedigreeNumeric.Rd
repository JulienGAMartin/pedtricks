\name{makePedigreeNumeric}
\alias{makePedigreeNumeric}
\title{Converts a pedigree with individuals specified as factors to a numeric pedigree}
\description{
  Some internal \code{pedantics} modules require that pedigrees be
specified only by numerical values, including numerical values
for missing data, this provides that conversion
}
\usage{
makePedigreeNumeric(id, sire, dam, missingVal = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{Individual identifiers - pass using \code{as.character()}}
  \item{sire}{Sire codes - pass \code{using as.character()}}
  \item{dam}{Dam codes - pass \code{using as.character()}}
  \item{missingVal}{the indicator that should be substituted for missing values}
}
\value{
  \item{numericPedigree}{The factor pedigree in numeric form}
  \item{idKey}{A key to facilitate conversion back to the original identifiers}
}
\references{Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@st-andrews.ac.uk}}
\seealso{\code{\link{makePedigreeFactor}}}
\examples{
pedigree<-as.data.frame(matrix(c(
"m1",	NA,	NA,
"m2",	NA,	NA,
"m3",	NA,	NA,
"d4",	NA,	NA,
"d5",	NA,	NA,
"o6",	"m1",	"d4",
"o7",	"m1",	"d4",
"o8",	"m1",	"d4",
"o9",	"m1",	"d4",
"o10",	"m2",	"d5",
"o11",	"m2",	"d5",
"o12",	"m2",	"d5",
"o13",	"m2",	"d5",
"o14",	"m3",	"d5",
"o15",	"m3",	"d5",
"o16",	"m3",	"d5",
"o17",	"m3",	"d5"),17,3,byrow=TRUE))
names(pedigree)<-c("id","dam","sire")
for(x in 1:3) pedigree[,x]<-as.factor(pedigree[,x])

## make the test pedigree numeric with NAs denoted by -1
makePedigreeNumeric(id=as.character(pedigree[,1]),
                          dam=as.character(pedigree[,2]),
                          sire=as.character(pedigree[,3]),
                          missingVal=-1)
}
\keyword{manip}
