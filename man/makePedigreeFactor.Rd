\name{makePedigreeFactor}
\alias{makePedigreeFactor}
\title{Converts a numeric pedigee to a pedigree with factors}
\description{
  Some internal \code{pedantics} modules require that pedigrees be
specified only by numerical values, including numerical values
for missing data, this converts them back to factors
}
\usage{
makePedigreeFactor(id, sire, dam, key)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{Numeric individual identifiers}
  \item{sire}{Numeric sire codes}
  \item{dam}{Numeric dam codes}
  \item{key}{A dataframe, as produced by \code{makePedigreeNumeric}, specifying factor codes for numeric values in is, sire, and dam}
}
\value{
  returns the pedigree with all ids specified as factors according to key
}
\references{Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@st-andrews.ac.uk}}
\seealso{\code{\link{makePedigreeNumeric}}}
\examples{

## first we'll implement the example from makePedigreeNumeric(), 
## and use makePedigreeFactor() to turn it back again:

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
test<-makePedigreeNumeric(id=as.character(pedigree[,1]),
                          dam=as.character(pedigree[,2]),
                          sire=as.character(pedigree[,3]),
                          missingVal=-1)

test$numericPedigree

test$idKey

## and turn it back again
makePedigreeFactor(id=test$numericPedigree$id,
                          dam=test$numericPedigree$dam,
                          sire=test$numericPedigree$sire,
                          key=test$idKey)

}
\keyword{manip}
