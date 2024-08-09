\name{fixPedigree}
\alias{fixPedigree}
\title{Manipulating pedigrees to prepare them for requirements of 
subsequent analyses}
\description{
  Prepares a pedigree to conform with requirements of many softwares
used in quantitative genetic analysis, as well as for many of the 
functions in pedantics.
}
\usage{
fixPedigree(Ped, dat = NULL)
}
\arguments{
  \item{Ped}{An ordered pedigree with 3 columns: id, dam, sire}
  \item{dat}{An optional data frame, the same length as the pedigree}
}
\value{
  Returns a pedigree in which all individuals that exsit in the dam 
and sire columns are represented by their own record lines, occurring
before the records of their first offspring.  If data are supplied,
then fixPedigree will return a dataframe, the first two columns are
the 'fixed' pedigree, and the following columns of which contain
appropriately reordered data.
}
\references{Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@ed.ac.uk}}
\examples{
##  a valid pedigree, i.e., no loops, no bisexuality, etc.,
## but where not all parents have a record line, and where 
## parents do not necessarily occur before their offspring:
pedigree<-as.data.frame(matrix(c(
  10,1,2,
  11,1,2,
  12,1,3,
  13,1,3,
  14,4,5,
  15,6,7,
  4,NA,NA,
  5,NA,NA,
  6,NA,NA,
  7,NA,NA),10,3,byrow=TRUE))
names(pedigree)<-c("id","dam","sire")
pedigree
fixedPedigree<-fixPedigree(Ped=pedigree)
fixedPedigree
}
\keyword{manip}
