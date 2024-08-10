#' Converts a pedigree with individuals specified as factors to a numeric
#' pedigree
#'
#' Some internal \code{pedantics} modules require that pedigrees be
#' specified only by numerical values, or including numerical values
#' for missing data. This function provides the conversion to numeric but also
#' back to factors if needed
#'
#' @param type define how to convert the pedigree so "numeric" or "factor"
#' @param id Individual identifiers - pass using \code{as.character()}
#' @param sire Sire codes - pass \code{using as.character()}
#' @param dam Dam codes - pass \code{using as.character()}
#' @param missingVal the indicator that should be substituted for missing values
#' @param key A dataframe, as produced by \code{makePedigreeNumeric}, specifying factor codes for numeric values in is, sire, and dam
#'
#' @return
#'   \item{numericPedigree}{The factor pedigree in numeric form}
#'   \item{idKey}{A key to facilitate conversion back to the original
#' identifiers}
#'
#' @examples
#' pedigree <- as.data.frame(matrix(c(
#'   "m1",   NA,     NA,
#'   "m2",   NA,     NA,
#'   "m3",   NA,     NA,
#'   "d4",   NA,     NA,
#'   "d5",   NA,     NA,
#'   "o6",   "m1",   "d4",
#'   "o7",   "m1",   "d4",
#'   "o8",   "m1",   "d4",
#'   "o9",   "m1",   "d4",
#'   "o10",  "m2",   "d5",
#'   "o11",  "m2",   "d5",
#'   "o12",  "m2",   "d5",
#'   "o13",  "m2",   "d5",
#'   "o14",  "m3",   "d5",
#'   "o15",  "m3",   "d5",
#'   "o16",  "m3",   "d5",
#'   "o17",  "m3",   "d5"
#' ), 17, 3, byrow = TRUE))
#' names(pedigree) <- c("id", "dam", "sire")
#' for (x in 1:3) pedigree[, x] <- as.factor(pedigree[, x])
#'
#' ## make the test pedigree numeric with NAs denoted by -1
#' convert_ped(
#'   type = "numeric",
#'   id = as.character(pedigree[, 1]),
#'   dam = as.character(pedigree[, 2]),
#'   sire = as.character(pedigree[, 3]),
#'   missingVal = -1
#' )
#' @keywords manipulation
#'
#' @export

convert_ped <- function(
    type = "numeric", id, sire, dam, missingVal = NA,
    key = NULL) {
  if (type == "numeric") {
    p <- cbind(id, sire, dam)
    pf <- as.factor(p)
    pn <- as.numeric(pf)
    for (i in 1:length(pn)) {
      if (is.null(missingVal) == FALSE & is.na(pn[i]) == TRUE) {
        pn[i] <- missingVal
      }
    }
    k <- as.data.frame(cbind(as.numeric(pn), as.character(pf)))
    k <- unique(k)
    k[, 1] <- as.numeric(as.character(k[, 1]))
    names(k) <- c("pn", "pf")
    ped <- as.data.frame(matrix(pn, length(id), 3, byrow = FALSE))
    names(ped) <- c("id", "sire", "dam")
    r <- list(numericPedigree = ped, idKey = k)
    return(r)
  }
  if (type == "factor") {
    p <- as.data.frame(as.factor(cbind(id, sire, dam)))
    p$ids <- key$pf[match(p[, 1], key$pn)]
    ped <- as.data.frame(matrix(p$ids, length(id), 3, byrow = FALSE))
    for (x in 1:3) ped[, x] <- as.factor(ped[, x])
    names(ped) <- c("id", "sire", "dam")
    return(ped)
  }
}
