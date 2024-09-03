#' Manipulating pedigrees to prepare them for requirements of
#' subsequent analyses

#' Prepares a pedigree to conform with requirements of many softwares
#' used in quantitative genetic analysis, as well as for many of the
#' functions in pedtricks.
#' @param ped An ordered pedigree with 3 columns: id, dam, sire
#' @param dat An optional data frame, the same length as the pedigree
#'
#' @return
#'   Returns a pedigree in which all individuals that exist in the dam
#' and sire columns are represented by their own record lines, occurring
#' before the records of their first offspring.  If data are supplied,
#' then fix_ped will return a dataframe, the first three columns are
#' the 'fixed' pedigree, and the following columns of which contain
#' appropriately reordered data.
#'
#' @examples
#' ##  a valid pedigree, i.e., no loops, no bisexuality, etc.,
#' ## but where not all parents have a record line, and where
#' ## parents do not necessarily occur before their offspring:
#' pedigree <- as.data.frame(matrix(c(
#'   10, 1, 2,
#'   11, 1, 2,
#'   12, 1, 3,
#'   13, 1, 3,
#'   14, 4, 5,
#'   15, 6, 7,
#'   4, NA, NA,
#'   5, NA, NA,
#'   6, NA, NA,
#'   7, NA, NA
#' ), 10, 3, byrow = TRUE))
#' names(pedigree) <- c("id", "dam", "sire")
#' pedigree
#' fixed_pedigree <- fix_ped(ped = pedigree)
#' fixed_pedigree
#'
#' @keywords manipulation
#'
#' @export


fix_ped <- function(ped, dat = NULL) {
  if (is.null(dat) == FALSE){
    if (
      is.null(dim(dat)) == FALSE &&
      length(ped[, 1]) != length(dat[, 1])
  ) {
    stop(paste("pedigree and cohorts differ in length.", "\n"))
    flush.console()
    stop()
  }
  if (
    is.null(dim(dat)) &&
      length(ped[, 1]) != length(dat)
  ) {
    stop(paste("pedigree and cohorts differ in length.", "\n"))
    flush.console()
  }
  }

  names(ped) <- c("id", "dam", "sire")
  ntotal <- length(ped$id) * 3
  IDs <- array(dim = ntotal)
  for (x in seq_along(ped$id)) {
    IDs[x] <- as.character(ped$id[x])
    IDs[x + ntotal] <- as.character(ped$dam[x])
    IDs[x + ntotal * 2] <- as.character(ped$sire[x])
  }
  IDs <- as.data.frame(IDs)
  IDs <- unique(IDs)
  IDs <- subset(IDs, is.na(IDs) == FALSE)
  names(IDs) <- "id"
  IDs$dam <- ped$dam[match(IDs$id, ped$id)]
  IDs$sire <- ped$sire[match(IDs$id, ped$id)]

  fixed_ped <- orderPed(IDs)
  if (is.null(dat) == FALSE) {
    if (
        names(dat)[1] == "id" ||
        names(dat)[1] == "ID" ||
        names(dat)[1] == "ids" ||
        names(dat)[1] == "IDS") {
      for (x in 2:length(dat[1, ])) {
        fixed_ped[, (3 + x - 1)] <- dat[match(fixed_ped[, 1], dat[, 1]), x]
      }
    } else {
      warning(
        paste(
          "No id column detected in dat, assuming same order as ped.",
          "\n"
        )
      )
      flush.console()
      dat$id <- ped[, 1]
      for (x in 1:(length(dat[1, ]) - 1)) {
        fixed_ped[, (3 + x - 1)] <- dat[match(fixed_ped[, 1], dat$id), x]
      }
    }
  }
  for (x in 1:3) fixed_ped[, x] <- as.factor(fixed_ped[, x])
  fixed_ped
}

#' @rdname pedantics-deprecated
#' @section \code{fixPedigree}: the function has just been renamed with no
#' other changes for the moment
#' @export
fixPedigree <- function(ped, dat = NULL) {
  .Deprecated(fix_ped,
    msg = "this function from pedantics is deprecated, please use the new
    'fix_ped()' instead",
  )
  fix_ped(ped, dat)
}
