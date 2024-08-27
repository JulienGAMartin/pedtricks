#' Produce a graphical representation of a pedigree
#'
#' Plots a pedigree, with options specific to considerations for pedigrees used
#' to for quantitative genetic inference in natural populations.  Pedigrees
#' containing only those individuals that are informative with respect to
#' (genetic) variation in an arbitrary trait can be plotted, potentially
#' overlain on a complete pedigree.
#' Functions also exist to plot various types of pedigree links associated with
#' focal individuals.
#'
#' @param Ped  An ordered pedigree with 3 columns: id, dam, sire
#' @param cohorts An optional numeric vector of the same length as the pedigree
#'   designating, for example cohort affinities or birth years
#' @param sex An optional numeric vector of the same length as the pedigree
#'   containing the sexes (may be unknown) of all individuals with entries in the
#'   pedigree.  Defaults (modifiable with \code{sexInd}) are 0=male and 1=female
#' @param dat An optional vector or data frame containing indicators of data availability.  If dat contains only ones and zeros, then any individual with any entry of one will be considered as having data records.  If data contains values other than ones and zeros, individuals in the pedigree with rows in data that contain at least one available record, i.e., one data record is not NA, will be treated as having data.
#' @param dots If 'y', then a dot will be printed representing each individual in the pedigree.  If sexes are available, dots will be colour coded by sex.
#' @param plotfull To be used when dat is supplied.  If 'y' (the default), individuals in the pedigree that are uninformative with respect to the available data have their pedigree links plotted in gray.
#' @param writeCohortLabels To be used when cohorts is used.  Will plot the cohort values on the left hand side of the pedigree image.
#' @param links Default is 'all', other values are 'mums' to print only maternal pedigree links and 'dads' to print only paternal pedigree links.
#' @param sexInd To be used with if sex is supplied and if the vector of sex specifiers differs from the default.
#' @param dotSize Set the dot size bigger or smaller
#' @param dataDots Will print dots over the dots denoting individuals, but denoting individuals with available data as indicated by dat.
#' @param dataDots.cex controls the size of dataDots relative to dots.
#' @param cohortLabs.cex controls the size of cohort labels.
#' @param retain When those pedigree links only informative relative to phenotypic data availability are to be plotted, this controls whether or not a pruned pedigree based on phenotypic data is plotted (if set to "pruned"), or whether strictly only those informative pedigree links are plotted (the default)
#' @param focal An optional list containing the id of an individual and the kinds of relatives of the focal individual to which to plot pedigree links.  Available types are 'offspring','descendants','parents',,ancestors', and 'kin'.
#' @param sexColours The colours that will be used to draw points and or lines associated with males and females.
#' @param ... Additional graphical parameters.
#'
#'
#' @author Michael Morrissey \email{michael.morrissey@@st-andrews.ac.uk}
#'
#' @seealso \code{\link{fix_ped}} to prepare pedigrees that may not explicitly contain records for all individuals (i.e., where founding individuals may only appear in the dam or sire column).)
#'
#' @examples
#' data(gryphons)
#' pedigree <- fix_ped(gryphons[, 1:3])
#'
#' ## draw the gryphon pedigree by pedigree depth
#' draw_ped(pedigree)
#'
#' ## draw the gryphon pedigree by cohort
#' draw_ped(pedigree,
#'   cohorts = gryphons$cohort, writeCohortLabels = "y",
#'   cohortLabs.cex = 1
#' )
#'
#' ## draw the gryphon pedigree by cohort with only maternal links
#' draw_ped(pedigree, cohorts = gryphons$cohort, links = "mums")
#'
#' ## draw the gryphon pedigree by cohort with colour only for those
#' ## indiduals that are informative relative to the quantitative
#' ## genetics of a hypothetical trait only measured for individuals
#' ## in the last two cohorts, emphasize the phenotyped individuals
#' ## with large black dots, and all other individuals with dots
#' ## colour coded by sex:
#' dataAvailability <- (gryphons$cohort >= (max(gryphons$cohort) - 1)) + 0
#'
#' draw_ped(pedigree,
#'   cohorts = gryphons$cohort, sex = gryphons$sex,
#'   dots = "y", dat = dataAvailability, writeCohortLabels = "y", dataDots = "y"
#' )
#'
#' @keywords plot
#' @export


draw_ped <- function(
    Ped, cohorts = NULL, sex = NULL, dat = NULL, dots = "n", plotfull = "y",
    writeCohortLabels = "n", links = "all", sexInd = c(0, 1), dotSize = 0.001,
    dataDots = "n", dataDots.cex = 2, cohortLabs.cex = 1,
    retain = "informative", focal = NULL,
    sexColours = c("red", "blue"), ...) {
  for (x in 1:3) Ped[, x] <- as.character(Ped[, x])

  names(Ped)[1] <- "id"
  if (names(Ped)[2] != "dam" | names(Ped)[3] != "sire") {
    if (names(Ped)[3] %in% c("mum", "mom", "mother", "Mum", "Mmom", "Dam", "Mother", "MUM", "MOM", "DAM", "MOTHER")) {
      cat(paste("'mum' appears to be in third column, reordering to 'id','dam','sire'"))
      flush.console()
      Ped <- Ped[, c(1, 3, 2)]
    }
    if (names(Ped)[2] %in% c("mum", "mom", "mother", "Mum", "Mom", "Dam", "Mother", "MUM", "MOM", "DAM", "MOTHER")) {
      names(Ped)[2] <- "dam"
      names(Ped)[3] <- "sire"
    }
    if (names(Ped)[2] != "dam" | names(Ped)[3] != "sire") {
      stop("Unable to identify column names, expecting 'id','dam','sire'")
    }
  }



  if (is.null(cohorts) == FALSE && length(Ped[, 1]) != length(cohorts)) {
    stop("Pedigree and cohorts differ in length.")
  }

  if (is.null(dat) == FALSE && is.numeric(dat) && length(Ped[, 1]) != dim(as.data.frame(dat))[1]) {
    stop("Pedigree and available data differ in length.")
  }

  if (is.null(sex) == FALSE && length(Ped[, 1]) != length(sex)) {
    stop("Pedigree and sex differ in length.")
  }



  if (is.null(cohorts)) cohorts <- kindepth(Ped[, 1], Ped[, 2], Ped[, 3])
  cohorts <- as.numeric(as.character(cohorts))
  Ped$cohorts <- cohorts
  plotWidth <- 1
  plotHeight <- 1
  cohortSizes <- table(Ped$cohorts)
  names(cohortSizes) <- names(table(Ped$cohorts))

  cohortIndex <- array(1, dim = length(cohortSizes))
  names(cohortIndex) <- names(table(Ped$cohorts))

  Ped$xlocs <- NA
  Ped$ylocs <- NA

  scaledCohorts <- NULL

  for (x in 1:length(Ped[, 1])) {
    Ped$xlocs[x] <- cohortIndex[as.character(Ped$cohorts[x])] / (plotWidth * cohortSizes[as.character(Ped$cohort[x])]) - 1 / (2 * plotWidth * cohortSizes[as.character(Ped$cohorts[x])])
    scaledCohorts <- cohorts - min(cohorts)
    scaledCohorts <- scaledCohorts / max(scaledCohorts)
    scaledCohorts <- (scaledCohorts / 1.1) + 0.05
    Ped$ylocs[x] <- plotHeight - scaledCohorts[x]
    cohortIndex[as.character(Ped$cohorts[x])] <- cohortIndex[as.character(Ped$cohorts[x])] + 1
  }
  Ped$xlocs <- Ped$xlocs * 0.94 + 0.03
  if (writeCohortLabels == "y") Ped$xlocs <- Ped$xlocs * 0.95 + 0.05

  Ped$matxlocs <- Ped$xlocs[match(Ped$dam, Ped$id)]
  Ped$matylocs <- Ped$ylocs[match(Ped$dam, Ped$id)]

  Ped$patxlocs <- Ped$xlocs[match(Ped$sire, Ped$id)]
  Ped$patylocs <- Ped$ylocs[match(Ped$sire, Ped$id)]

  if (dots == "y") {
    Ped$dots <- "black"
    if (is.null(sex) == FALSE) {
      for (x in 1:length(sex)) {
        if (sex[x] == sexInd[2]) Ped$dots[x] <- sexColours[1]
        if (sex[x] == sexInd[1]) Ped$dots[x] <- sexColours[2]
      }
    }
  }

  if (dataDots == "y") {
    if (is.numeric(dat)) dataPed <- subset(Ped, rowSums(as.data.frame(as.numeric(dat))) > 0)
    if (is.character(dat)) dataPed <- subset(Ped, Ped[, 1] %in% dat)
  }

  Ped.subset <- NULL

  if (is.null(dat) == FALSE) {
    if (is.numeric(dat)) {
      avail <- rowSums(as.data.frame(dat))
      keep <- subset(Ped$id, avail > 0)
    } else {
      keep <- dat
    }

    if (retain == "informative") {
      Ped.subset <- prune(Ped, keep, make.base = TRUE)
    } else {
      Ped.subset <- prune(Ped, keep, make.base = FALSE)
    }

    Ped.subset$xlocs <- Ped$xlocs[match(Ped.subset$id, Ped$id)]
    Ped.subset$ylocs <- Ped$ylocs[match(Ped.subset$id, Ped$id)]
    Ped.subset$matxlocs <- Ped.subset$xlocs[match(Ped.subset$dam, Ped.subset$id)]
    Ped.subset$matylocs <- Ped.subset$ylocs[match(Ped.subset$dam, Ped.subset$id)]
    Ped.subset$patxlocs <- Ped.subset$xlocs[match(Ped.subset$sire, Ped.subset$id)]
    Ped.subset$patylocs <- Ped.subset$ylocs[match(Ped.subset$sire, Ped.subset$id)]
    Ped.subset$dots <- Ped$dots[match(Ped.subset$id, Ped$id)]
  }

  Ped.focal <- NULL
  if (is.null(focal) == FALSE) {
    Ped.focal <- subset(Ped, Ped$id == as.character(focal[1]))
    if (focal[2] == "offspring") {
      o <- subset(Ped, as.character(Ped$dam) == as.character(focal[1]) |
        as.character(Ped$sire) == as.character(focal[1]))
      Ped.focal <- rbind(Ped.focal, Ped[which(Ped$id %in% o$id), ])
    }
    if (focal[2] == "parents") {
      f <- subset(Ped, Ped$id == as.character(focal[1]))
      if (is.na(f$sire) == FALSE) Ped.focal <- rbind(Ped.focal, subset(Ped, Ped$id == f$sire))
      if (is.na(f$dam) == FALSE) Ped.focal <- rbind(Ped.focal, subset(Ped, Ped$id == f$dam))
    }


    if (focal[2] == "ancestors" | focal[2] == "kin") {
      for (x in 1:10) {
        Ped.focal <- rbind(Ped.focal, Ped[which(Ped$id %in% Ped.focal$dam), ])
        Ped.focal <- rbind(Ped.focal, Ped[which(Ped$id %in% Ped.focal$sire), ])
      }
    }


    if (focal[2] == "descendants" | focal[2] == "kin") {
      for (x in 1:10) {
        Ped.focal <- rbind(Ped.focal, Ped[which(Ped$dam %in% Ped.focal$id), ])
        Ped.focal <- rbind(Ped.focal, Ped[which(Ped$sire %in% Ped.focal$id), ])
      }
    }




    Ped.focal$matxlocs.focal <- Ped.focal$xlocs[match(Ped.focal$dam, Ped.focal$id)]
    Ped.focal$matylocs.focal <- Ped.focal$ylocs[match(Ped.focal$dam, Ped.focal$id)]

    Ped.focal$patxlocs.focal <- Ped.focal$xlocs[match(Ped.focal$sire, Ped.focal$id)]
    Ped.focal$patylocs.focal <- Ped.focal$ylocs[match(Ped.focal$sire, Ped.focal$id)]
  }

  if (is.null(dat) == FALSE) {
    cat(paste("Individuals in full pedigree:", length(Ped[, 1]), "\n"))
    flush.console()
    cat(paste("Individuals in informative pedigree subset:", length(Ped.subset[, 1]), "\n"))
    flush.console()
  }

  if (is.null(dat) & is.null(focal)) {
    if (links == "all" | links == "mums") {
      grid.segments(Ped$xlocs, Ped$ylocs, Ped$matxlocs, Ped$matylocs, gp = gpar(col = sexColours[1]))
    }
    if (links == "all" | links == "dads") {
      grid.segments(Ped$xlocs, Ped$ylocs, Ped$patxlocs, Ped$patylocs, gp = gpar(col = sexColours[2]))
    }
    if (dots == "y") for (x in 10:0) grid.circle(Ped$xlocs, Ped$ylocs, r = dotSize, gp = gpar(col = Ped$dots, fill = Ped$dots))
  } else {
    if (plotfull == "y" & sexColours[1] == "red") {
      if (links == "all" | links == "mums") {
        grid.segments(Ped$xlocs, Ped$ylocs, Ped$matxlocs, Ped$matylocs, gp = gpar(col = "gray"))
      }
      if (links == "all" | links == "dads") {
        grid.segments(Ped$xlocs, Ped$ylocs, Ped$patxlocs, Ped$patylocs, gp = gpar(col = "gray"))
      }
    }
    if (plotfull == "y" & sexColours[1] != "red") {
      if (links == "all" | links == "mums") {
        grid.segments(Ped$xlocs, Ped$ylocs, Ped$matxlocs, Ped$matylocs, gp = gpar(col = colours()[354]))
      }
      if (links == "all" | links == "dads") {
        grid.segments(Ped$xlocs, Ped$ylocs, Ped$patxlocs, Ped$patylocs, gp = gpar(col = colours()[354]))
      }
    }
    if (is.null(dat) == FALSE & (links == "all" | links == "mums")) {
      grid.segments(Ped.subset$xlocs, Ped.subset$ylocs, Ped.subset$matxlocs, Ped.subset$matylocs, gp = gpar(col = sexColours[1]))
    }
    if (is.null(dat) == FALSE & (links == "all" | links == "dads")) {
      grid.segments(Ped.subset$xlocs, Ped.subset$ylocs, Ped.subset$patxlocs, Ped.subset$patylocs, gp = gpar(col = sexColours[2]))
    }

    if (is.null(focal) == FALSE & (links == "all" | links == "mums")) {
      grid.segments(Ped.focal$xlocs, Ped.focal$ylocs, Ped.focal$matxlocs.focal, Ped.focal$matylocs.focal, gp = gpar(col = sexColours[1]))
    }
    if (is.null(focal) == FALSE & (links == "all" | links == "dads")) {
      grid.segments(Ped.focal$xlocs, Ped.focal$ylocs, Ped.focal$patxlocs.focal, Ped.focal$patylocs.focal, gp = gpar(col = sexColours[2]))
    }

    if (dots == "y") grid.circle(Ped$xlocs, Ped$ylocs, r = dotSize, gp = gpar(col = Ped$dots, fill = Ped$dots))
    if (dataDots == "y") grid.circle(dataPed$xlocs, dataPed$ylocs, r = dotSize * dataDots.cex, gp = gpar(col = "black", fill = "black"))
  }
  if (writeCohortLabels == "y") {
    grid.text(names(table(cohorts)), x = rep(0.05, length(table(cohorts))), y = plotHeight - as.numeric(names(table(scaledCohorts))), gp = gpar(cex = cohortLabs.cex))
  }
}

#' @rdname pedantics-deprecated
#' @section \code{drawPedigree}: the function has just been renamed with no other changes for the moment, but will soon be replace by "ggpedigree" function
#' @export
drawPedigree <- function(
    Ped, cohorts = NULL, sex = NULL, dat = NULL, dots = "n", plotfull = "y",
    writeCohortLabels = "n", links = "all", sexInd = c(0, 1), dotSize = 0.001,
    dataDots = "n", dataDots.cex = 2, cohortLabs.cex = 1,
    retain = "informative", focal = NULL,
    sexColours = c("red", "blue"), ...) {
  .Deprecated(fix_ped,
    msg = "this function from pedantics is deprecated, please use the new 'draw_ped()' instead",
  )
  draw_ped(
    Ped, cohorts, sex, dat, dots, plotfull, writeCohortLabels, links,
    sexInd, dotSize,
    dataDots, dataDots.cex, cohortLabs.cex,
    retain, focal,
    sexColours, ...
  )
}
