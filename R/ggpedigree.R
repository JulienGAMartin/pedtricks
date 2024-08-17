#' ggpedigree: Pedigree plotting tool for complex pedigrees.
#' @param ids A vector of individual identifiers
#' @param mothers A vector of mothers corresponding to ids. Missing values are 0
#'   or NA.
#' @param fathers A vector of fathers corresponding to ids. Missing values are 0
#'   or NA.
#' @param cohort integer. Default NULL. An optional vector assigning a cohort to
#'   each id. If NULL, then `kinship2::kindepth` is used to assign cohorts to
#'   ids.
#' @param sex integer or character. Default NULL. An optional vector assigning a
#'   sex to each id. This can be any value, but the first level numerically or
#'   alphabetically (e.g. 0 or "F") will be plotted with a circle (traditionally
#'   denoting a female) and the second level will be plotted with a square
#'   (traditionally denoting a male). Any negative values, NA values, or third
#'   alphabetically values will be plotted with a triangle. NOTE: These can be
#'   overridden by specifying `sex_female` or `sex_male` values.
#' @param sex_female integer or character. Default NULL. Indicates the value
#'   used in `sex` for females. Will be plotted as circles.
#' @param sex_male integer or character. Default NULL. Indicates the value used
#'   in `sex` for males. Will be plotted as squares.
#' @param id_labels logical. Default FALSE. Print the ids on the pedigree plot.
#' @param remove_singletons logical. Default TRUE. Remove ids with no relatives
#'   i.e., no offspring or parents assigned.
#' @param plot_unknown_cohort logical. Default TRUE. Plots ids of unknown
#'   cohorts. These are plotted in an "Unknown" cohort at the top of the
#'   pedigree. Be aware that any mothers and fathers of these individuals will
#'   be plotted below them.
#' @param randomise_x_coordinates logical. Default TRUE. Randomise the position of
#'   individuals within each cohort.
#' @param print_cohort_labels logical. Default TRUE. Prints cohort levels on the
#'   left hand side of plot.
#' @param return_plot_tables logical. Default FALSE. Returns an object with the
#'   line and point data used for the plot.
#' @param suppress_plot logical. Default FALSE. Will stop the plot from being
#'   output e.g., if using return_plot_tables to retrieve the plot objects.
#' @param line_col_mother Default = "#E41A1C". Line colour for maternal links.
#' @param line_col_father Default = "#377EB8". Line colour for paternal links.
#' @param line_alpha Default = 0.3. Line alpha (transparency) value for maternal
#'   and paternal links.
#' @param point_size Default = 1. Point size for ids.
#' @param point_colour Default = "black". Point colour for ids.
#' @param point_alpha Default = 1. Point alpha for ids.
#' @param xlab Default = blank. x-axis label.
#' @param ylab  Default = blank. y-axis label
#' @param gg_theme Alternate theme information specified as a ggplot::theme()
#'   object.
#' @import dplyr
#' @import kinship2
#' @import ggplot2
#' @import tidyr
#' @keywords plot
#' @export


ggpedigree <- function(ids = NULL,
                       mothers = NULL,
                       fathers = NULL,
                       cohort = NULL,
                       sex = NULL,
                       sex_female = NULL,
                       sex_male = NULL,
                       id_labels = FALSE,
                       remove_singletons = TRUE,
                       plot_unknown_cohort = TRUE,
                       randomise_x_coordinates = FALSE,
                       print_cohort_labels = TRUE,
                       return_plot_tables = FALSE,
                       suppress_plot = FALSE,
                       line_col_mother = "#E41A1C",
                       line_col_father = "#377EB8",
                       line_alpha = 0.3,
                       point_size = 1,
                       point_colour = "black",
                       point_alpha = 1,
                       xlab = "",
                       ylab = "",
                       gg_theme = NULL) {

  # Check that ids have not been duplicated

  if(any(tabulate(factor(ids)) > 1)) stop("Duplicated values in ids")

  # Check that cohort is an integer

  if(!is.null(cohort)) if(!is.integer(cohort)) stop("Cohort must be an integer")

  # Check that ids, mother, father, cohort (if specified) and sex (if specified) are the same length

  if(length(na.omit(unique(c(length(ids),
                             length(mothers),
                             length(fathers),
                             ifelse(is.null(cohort), NA, length(cohort)),
                             ifelse(is.null(sex), NA, length(sex)))))) > 1){
    stop("ids, mothers, fathers, cohort (if specified), and sex (if specified) must be the same length.")
  }

  # if sex_female is defined, check that sex_male is defined, and vice versa

  if(!is.null(sex_female) &  is.null(sex_male)) stop("if sex_female is defined, sex_male must also be defined (and vice versa)")
  if( is.null(sex_female) & !is.null(sex_male)) stop("if sex_female is defined, sex_male must also be defined (and vice versa)")

  # Format the pedigree to have ID, MOTHER, FATHER columns and recode NA to 0.

  ped <- data.frame(ID = as.character(ids),
                    MOTHER = as.character(mothers),
                    FATHER = as.character(fathers))

  for (i in 1:3) ped[which(is.na(ped[, i])), i] <- 0



  # Add in parents that are not in ids as founders

  baseped <- rbind(data.frame(ID = ped[which(!ped$FATHER %in% ped$ID), "FATHER"],
                              MOTHER = 0,
                              FATHER = 0),
                   data.frame(ID = ped[which(!ped$MOTHER %in% ped$ID), "MOTHER"],
                              MOTHER = 0,
                              FATHER = 0),
                   ped)

  baseped <- unique(subset(baseped, ID != 0))

  # Remove Singletons

  if (remove_singletons) {
    singleton_vec <- which(baseped$MOTHER == 0 & baseped$FATHER == 0 & !baseped$ID %in% c(baseped$MOTHER, baseped$FATHER))
    if (length(singleton_vec) > 0) baseped <- baseped[-singleton_vec, ]
  }

  # Determine cohorts and add to the baseped data.frame if not specified

  if (is.null(cohort)) {
    cohortvec <- data.frame(ID = baseped[, "ID"],
                            Cohort = kindepth(baseped[,"ID"], baseped[,"FATHER"], baseped[,"MOTHER"]))
  } else {
    cohortvec <- data.frame(ID = as.character(ids), Cohort = cohort)
  }

  suppressMessages(baseped <- left_join(baseped, cohortvec))

  # Collate information for axis labels and plotting unknown cohort individuals

  cohort_order <- sort(unique(as.numeric(baseped$Cohort)))
  cohort_labels <- min(cohort_order):max(cohort_order)

  if (plot_unknown_cohort) {
    baseped$Cohort[which(is.na(baseped$Cohort))] <- min(cohort_order) - 1

    if (length(which(baseped$Cohort == min(cohort_order) - 1)) > 0) {
      cohort_order <- c(min(cohort_order) - 1, cohort_order)
      cohort_labels <- c("Unknown", cohort_labels)
    }
  }

  if(!print_cohort_labels) cohort_labels <- rep("", length(cohort_order))

  # Assign x coordinates

  baseped$xcoord <- NA

  for(i in unique(baseped$Cohort)){

    x <- which(baseped$Cohort == i)

    if(length(x) > 1){
      if(randomise_x_coordinates){
        baseped$xcoord[x] <- sample(seq(0, 1, 1/(length(x)-1)))
      } else {
        baseped$xcoord[x] <- seq(0, 1, 1/(length(x)-1))

      }
    } else {
      baseped$xcoord[x] <- 0.5
    }

  }

  # Pivot ped and get rid of connections where value = 0 (means parental connection is unknown)

  baseped2 <- pivot_longer(ped, cols = c("MOTHER", "FATHER"), names_to = "ParentSex", values_to = "ParentID")

  baseped2 <- subset(baseped2, ParentID != 0)

  # Create a group vector for parent/offspring relationship

  baseped2$Group <- 1:nrow(baseped2)

  # Pivot again to create a single line per ID with Group specified

  baseped3 <- pivot_longer(baseped2, cols = c("ID", "ParentID"), names_to = "Relationship", values_to = "ID")

  # Add cohort information

  suppressMessages(baseped3 <- left_join(baseped3, subset(baseped, select = c(ID, Cohort, xcoord))))

  # baseped3 is used for the parental links. Create baseped4 which is a unique
  # value for each individual to avoid overplotting.

  baseped4 <- droplevels(unique(subset(baseped3, select = c(ID, xcoord, Cohort))))

  # If sex is defined, make a data.frame with sex information.

  if(!is.null(sex)){

    sextab <- data.frame(ID = as.character(ids),
                         SEX = sex)

    if(!is.null(sex_female)){

      sex2 <- data.frame(SEX = c(sex_female, sex_male),
                         GraphSEX = c(1, 2))

      suppressMessages(sextab <- left_join(sextab, sex2))
      sextab$GraphSEX[which(is.na(sextab$GraphSEX))] <- 3

      sextab$SEX <- as.character(sextab$GraphSEX)
      sextab$GraphSEX <- NULL

      suppressMessages(baseped4 <- left_join(baseped, sextab))

    } else {

      suppressWarnings(x <- na.omit(unique(sex[as.numeric(sex) < 0])))

      sexlevels <- sort(unique(sex))

      sexlevels <- sexlevels[-which(sexlevels %in% x)]
      if(length(x) > 0) sexlevels <- c(sexlevels, x)

      suppressMessages(baseped4 <- left_join(baseped, sextab))
      baseped4$SEX <- factor(baseped4$SEX, levels = sexlevels)

    }

  } else {

    baseped4$SEX <- "1"

  }

  # Plot the pedigrees

  if (is.null(gg_theme)) {
    gg_theme <- theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = "black"),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white"),
      panel.background = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm")
    )
  }

  if (!suppress_plot) {

    p <- ggplot() +
      geom_line(data = baseped3, aes(x = xcoord,
                                     y = -Cohort,
                                     group = Group,
                                     colour = ParentSex), alpha = line_alpha) +
      scale_y_continuous(
        breaks = -seq(min(cohort_order), max(cohort_order), 1),
        labels = cohort_labels
      ) +
      scale_colour_manual(values = c(line_col_father, line_col_mother)) +
      labs(x = xlab, y = ylab) +
      gg_theme

    if (id_labels) {

      print(p +
              geom_text(
                data = baseped4, aes(x = xcoord, y = -Cohort, label = ID),
                size = point_size, colour = point_colour, alpha = point_alpha
              )
      )

    } else {

      print(p +
              geom_point(
                data = baseped4, aes(x = xcoord, y = -Cohort, shape = SEX),
                size = point_size, colour = point_colour, alpha = point_alpha
              ) +
              scale_shape_manual(values = c(16, 15, 17))
      )
    }


  }

  # Create a return object if return.plot.tables

  if (return_plot_tables) {
    return(list(
      LinePlotFrame = baseped3,
      PointPlotFrame = baseped4
    ))
  }

}
