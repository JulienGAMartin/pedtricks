#' ggpedigree: Pedigree plotting tool for complex pedigrees.
#'
#' @param .data an optional data frame object with all the pedigree information
#' @param ids a column of .data or a vector of individual identifiers
#' @param mothers A column of .data or a vector of mothers corresponding to ids. Missing values are 0
#'   or NA.
#' @param fathers A column of .data or a vector of fathers corresponding to ids. Missing values are 0
#'   or NA.
#' @param cohort integer. Default NULL. A column of .data or an optional vector assigning a cohort to
#'   each id. If NULL, then `kinship2::kindepth` is used to assign cohorts to
#'   ids.
#' @param sex integer or character. Default NULL. An optional column of .data or a vector assigning a
#'   sex to each id. This can be any value, but the first level numerically or
#'   alphabetically (e.g. 0 or "F") will be plotted with a circle (traditionally
#'   denoting a female) and the second level will be plotted with a square
#'   (traditionally denoting a male). Any negative values, NA values, or third
#'   alphabetically values will be plotted with a triangle. NOTE: These can be
#'   overridden by specifying `sex_female` or `sex_male` values.
#' @param sex_code Default NULL. A vector of length 2, indicating the value
#'   used in `sex` for females and males respectively. Females will be plotted
#'   as circles and males as squares.
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
#'   line and point data used for the plot, but the plot will not be generated
#' @param line_col_mother Default = "#E41A1C". Line colour for maternal links.
#' @param line_col_father Default = "#377EB8". Line colour for paternal links.
#' @param line_alpha Default = 0.3. Line alpha (transparency) value for maternal
#'   and paternal links.
#' @param point_size Default = 1. Point size for ids.
#' @param point_colour Default = "black". Point colour for ids.
#' @param point_alpha Default = 1. Point alpha for ids.
#'
#' @examples
#' data(gryphons)
#' pedigree <- fix_ped(gryphons[, 1:3])
#'
#' ## draw the gryphon pedigree by pedigree depth
#' ggpedigree(pedigree)
#' \dontrun{
#' # specifying the column names for id, mother and father
#' ggpedigree(pedigree, id, dam, sire)
#'
#' # with cohort and sex
#' ggpedigree(gryphons, cohort = cohort, sex = sex, sex_code = c(0, 1))
#' }
#' @keywords plot
#' @export


ggpedigree <- function(.data,
                       ids,
                       mothers,
                       fathers,
                       cohort,
                       sex,
                       sex_code = NULL,
                       id_labels = FALSE,
                       remove_singletons = TRUE,
                       plot_unknown_cohort = TRUE,
                       randomise_x_coordinates = FALSE,
                       print_cohort_labels = TRUE,
                       return_plot_tables = FALSE,
                       line_col_mother = "#E41A1C",
                       line_col_father = "#377EB8",
                       line_alpha = 0.3,
                       point_size = 1,
                       point_colour = "black",
                       point_alpha = 1) {
  if (hasArg(.data)) {
    if (!hasArg(ids) || !hasArg(mothers) || !hasArg(fathers)) {
      warning(
        "the first 3 columns were used for id, dam and sire identity, please specify if not correct or to remove the warning"
      )
    }
    if (hasArg(ids)) {
      ids <- as.vector(select(.data, {{ ids }}))[[1]]
    } else {
      ids <- as.vector(select(.data, 1))[[1]]
    }
    if (hasArg(mothers)) {
      mothers <- as.vector(select(.data, {{ mothers }}))[[1]]
    } else {
      mothers <- as.vector(select(.data, 2))[[1]]
    }
    if (hasArg(fathers)) {
      fathers <- as.vector(select(.data, {{ fathers }}))[[1]]
    } else {
      fathers <- as.vector(select(.data, 3))[[1]]
    }
    if (hasArg(sex)) {
      sex <- as.vector(select(.data, {{ sex }}))[[1]]
    } else {
      sex <- NULL
    }
    if (hasArg(cohort)) {
      cohort <- as.vector(select(.data, {{ cohort }}))[[1]]
    } else {
      cohort <- NULL
    }
  }
  # Check that ids have not been duplicated

  if (any(tabulate(factor(ids)) > 1)) stop("Duplicated values in ids")

  # Check that cohort is an integer

  if (!is.null(cohort)) if (!is.integer(cohort)) stop("Cohort must be an integer")

  # Check that ids, mother, father, cohort (if specified) and sex (if specified) are the same length

  if (length(na.omit(unique(c(
    length(ids),
    length(mothers),
    length(fathers),
    ifelse(is.null(cohort), NA, length(cohort)),
    ifelse(is.null(sex), NA, length(sex))
  )))) > 1) {
    stop("ids, mothers, fathers, cohort (if specified), and sex (if specified) must be the same length.")
  }

  # if sex_female is defined, check that sex_male is defined, and vice versa

  if (!is.null(sex) & is.null(sex_code)) stop("sex_code should be provided when sex is specified")
  if (!is.null(sex_code) & length(sex_code) != 2) stop("sex_code should be a vector of length 2 providing the coding for females and then males")

  # Format the pedigree to have ID, MOTHER, FATHER columns and recode NA to 0.

  ped <- data.frame(
    ID = as.character(ids),
    MOTHER = as.character(mothers),
    FATHER = as.character(fathers)
  )

  for (i in 1:3) ped[which(is.na(ped[, i])), i] <- 0



  # Add in parents that are not in ids as founders

  baseped <- rbind(
    data.frame(
      ID = ped[which(!ped$FATHER %in% ped$ID), "FATHER"],
      MOTHER = 0,
      FATHER = 0
    ),
    data.frame(
      ID = ped[which(!ped$MOTHER %in% ped$ID), "MOTHER"],
      MOTHER = 0,
      FATHER = 0
    ),
    ped
  )

  baseped <- unique(filter(baseped, .data$ID != 0))

  # Remove Singletons

  if (remove_singletons) {
    singleton_vec <- which(baseped$MOTHER == 0 & baseped$FATHER == 0 & !baseped$ID %in% c(baseped$MOTHER, baseped$FATHER))
    if (length(singleton_vec) > 0) baseped <- baseped[-singleton_vec, ]
  }

  # Determine cohorts and add to the baseped data.frame if not specified

  if (is.null(cohort)) {
    cohortvec <- data.frame(
      ID = baseped[, "ID"],
      Cohort = kindepth(baseped[, "ID"], baseped[, "FATHER"], baseped[, "MOTHER"])
    )
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

  if (!print_cohort_labels) cohort_labels <- rep("", length(cohort_order))

  # Assign x coordinates

  baseped$xcoord <- NA

  for (i in unique(baseped$Cohort)) {
    x <- which(baseped$Cohort == i)

    if (length(x) > 1) {
      if (randomise_x_coordinates) {
        baseped$xcoord[x] <- sample(seq(0, 1, 1 / (length(x) - 1)))
      } else {
        baseped$xcoord[x] <- seq(0, 1, 1 / (length(x) - 1))
      }
    } else {
      baseped$xcoord[x] <- 0.5
    }
  }

  # Pivot ped and get rid of connections where value = 0 (means parental connection is unknown)

  baseped2 <- pivot_longer(ped, cols = c("MOTHER", "FATHER"), names_to = "ParentSex", values_to = "ParentID")

  baseped2 <- filter(baseped2, .data$ParentID != 0)

  # Create a group vector for parent/offspring relationship

  baseped2$Group <- 1:nrow(baseped2)

  # Pivot again to create a single line per ID with Group specified

  baseped3 <- pivot_longer(baseped2, cols = c("ID", "ParentID"), names_to = "Relationship", values_to = "ID")

  # Add cohort information

  suppressMessages(baseped3 <- left_join(baseped3, select(baseped, "ID", "Cohort", "xcoord")))

  # baseped3 is used for the parental links. Create baseped4 which is a unique
  # value for each individual to avoid overplotting.

  baseped4 <- droplevels(unique(select(baseped3, "ID", "xcoord", "Cohort")))

  # If sex is defined, make a data.frame with sex information.

  if (!is.null(sex)) {
    sextab <- data.frame(
      ID = as.character(ids),
      SEX = sex
    )

    if (!is.null(sex_code)) {
      sex_female <- sex_code[1]
      sex_male <- sex_code[2]
      sex2 <- data.frame(
        SEX = c(sex_female, sex_male),
        GraphSEX = c(1, 2)
      )

      suppressMessages(sextab <- left_join(sextab, sex2))
      sextab$GraphSEX[which(is.na(sextab$GraphSEX))] <- 3

      sextab$SEX <- as.character(sextab$GraphSEX)
      sextab$GraphSEX <- NULL

      suppressMessages(baseped4 <- left_join(baseped, sextab))
    } else {
      suppressWarnings(x <- na.omit(unique(sex[as.numeric(sex) < 0])))

      sexlevels <- sort(unique(sex))

      sexlevels <- sexlevels[-which(sexlevels %in% x)]
      if (length(x) > 0) sexlevels <- c(sexlevels, x)

      suppressMessages(baseped4 <- left_join(baseped, sextab))
      baseped4$SEX <- factor(baseped4$SEX, levels = sexlevels)
    }
  } else {
    baseped4$SEX <- "1"
  }

  # Plot the pedigrees

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

  if (!return_plot_tables) {
    p <- ggplot() +
      geom_line(data = baseped3, aes(
        x = .data$xcoord,
        y = -.data$Cohort,
        group = .data$Group,
        colour = .data$ParentSex
      ), alpha = line_alpha) +
      scale_y_continuous(
        breaks = -seq(min(cohort_order), max(cohort_order), 1),
        labels = cohort_labels
      ) +
      scale_colour_manual(values = c(line_col_father, line_col_mother)) +
      labs(x = "", y = "") +
      gg_theme

    if (id_labels) {
      return(p +
        geom_text(
          data = baseped4, aes(x = .data$xcoord, y = -.data$Cohort, label = .data$ID),
          size = point_size, colour = point_colour, alpha = point_alpha
        ))
    } else {
      return(p +
        geom_point(
          data = baseped4, aes(x = .data$xcoord, y = -.data$Cohort, shape = .data$SEX),
          size = point_size, colour = point_colour, alpha = point_alpha
        ) +
        scale_shape_manual(values = c(16, 15, 17)))
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
