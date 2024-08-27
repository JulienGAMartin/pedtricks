#' Produce a graphical representation of the relatedness matrix of a pedigree
#'
#' Creates the object needed to plot a pedigree's numerator relatedness matrix
#' given a few different choices for ordering. The resulting image for a pedigree
#' of size n can be visualized as a n x n grid of colored squares based on
#' values of the numerator relatedness matrix.
#'
#' @param pedigree A data.frame of a pedigree with 3 columns: id, dam, sire
#' @param order A character expression for how the pedigree should be ordered
#'   for visualization. See Details.
#' @param grp A character expression for the column name in pedigree indicating
#'   how to order the pedigree for visualization.
#' @param \dots Additional plotting arguments passed to \code{\link[Matrix]{image}}.
#'
#' @return A list of class \dQuote{trellis}.
#' @author \email{mew0099@@auburn.edu}
#' @examples
#'
#'  data(gryphons)
#'  pedigree <- fix_ped(gryphons[, 1:3])
#' 
#'   ## draw the gryphon pedigree
#'   draw_pedA(pedigree, order = "original")
#'
#'   ## draw the gryphon pedigree by function assigned generation
#'   (Agen <- draw_pedA(pedigree, order = "generation"))
#'
#'   ## draw the gryphon pedigree by cohort in the dataset
#'     ## add cohort back from original data
#'     pedigree$cohort <- NA
#'     pedigree$cohort[match(gryphons$id, pedigree[, 1])] <- gryphons$cohort
#'   (Achrt <- draw_pedA(pedigree, order = "user", grp = "cohort"))
#' \dontrun{
#'   ## show two images of the same pedigree in different orders
#'   ### (i.e., plotting multiple trellis objects in the same figure)
#'   plot(Agen, position = c(xmin = 0, ymin = 0, xmax = 0.45, ymax = 1),
#'     more = TRUE)
#'   plot(Achrt, position = c(xmin = 0.55, ymin = 0, xmax = 1, ymax = 1))
#' }
#' 
#' @keywords plot
#' @export
draw_pedA <- function(pedigree,
  order = c("original", "generation", "user"), grp = NULL, ...){
  
  ord <- match.arg(order)
  if(ord == "original"){
    ped <- pedigree
    #TODO make axis 1 and 2 labels (e.g., "pedigree row")
  }
    
  if(ord == "generation"){
    pedigree$gen <- genAssign(pedigree[, 1:3])
    ped <- pedigree[order(pedigree$gen), ]
    #TODO make axis 1 and 2 labels (e.g., "cohort/generation")
  }

  if(ord == "user"){
    if(is.null(grp)){
      warning(cat("No value in 'grp' argument when order='user'.",
        "\n Using the original order instead\n"))
      ped <- pedigree
    } else{
        if(grp %in% names(pedigree)){
          ped <- pedigree[order(pedigree[, grp]), ]
          #TODO make axis 1 and 2 labels (i.e., grp)
        } else{
            warning(cat("No column named", grp, "in 'pedigree'",
              "\n Using the original order instead\n")) 
            ped <- pedigree  
          }
      }  #<-- end if/else grp=NULL
  }  #<-- end if ord=user
      
  ##
  A <- makeA(ped[, 1:3])
  Matrix::image(A, colorkey = TRUE, ...)
}

