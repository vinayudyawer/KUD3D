#' Calculate 3DKUD volume
#'
#' @description Helper function to calculate 3DKUD volume for a `kde` object
#'
#' @param fhat a kernel density estimate object of class `kde`. see \code{\link{?kde}} using the `ks` package
#' @param cont voxel contour for which volume will be calulated (e.g. 50, 95)
#'
#' @return A numeric vector with volume of kde object
#'
#' @import ks
#' @importFrom ks contourLevels
#'
#' @export vol3d
#'
#' @examples
#' ## Import detection data
#' data(exampletag)
#'
#' ## Estimate 3D kernel density using the ks package
#' library(ks)
#' H.pi <- Hpi(exampletag, binned = TRUE) * 3
#' fhat <- kde(exampletag, H = H.pi)
#'
#' ## Calculate 3D KUD volume
#' vol3d(fhat = fhat, cont = 50)   ## volume of 50% contour in m3
#' vol3d(fhat = fhat, cont = 95)   ## volume of 95% contour in m3
#'
#'

vol3d <- function(fhat, cont = NULL) {
  if (!inherits(fhat, "kde"))
    stop(
      "Oops! Input fhat data needs to be a 'kde' object.\nUse the kde() function before running this operation. see ?kde"
    )

  if (is.null(cont))
    stop("You havent supplied a contour level. Provide a single value.")

  if (length(cont) != 1)
    stop("Whoah! one at a time please. I can only calculate volume for one contour at a time.")

  ct <- contourLevels(fhat, cont = cont, approx = TRUE)

  vol.voxel <- prod(sapply(fhat$eval.points, diff)[1, ])

  no.voxel <- sum(fhat$estimate > ct)

  out <- no.voxel * vol.voxel

  return(out)

}

