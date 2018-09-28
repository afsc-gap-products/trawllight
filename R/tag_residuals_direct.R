#' Direct method for detecting tag obstruction
#'
#' \code{tag_residuals_direct} uses a linear model to estimate the relationship between surface light and light below the surface for the specified depth bins.
#'
#' @param x Data frame containing surface light, light at depth, depth, and unique identifiers for survey vessels, years, etc.
#' @param formula an object of class formula, passed to \code{lm}
#' @param water.col A character vector with the name of the column containing water column light measurements.
#' @param surface.col A character vector with the name of the column containing surface light measurements.
#' @param depth.col A character vector with the name of the column containing depths.
#' @param depth.bins Depth bins for which linear regression should be conducted.
#' @param ... Additional arguments
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}

tag_residuals_direct <- function(x, formula = log10(trans_llight) ~ log10(surf_trans_llight) + interaction(vessel, cruise), water.col = "trans_llight", surface.col = "surf_trans_llight", depth.col = "cdepth", depth.bins = c(1, 3, 5, 7, 9), ...) {
  names(x)[names(x) == water.col] <- "trans_llight"
  names(x)[names(x) == surface.col] <- "surf_trans_llight"
  names(x)[names(x) == depth.col] <- "cdepth"
  lout <- list()

  if(mean(c(depth.bins) %in% x$cdepth) < 1) {
    stop(paste0("tag_residuals_direct: Cannot calculate residuals. Some depth.bins not found in ", depth.col))
  }

  for(i in 1:length(depth.bins)) {
    x_sub <- subset(x, cdepth == depth.bins[i])

    # Linear regression
    DIRECT_LM <- lm(formula, data = x_sub, ...)

    # Append lin. reg. residuals to subset data frame
    x_sub$direct_residual <- residuals(DIRECT_LM)
    if(depth.bins[1] == depth.bins[i]) {
      output.df <- x_sub
    } else {
      output.df <- plyr::rbind.fill(output.df, x_sub, ...)
    }
   lout[[i]] <- DIRECT_LM
   names(lout)[i] <- paste0("lm_dbin", depth.bins[i])

    names(x)[names(x) == "trans_llight"] <- water.col
    names(x)[names(x) == "surf_trans_llight"] <- surface.col
    names(x)[names(x) == "cdepth"] <- depth.col
  }

  lout$resid_df <- output.df
  return(lout)
}
