#' Create a plot of the observations, sample means and fitted isotonic
#' regression curve for one nanomaterial
#'
#' This function generates a dose-response plot (scatter plot) of the
#' observations, sample means and fitted isotonic regression curve for one
#' nanomaterial
#' @usage IsoPlot.nm(data.nm, dose, response, type = c("continuous", "ordinal"),
#'   add.curve = TRUE, nano.cat = NULL, xlabel="Dose", ylabel="Response")
#' @param data.nm Dataset of a particular nanomaterial
#' @param dose Dose or concentration (with the same unit of measurement)
#' @param response Response (a certain endpoint value)
#' @param type Type of the dose
#' @param add.curve Adding curve to the plot
#' @param nano.cat Title of the plot (referring to the name of the nanomaterial)
#' @param xlabel label for the x-axis of the plot
#' @param ylabel label for the y-axis of the plot
#' @details This function is intended to be used inside the function
#'   \code{\link{Isoplot}}. However, it can also be used to generate a plot
#'   for one nanomaterial, with a particular unit of measurement of the dose and
#'   for a certain toxicity endpoint.
#' @return This function produces a plot of the observations, sample means and
#'   fitted isotonic regression curve for one nanomaterial
#' @references Lin D., Pramana, S., Verbeke, T., and Shkedy, Z. (2015). IsoGene:
#'   Order-Restricted Inference for Microarray Experiments. R package version
#'   1.0-24. \url{https://CRAN.R-project.org/package=IsoGene}
#'
#'   Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, L. (editors).
#'   (2012) Modeling Doseresponse Microarray Data in Early Drug Development
#'   Experiments Using R. Springer.
#' @examples
#' #nm400 contains the result of genetic toxicity in vitro study of NM-400
#' #(Multi-walled carbon nanotubes) with associated controls
#' IsoPlot.nm(data.nm=nm400, dose="concentration", response="value",
#'            nano.cat="Multi-walled carbon nanotubes", xlabel="Concentration",
#'            ylabel="DNA STRAND BREAKS")
#' @import ggplot2
#' @importFrom Iso pava
#' @export
IsoPlot.nm<-function(data.nm, dose, response, type = c("continuous",
                    "ordinal"), add.curve = TRUE, nano.cat = NULL,
                    xlabel="Dose", ylabel="Response")
{
  type<-match.arg(type)
  x <- as.numeric(as.character(data.nm[[dose]]))
  y <- as.numeric(as.character(data.nm[[response]]))
  type <- match.arg(type)
  if (!(type %in% c("continuous", "ordinal")))
    stop("The dose can be only continuous or ordinal")
  miny <- min(y)
  maxy <- max(y)
  ordx <- order(x)
  unx <- sort(unique(x))
  y1 <- as.numeric(y)[ordx]
  x1 <- x[ordx]
  ybar.x <- tapply(y1, as.factor(x1), mean)
  ybar.tot <- rep(mean(y), length(unx))
  rep.x <- table(x)
  k <- length(rep.x)
  yiso.u <- pava(y = ybar.x, w = rep.x)
  yiso.d <- pava(y = ybar.x, w = rep.x, decreasing = TRUE)
  dir <- IsoTest.nm(data.nm = data.nm, dose = dose, response = response,
                    stat = "E2")[3]
  sort.x <- sort(x)
  unx.fac <- as.factor(unx)
  dt1 <- data.frame(cbind(sort.x, y1))
  dt2 <- data.frame(cbind(unx, ybar.x, yiso.u, yiso.d))
  if (type == "continuous") {
    if (dir == "u") {
      p <- ggplot() + geom_point(data = dt1, aes(x = sort.x,
                                                 y = y1), size = 2) + geom_point(data = dt2, aes(x = unx,
                                                                                                 y = ybar.x), shape = 24, size = 2, col = 2, fill = "red") +
        geom_point(data = dt2, aes(x = unx, y = yiso.u),
                   shape = 4, size = 2) + xlab(xlabel) +
        ylab(ylabel) + ggtitle(nano.cat)
      if (add.curve) {
        plot.iso <- p + geom_line(data = dt2, aes(x = unx,
                                                  y = yiso.u), col = 4, lwd = 1) + theme_bw()
      }
      else {
        plot.iso <- p + theme_bw()
      }
    }
    else {
      p <- ggplot() + geom_point(data = dt1, aes(x = sort.x,
                                                 y = y1), size = 2) + geom_point(data = dt2, aes(x = unx,
                                                                                                 y = ybar.x), shape = 24, size = 2, col = 2, fill = "red") +
        geom_point(data = dt2, aes(x = unx, y = yiso.d),
                   shape = 4, size = 2) + xlab(xlabel) +
        ylab(ylabel) + ggtitle(nano.cat)
      if (add.curve) {
        plot.iso <- p + geom_line(data = dt2, aes(x = unx,
                                                  y = yiso.d), col = 4, lwd = 1) + theme_bw()
      }
      else {
        plot.iso <- p + theme_bw()
      }
    }
  }
  else {
    if (dir == "u") {
      p <- ggplot() + geom_point(data = dt1, aes(x = as.factor(sort.x),
                                                 y = y1), size = 2) + geom_point(data = dt2, aes(x = as.factor(unx),
                                                                                                 y = ybar.x), shape = 24, size = 2, col = 2, fill = "red") +
        geom_point(data = dt2, aes(x = as.factor(unx),
                                   y = yiso.u), shape = 4, size = 2) + xlab(xlabel) +
        ylab(ylabel) + ggtitle(nano.cat)
      if (add.curve) {
        plot.iso <- p + geom_line(data = dt2, aes(x = as.factor(unx),
                                                  y = yiso.u, group = 1), col = 4, lwd = 1) + theme_bw()
      }
      else {
        plot.iso <- p + theme_bw()
      }
    }
    else {
      p <- ggplot() + geom_point(data = dt1, aes(x = as.factor(sort.x),
                                                 y = y1), size = 2) + geom_point(data = dt2, aes(x = as.factor(unx),
                                                                                                 y = ybar.x), shape = 24, size = 2, col = 2, fill = "red") +
        geom_point(data = dt2, aes(x = as.factor(unx),
                                   y = yiso.d), shape = 4, size = 2) + xlab(xlabel) +
        ylab(ylabel) + ggtitle(nano.cat)
      if (add.curve) {
        plot.iso <- p + geom_line(data = dt2, aes(x = as.factor(unx),
                                                  y = yiso.d, group = 1), col = 4, lwd = 1) + theme_bw()
      }
      else {
        plot.iso <- p + theme_bw()
      }
    }
  }
  return(plot.iso)
}
