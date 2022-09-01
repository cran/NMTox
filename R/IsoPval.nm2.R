#' Calculate p-values based on permutation for both increasing and
#' decreasing ordered alternatives
#'
#' This function calculates p-values based on permutation.
#'
#' @usage IsoPval.nm2(data.nm, dose, response, stat=c("E2", "Williams","Marcus",
#'   "M","ModM"), niter, nano.cat=NULL)
#' @param data.nm Nanomaterial dataset
#' @param dose Dose or concentration (with the same unit of measurement)
#' @param response Response (a certain endpoint value)
#' @param stat Test statistics ("\code{E2}" for the global likelihood test,
#'   "\code{Williams}" for Williams test, "\code{Marcus}" for Marcus test,
#'   "\code{M}" for M test, or "\code{ModM} for modified M test")
#' @param niter Number of permutations
#' @param nano.cat Name of the nanomaterial
#' @return This function provides p-values based on the permutation
#' @details This function is intended to be used inside the function
#'   \code{\link{Isotest}}. However, it can also be used to calculate p-values
#'   for one nanomaterial, with a particular unit of
#'   measurement of the dose, for a certain toxicity endpoint.
#' @references Lin D., Pramana, S., Verbeke, T., and Shkedy, Z. (2015). IsoGene:
#'   Order-Restricted Inference for Microarray Experiments. R package version
#'   1.0-24. \url{https://CRAN.R-project.org/package=IsoGene}
#'
#'   Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, L. (editors).
#'   (2012) Modeling Doseresponse Microarray Data in Early Drug Development
#'   Experiments Using R. Springer.
#' @seealso \code{\link{Isotest}}
#' @examples
#' #nm400 contains the result of genetic toxicity in vitro study of NM-400
#' #(Multi-walled carbon nanotubes) with associated controls
#' IsoPval.nm2(data.nm=nm400, dose="concentration", response="value",
#'             stat="E2", niter=1000)
#' @importFrom stats density
#' @importFrom grDevices dev.new
#' @import graphics
#' @noRd
IsoPval.nm2=function(data.nm, dose, response,
                     stat = c("E2", "Williams","Marcus","M","ModM"),
                     niter, nano.cat = NULL)
{
  stat<-match.arg(stat)
  x <- as.numeric(as.character(data.nm[[dose]]))
  y <- as.numeric(as.character(data.nm[[response]]))
  obs <- IsoTest.nm(data.nm = data.nm, dose = dose, response = response,
                    stat = stat)
  obs.u <- as.numeric(obs[1])
  obs.d <- as.numeric(obs[2])
  exp.u <- exp.d <- 1:niter
  x.niter <- t(sapply(1:niter, function(i) sample(x)))
  for (j in 1:niter) {
    exps <- IsoTest.nm(dose = x.niter[j, ], response = y,
                       stat = stat)
    exp.u[j] <- as.numeric(exps[1])
    exp.d[j] <- as.numeric(exps[2])
  }
  rawp.u <- sum(obs.u < exp.u)/niter
  rawp.d <- sum(obs.d > exp.d)/niter
  if (stat == "E2") {
    rawp.d <- sum(obs.d < exp.d)/niter
  }
  rawp <- data.frame(pvalue.up = rawp.u, pvalue.down = rawp.d)
  return(rawp)
}
