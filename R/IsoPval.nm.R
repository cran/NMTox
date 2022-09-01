#' Calculate p-values based on permutation for both increasing and
#' decreasing ordered alternatives
#'
#' This function calculates p-values based on permutation. The plots of
#' the null distribution and the observed test statistic under increasing and
#' decreasing ordered alternatives are also given.
#'
#' @usage IsoPval.nm(data.nm, dose, response, stat=c("E2", "Williams","Marcus",
#'   "M", "ModM"), niter, nano.cat=NULL)
#' @param data.nm Nanomaterial dataset
#' @param dose Dose or concentration (with the same unit of measurement)
#' @param response Response (a certain endpoint value)
#' @param stat Test statistics ("\code{E2}" for the global likelihood test,
#'   "\code{Williams}" for Williams test, "\code{Marcus}" for Marcus test,
#'   "\code{M}" for M test or "\code{ModM}" for modified M test)
#' @param niter Number of permutations
#' @param nano.cat Name of the nanomaterial
#' @return This value provides p-values based on the permutation, the plot of
#'   the null distribution and the observed test statistics.
#' @details This function is intended to be used inside the function
#'   \code{\link{Isotest}}. However, it can also be used to calculate p-values
#'   and generate the plot for one nanomaterial, with a particular unit of
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
#' \donttest{IsoPval.nm(data.nm=nm400, dose="concentration", response="value",
#'            stat="E2", niter=1000)}
#' @importFrom stats density
#' @importFrom grDevices dev.new
#' @import graphics
#' @export
IsoPval.nm<-function(data.nm, dose, response,
                         stat=c("E2", "Williams","Marcus","M","ModM"),
                         niter, nano.cat=NULL){
  stat<-match.arg(stat)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  x <- as.numeric(as.character(data.nm[[dose]]))
  y <- as.numeric(as.character(data.nm[[response]]))
  obs<-IsoTest.nm(data.nm=data.nm,dose=dose,response=response, stat=stat)
  obs.u<-as.numeric(obs[1])
  obs.d<-as.numeric(obs[2])
  exp.u <- exp.d <- 1:niter
  x.niter <- t(sapply(1:niter, function(i) sample(x)))
  for (j in 1:niter) {
    exps <- IsoTest.nm(dose=x.niter[j, ], response=y, stat=stat)
    exp.u[j] <- as.numeric(exps[1])
    exp.d[j] <- as.numeric(exps[2])
  }
  rawp.u <- sum(obs.u < exp.u)/niter
  rawp.d <- sum(obs.d > exp.d)/niter
  if (stat == "E2") {
    rawp.d <- sum(obs.d < exp.d)/niter
  }
  par(mfrow = c(2, 1))
  hist(exp.u, main = "", nclass = 1000, col = 0, probability = TRUE,
       xlim = c(min(exp.u, obs.u), max(exp.u, obs.u)), xlab = paste(stat))
  dx <- density(exp.u, from = min(exp.u), to = max(exp.u))
  lines(dx$x, dx$y, lwd = 3, col = 5)
  abline(v = obs.u, col = 7, lwd = 3)
  title(paste(nano.cat,"\np-value=",rawp.u,"(up)"))
  hist(exp.d, main = "", nclass = 1000, col = 0, probability = TRUE,
       xlim = c(min(exp.d, obs.d), max(exp.d, obs.d)), xlab = paste(stat))
  dx <- density(exp.d, from = min(exp.d), to = max(exp.d))
  lines(dx$x, dx$y, lwd = 3, col = 5)
  abline(v = obs.d, col = 7, lwd = 3)
  title(paste(nano.cat,"\np-value=", rawp.d, "(down)"))
  rawp<-data.frame("pvalue.up"=rawp.u,"pvalue.down"=rawp.d)
  return(rawp)
}

