#' Calculate the value of the test statistics for testing monotonic trend
#'
#' This function provides the value of the test statistics (the global
#' likelihood test, Williams, Marcus, M or modified M test) for one nanomaterial
#'
#' @usage IsoTest.nm(data.nm, dose, response, stat=c("E2", "Williams","Marcus",
#'   "M", "ModM"))
#' @param data.nm Nanomaterial dataset
#' @param dose Dose or concentration (with the same unit of measurement)
#' @param response Response (a certain endpoint value)
#' @param stat Test statistics ("\code{E2}" for the global likelihood test,
#'   "\code{Williams}" for Williams test, "\code{Marcus}" for Marcus test,
#'   "\code{M}" for M test or "\code{ModM}" for modified M test)
#' @return This function calculates the value of the specified test statistics
#'   for one nanomaterial
#' @details This function is intended to be used inside the function
#'   \code{\link{Isotest}}. However, it can also be used to obtain the value of
#'   a specified test statistics for one nanomaterial, with a particular unit of
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
#' IsoTest.nm(data.nm=nm400, dose="concentration", response="value", stat="E2")
#' @importFrom Iso pava
#' @export
IsoTest.nm<-function(data.nm, dose, response,
                      stat=c("E2", "Williams", "Marcus","M","ModM")){
  stat<-match.arg(stat)
  if (missing(data.nm)) {
    x <- as.numeric(dose)
    y <- as.numeric(response)
  } else {
    x <- as.numeric(as.character(data.nm[[dose]]))
    y <- as.numeric(as.character(data.nm[[response]]))
  }
  ordx <- order(x)
  x <- x[ordx]
  y <- y[ordx]
  k<-length(unique(x))
  rep.x<-table(x)
  n0<-rep.x[1]
  nk<-rep.x[k]
  n<-sum(rep.x)
  ybar.x<-tapply(y, as.factor(x), mean)
  ybar.x.all<-rep(ybar.x, rep.x)
  ybar.0<-ybar.x[1]
  yiso.u<-pava(y = ybar.x,w = rep.x)
  yiso.u.0<-yiso.u[1]
  yiso.u.all<-rep(yiso.u, rep.x)
  yiso.d<-pava(y = ybar.x,w = rep.x,decreasing = TRUE)
  yiso.d.all<-rep(yiso.d, rep.x)
  yiso.d.0<-yiso.d[1]
  s2<-(sum((y-ybar.x.all)^2))/(n-k)
  tw.u<-(yiso.u[k]-ybar.0)/(sqrt(s2*(1/nk+1/n0)))
  tw.d<-(yiso.d[k]-ybar.0)/(sqrt(s2*(1/nk+1/n0)))
  tm.u<-(yiso.u[k]-yiso.u.0)/(sqrt(s2*(1/nk+1/n0)))
  tm.d<-(yiso.d[k]-yiso.d.0)/(sqrt(s2*(1/nk+1/n0)))
  rss1.u<-sum((y-yiso.u.all)^2)
  rss1.d<-sum((y-yiso.d.all)^2)
  tM.u <-(yiso.u[k]-yiso.u.0)/sqrt(rss1.u/(n-k))
  tM.d <-(yiso.d[k]-yiso.d.0)/sqrt(rss1.d/(n-k))
  tmodM.u <-(yiso.u[k]-yiso.u.0)/sqrt(rss1.u/(n-length(unique(yiso.u))))
  tmodM.d <-(yiso.d[k]-yiso.d.0)/sqrt(rss1.d/(n-length(unique(yiso.d))))
  rss0<-sum((y-mean(y))^2)
  E.u<-1-(rss1.u/rss0)
  E.d<-1-(rss1.d/rss0)
  dir <- if (rss1.u <= rss1.d) {
    "u"} else {
      "d"}
  res.E2<-data.frame("test.stat(up)"=E.u,"test.stat(down)"=E.d,"direction"=dir)
  res.Williams<-data.frame("test.stat(up)"=as.numeric(tw.u),"test.stat(down)"=as.numeric(tw.d),"direction"=dir)
  res.Marcus<-data.frame("test.stat(up)"=as.numeric(tm.u),"test.stat(down)"=as.numeric(tm.d),"direction"=dir)
  res.M<-data.frame("test.stat(up)"=as.numeric(tM.u),"test.stat(down)"=as.numeric(tM.d),"direction"=dir)
  res.modM<-data.frame("test.stat(up)"=as.numeric(tmodM.u),"test.stat(down)"=as.numeric(tmodM.d),"direction"=dir)
  res<-switch(stat,E2=res.E2,Williams=res.Williams,Marcus=res.Marcus, M=res.M, ModM=res.modM)
  return(res)
}
