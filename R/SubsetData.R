#' Create a subset of data
#'
#' This function creates a subset of dataset according to specified criteria.
#'
#' @param data Data, structured in a dataframe
#' @param x Variable(s) used to subset the data
#' @param x.cat Specific criteria (value(s)) of \code{x} used to subset the data
#' @param include Include/exclude value specified in \code{x.cat}. If
#'   \code{include = TRUE} (default value), then observations with
#'   \code{x} = \code{x.cat} are selected. If \code{include = FALSE}, then
#'   observations with value specified in \code{x.cat} will be omitted from the
#'   subset of the data.
#' @details
#' \itemize{
#' \item{If there are several variable \code{x} used as criteria to subset the
#' data, \code{x.cat} can be written as \code{list(... , ...)}}
#' \item{Values in \code{x.cat} should be specified in the same order as the
#' \code{x}'s}
#' }
#' @return This function returns a subset of data
#' @examples
#'
#' # Create data of NM-400 (Multi-walled carbon nanotubes) from geninvitro dataset
#' data.sub<-SubsetData(data=geninvitro, x="name",
#' x.cat="NM-400 (Multi-walled carbon nanotubes)", include = TRUE)
#'
#' # Create data of NM-400 (Multi-walled carbon nanotubes)
#' # with DNA STRAND BREAKS as the endpoint from geninvitro dataset
#' data.sub<-SubsetData(data=geninvitro, x=c("name","endpoint"),
#' x.cat=list("NM-400 (Multi-walled carbon nanotubes)","DNA STRAND BREAKS"),
#' include=TRUE)
#'
#' # Exclude NM-400 (Multi-walled carbon nanotubes) from geninvitro dataset
#' data.sub<-SubsetData(data=geninvitro, x="name",
#' x.cat="NM-400 (Multi-walled carbon nanotubes)", include = FALSE)
#'
#' # Create data of NM-400 (Multi-walled carbon nanotubes)
#' # and NM-110 (Zinc Oxide, uncoated) from geninvitro dataset
#' data.sub<-SubsetData(data=geninvitro, x="name",
#' x.cat=c("NM-400 (Multi-walled carbon nanotubes)",
#' "NM-110 (Zinc Oxide, uncoated)"), include = TRUE)
#'
#' # Create data of NM-400 (Multi-walled carbon nanotubes)
#' # and NM-110 (Zinc Oxide, uncoated), with DNA STRAND BREAKS as the endpoint
#' data.sub<-SubsetData(data=geninvitro, x=c("name","endpoint"),
#' x.cat=list(c("NM-400 (Multi-walled carbon nanotubes)",
#' "NM-110 (Zinc Oxide, uncoated)"),"DNA STRAND BREAKS"),include = TRUE)
#'
#' # Create a new dataset containing only control values from geninvitro dataset
#' controldata<-SubsetData(data=geninvitro, x="name", x.cat=c("control", "Control",
#' "medium", "medium + BSA", "untreated"))
#' @export
SubsetData<-function(data,x,x.cat,include=TRUE){
  n<-length(x)
  if (n==1) {
    if (include==TRUE) {
      data<-subset(data, data[[x]] %in% x.cat)
    } else if (include==FALSE) {
      data<-subset(data, !(data[[x]] %in% x.cat))
    }
  } else {
    for (i in 1:n){
      if (include==TRUE){
        data<-subset(data, data[[x[i]]] %in% unlist(c(x.cat[i])))
      } else if (include==FALSE){
        data<-subset(data, !(data[[x[i]]] %in% unlist(c(x.cat[i]))))
      }
    }
  }
  return(data)
}


