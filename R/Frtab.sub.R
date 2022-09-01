#' Count the number of observations for each unique value of a variable on a
#' subset of data
#'
#' This function counts the number of observations for each unique value of a
#' variable on a subset of data. The subset of data is created according to the
#' specified endpoint and unit of measurement of the dose.
#'
#' @usage Frtab.sub(data.nm, data.control, id, nano, dose, end, end.cat, unit,
#'   unit.cat, control.opt=c("same","all"), x)
#' @param data.nm Data containing the result of toxicity study
#' @param data.control Data of control values
#' @param id Identifier of the experiment
#' @param nano Name of the nanomaterial
#' @param dose Dose or concentration
#' @param end Toxicity endpoint
#' @param end.cat Specific toxicity endpoint of interest
#' @param unit Unit of measurement of the dose
#' @param unit.cat Specific unit of measurement of the dose
#' @param control.opt Option for the control doses. If only control doses with
#'   the same unit of measurement as the non-control ones are included, then
#'   specify "\code{same}" in the \code{control.opt}. If all control doses with
#'   any units of measurement are included, then specify "\code{all}".
#' @param x Variable to be explored
#' @return This function generates a table containing the number of observations
#'   for each unique value of a variable on a subset of data
#' @details \itemize{
#' \item{This function counts for each nanomaterial in the dataset. The
#' different types of nanomaterials are identified by their names. Therefore, if
#' some control values are named differently (see: \code{\link{geninvitro}}
#' dataset and the \code{Examples}), a separate dataset containing only these
#' values first needs to be created. Controls in the new dataset can be linked
#' to the non-control observations belonging to the same experiment through the
#' identifier of the experiment (the linking is performed inside this function).
#' In this situation, it is necessary to have an indicator that can identify
#' different experiments (such as experiment ID).}
#' \item{If all controls in the dataset are named according to the related
#' nanomaterial names, \code{data.control} and \code{id} do not need to be
#' specified.}
#' \item{If doses used in the experiment are all measured in the same unit of
#' measurement, then specify "\code{same}" in \code{control.opt}}.
#' }
#' @examples
#' # Create a dataset containing controls (which are named differently)
#' # from geninvitro dataset:
#' controldata<-SubsetData(data=geninvitro, x="name", x.cat=c("control",
#'              "Control", "medium", "medium + BSA", "untreated"))
#'
#' # Exclude controls (which are named differently) from geninvitro dataset:
#' invitrodata<-SubsetData(data=geninvitro, x="name", x.cat=c("control",
#'              "Control", "medium", "medium + BSA", "untreated"), include=FALSE)
#'
#' # Frequency of each unique value of the dose for each nanomaterial
#' # in geninvitro dataset, with DNA STRAND BREAKS as the toxicity endpoint
#' # and only observations with concentration measured in "ug/cm2" are considered:
#' #
#' Frtab.sub(data.nm=invitrodata, data.control=controldata, nano="name",
#'            end="endpoint", end.cat="DNA STRAND BREAKS", id="experimentID",
#'            dose="concentration", unit="concentration_unit", unit.cat="ug/cm2",
#'            control.opt="same", x="concentration")
#'
#' @export
Frtab.sub<-function(data.nm, data.control, id, nano, dose,
                     end, end.cat, unit, unit.cat,
                     control.opt=c("same","all"), x){
  control.opt<-match.arg(control.opt)
  fr.nano<-vector(mode = "list")
  for (nano.cat in unique(data.nm[[nano]])) {
    df<-data.nm[data.nm[[nano]]==nano.cat &
                  data.nm[[end]]==end.cat &
                  !is.na(data.nm[[nano]]) &
                  !is.na(data.nm[[end]]) &
                  !is.na(data.nm[[dose]]),]
    df.un<-df[df[[unit]]==unit.cat
              & !is.na(df[[unit]]),]
    df.sub1<-df[as.numeric(as.character(df[[dose]]))!=0.0 &
                  df[[unit]]==unit.cat
                & !is.na(df[[unit]]),]
    df.sub2<-df[as.numeric(as.character(df[[dose]]))==0.0,]
    df.com<-rbind(df.sub1,df.sub2)
    if (missing(data.control)) {
      df.same<-df.un
      df.all<-df.com
    } else {
      c.end<-data.control[data.control[[end]]==end.cat &
                            !is.na(data.control[[end]]) &
                            !is.na(data.control[[dose]]),]
      c.exp<-c.end[c.end[[id]] %in%  df[[id]],]
      c.un<-c.exp[c.exp[[unit]]==unit.cat
                  & !is.na(c.exp[[unit]]),]
      df.same<-rbind(df.un,c.un)
      df.all<-rbind(df.com,c.exp)
    }
    df.nm<-switch(control.opt,same=df.same,all=df.all)
    if (nrow(df.nm)==0) next
    fr.nano[[nano.cat]]<-table(droplevels(as.factor(as.character(df.nm[[x]]))),useNA = "ifany")
  }
  return(fr.nano)
}
