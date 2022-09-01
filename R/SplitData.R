#' Split the data of each nanomaterial into different subsets of data
#'
#' This function splits the data of each nanomaterial into different subsets of
#' data according to the unique values of selected variable(s)
#'
#' @usage SplitData(data.nm, data.control, id, nano, dose, end, end.cat, unit,
#'   unit.cat, control.opt=c("same","all"), vars)
#' @param data.nm Data containing the result of toxicity study
#' @param data.control Data of control values
#' @param id Identifier of the experiment
#' @param nano Name of the nanomaterial
#' @param dose Dose or concentration
#' @param end Toxicity endpoint
#' @param end.cat Specific toxicity endpoint of interest
#' @param unit Unit of measurement of the dose
#' @param unit.cat Specific unit of measurement of the dose
#' @param control.opt Option for the control doses if \code{unit} and
#'   \code{unit.cat} are specified. If only control doses with the same unit of
#'   measurement as the non-control ones are included, then specify
#'   "\code{same}" in the \code{control.opt}. If all control doses with any
#'   units of measurement are included, then specify "\code{all}".
#' @param vars Variables used to split the data
#' @return This function splits the data of each nanomaterial into different
#'   subsets of data according to the unique values of selected variable(s)
#' @examples
#' # Example 1:
#' # Create a dataset containing controls (which are named differently)
#' # from geninvitro dataset:
#' controldata<-SubsetData(data=geninvitro, x="name", x.cat=c("control", "Control",
#'              "medium", "medium + BSA", "untreated"))
#'
#' # Exclude controls (which are named differently) from geninvitro dataset:
#' invitrodata<-SubsetData(data=geninvitro, x="name", x.cat=c("control", "Control",
#'              "medium", "medium + BSA", "untreated"), include=FALSE)
#' #
#' # Split geninvitro data according to the cell type, method, study provider,
#' # unit of the concentration and the type of the endpoint:
#' datasub<-SplitData(data.nm=invitrodata, data.control=controldata,
#'                    id="experimentID", nano="name", dose="concentration",
#'                    vars=c("celltype", "method","studyprovider",
#'                    "concentration_unit","endpoint"))
#'
#' # Example 2:
#' # Split geninvitro data with DNA STRAND BREAKS as the endpoint, according
#' # to the cell type, method, study provider, and unit of the concentration:
#' datasub2<-SplitData(data.nm=invitrodata, data.control=controldata,
#'                     id="experimentID", nano="name", dose="concentration",
#'                     end="endpoint", end.cat="DNA STRAND BREAKS",
#'                     vars=c("celltype","method","studyprovider",
#'                     "concentration_unit"))
#'
#' # Example 3:
#' # Split geninvitro data with DNA STRAND BREAKS as the endpoint and
#' # concentration measured in ug/cm2, according to the cell type:
#' datasub3<-SplitData(data.nm=invitrodata, data.control=controldata,
#'                     id="experimentID", nano="name", dose="concentration",
#'                     end="endpoint", end.cat="DNA STRAND BREAKS",
#'                     unit="concentration_unit", unit.cat="ug/cm2",
#'                     control.opt="same", vars="celltype")
#' @importFrom forcats fct_explicit_na
#' @export
SplitData<-function(data.nm, data.control, id, nano, dose,
                        end, end.cat, unit, unit.cat,
                        control.opt=c("same","all"), vars){
  control.opt<-match.arg(control.opt)
  df.nm<-vector(mode="list")
  subdata<-vector(mode="list")
  for (nano.cat in unique(data.nm[[nano]])) {
    if (missing(end)==FALSE && missing(end.cat)==FALSE) {
      df<-data.nm[data.nm[[nano]]==nano.cat &
                    data.nm[[end]]==end.cat &
                    !is.na(data.nm[[nano]]) &
                    !is.na(data.nm[[end]]) &
                    !is.na(data.nm[[dose]]),]

    }
    else
    {
      df<-data.nm[data.nm[[nano]]==nano.cat &
                    !is.na(data.nm[[nano]]) &
                    !is.na(data.nm[[dose]]),]
    }
    if (missing(unit)==FALSE) {
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
        if (missing(end)==FALSE && missing(end.cat)==FALSE) {
          c.end<-data.control[data.control[[end]]==end.cat &
                                !is.na(data.control[[end]]) &
                                !is.na(data.control[[dose]]),]
          c.exp<-c.end[c.end[[id]] %in%  df[[id]],]
          c.un<-c.exp[c.exp[[unit]]==unit.cat
                      & !is.na(c.exp[[unit]]),]
          df.same<-rbind(df.un,c.un)
          df.all<-rbind(df.com,c.exp)
        } else {
          c.exp<-data.control[data.control[[id]] %in%  df[[id]] &
                                !is.na(data.control[[dose]]),]
          #additional here 1: w data.control, w/o end, w unit
          c.un<-c.exp[c.exp[[unit]]==unit.cat
                      & !is.na(c.exp[[unit]]),]
          df.same<-rbind(df.un,c.un)
          df.all<-rbind(df.com,c.exp)
        }
      }
      df.nm[[nano.cat]]<-switch(control.opt,same=df.same,all=df.all)
    } else {
      if (missing(data.control)) {
        df.nm[[nano.cat]]<-df
      } else {
        if (missing(end)==FALSE && missing(end.cat)==FALSE) {
          c.end<-data.control[data.control[[end]]==end.cat &
                                !is.na(data.control[[end]]) &
                                !is.na(data.control[[dose]]),]
          c.exp<-c.end[c.end[[id]] %in%  df[[id]],]
        } else {
          c.exp<-data.control[data.control[[id]] %in%  df[[id]] &
                                !is.na(data.control[[dose]]),]
        }
        df.nm[[nano.cat]]<-rbind(df,c.exp)
      }
    }
    if (nrow(df.nm[[nano.cat]])==0) next
    df.nm[[nano.cat]]$nanomaterial<-nano.cat
    if (missing(vars)==TRUE) {
      subdata[[nano.cat]]<-df.nm[[nano.cat]]
    } else if (missing(vars)==FALSE) {
      if (length(vars)==1) {
        subdata[[nano.cat]]<-split(df.nm[[nano.cat]],
                                   fct_explicit_na(as.character(df.nm[[nano.cat]][[vars]]), "NA"),drop=T)
        # lapply(df.nm[[nano.cat]][vars], addNA),drop=T)
      } else {
        subdata[[nano.cat]]<-split(df.nm[[nano.cat]],
                                   lapply(df.nm[[nano.cat]][,vars], addNA),drop=T)
      }
    }
  }
  return(subdata)
}

