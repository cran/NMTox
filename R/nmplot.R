#' Create plots of the dose and the response for each nanomaterial in the
#' dataset
#'
#' This function generates scatter plots of the dose and the response for each
#' nanomaterial in the dataset
#'
#' @usage nmplot(data.nm, data.control, id, nano, response, dose, end, end.cat,
#'   unit, unit.cat, type=c("dose","log"), control.opt=c("same","all"), vars,
#'   nrow=1, ncol=1)
#' @param data.nm Data containing the result of toxicity study
#' @param data.control Data of control values
#' @param id Identifier of the experiment
#' @param nano Name of the nanomaterial
#' @param response Response (endpoint value)
#' @param dose Dose or concentration
#' @param end Toxicity endpoint
#' @param end.cat Specific toxicity endpoint of interest
#' @param unit Unit of measurement of the dose
#' @param unit.cat Specific unit of measurement of the dose
#' @param type Type of the dose to be plotted ("\code{dose}" for dose and
#'   "\code{log}" for log(dose))
#' @param control.opt Option for the control doses if \code{unit} and
#'   \code{unit.cat} are specified. If only control doses with the same unit of
#'   measurement as the non-control ones are included, then specify
#'   "\code{same}" in the \code{control.opt}. If all control doses with any
#'   units of measurement are included, then specify "\code{all}".
#' @param vars Variables used to subset the data
#' @param nrow Number of row in the plotting space (default is 1)
#' @param ncol Number of column in the plotting space (default is 1)
#' @return This function produces plots of the dose and the response for each
#'   nanomaterial in the dataset
#' @details
#' \itemize{
#' \item{This function generates plots for each nanomaterial in the dataset (or
#' for each subset of data). The different types of nanomaterials are identified
#' by their names. Therefore, if some control values are named differently (see:
#' \code{\link{geninvitro}} dataset and the \code{Examples}), a separate dataset
#' containing only these values first needs to be created. Controls in the new
#' dataset can be linked to the non-control observations belonging to the same
#' experiment through the identifier of the experiment (the linking is performed
#' inside this function). In this situation, it is necessary to have an
#' indicator that can identify different experiments (such as experiment ID).}
#' \item{If all controls in the dataset are named according to the related
#' nanomaterial names, \code{data.control} and \code{id} do not need to be
#' specified.}
#' \item{If doses used in the experiment are all measured in the same unit of
#' measurement, then specify "\code{same}" in \code{control.opt}}.
#' \item{Dose-response plot can also be generated for subsets of data in each
#' nanomaterial, by specifying the variables used to split the data in
#' \code{vars}}.
#' }
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
#' # Generate dose-response plot for geninvitro, with DNA STRAND BREAKS as the
#' # endpoint, concentrations measured in ug/cm2 and controls measured in any
#' # units of measurement:
#' #
#' \donttest{nmplot(data.nm=invitrodata, data.control=controldata, id="experimentID",
#'        nano="name", response="value", dose="concentration", end="endpoint",
#'        end.cat="DNA STRAND BREAKS", unit="concentration_unit",
#'        unit.cat="ug/cm2", type="dose", control.opt="all")}
#'
#' # Example 2:
#' # Split geninvitro data according to the cell type and unit of the
#' # concentration and generate dose-response plot for each
#' # subset of data with DNA STRAND BREAKS as the endpoint:
#' #
#' \donttest{nmplot(data.nm=invitrodata, data.control=controldata, id="experimentID",
#'        nano="name", response="value", dose="concentration", end="endpoint",
#'        end.cat="DNA STRAND BREAKS", type="dose", nrow=2, ncol=2,
#'        vars=c("celltype","concentration_unit"))}
#'
#' @import ggplot2
#' @importFrom gridExtra marrangeGrob
#' @export
nmplot<-function(data.nm, data.control, id, nano, response, dose,
                     end, end.cat, unit, unit.cat, type=c("dose","log"),
                     control.opt=c("same","all"), vars, nrow=1, ncol=1){
  if (missing(nano)==TRUE) {stop("argument \"nano\" is missing")
  } else if (missing(response)==TRUE) {stop("argument \"response\" is missing")
  } else if (missing(dose)==TRUE) {stop("argument \"dose\" is missing")
  } else if (missing(unit)==TRUE && missing(unit.cat)==FALSE) {stop("argument \"unit\" is missing")
  } else
  {
    type<-match.arg(type)
    control.opt<-match.arg(control.opt)
    dose.plot<-vector(mode="list")
    dplot<-vector(mode="list")
    for (nano.cat in unique(data.nm[[nano]])) {
      if (missing(end)==FALSE && missing(end.cat)==FALSE) {
        df<-data.nm[data.nm[[nano]]==nano.cat &
                      data.nm[[end]]==end.cat &
                      !is.na(data.nm[[nano]]) &
                      !is.na(data.nm[[end]]) &
                      !is.na(data.nm[[dose]]),]
      } else {
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
            c.exp<-data.control[data.control[[id]] %in%  df[[id]],]
            c.un<-c.exp[c.exp[[unit]]==unit.cat
                        & !is.na(c.exp[[unit]]),]
            df.same<-rbind(df.un,c.un)
            df.all<-rbind(df.com,c.exp)
          }

        }
        df.nm<-switch(control.opt,same=df.same,all=df.all)
      }
      else {
        if (missing(data.control)) {
          df.nm<-df
        } else {
          if (missing(end)==FALSE && missing(end.cat)==FALSE) {
            c.end<-data.control[data.control[[end]]==end.cat &
                                  !is.na(data.control[[end]]) &
                                  !is.na(data.control[[dose]]),]
            c.exp<-c.end[c.end[[id]] %in%  df[[id]],]
          } else {
            c.exp<-data.control[data.control[[id]] %in%  df[[id]],]
          }
          df.nm<-rbind(df,c.exp)
        }
      }
      df.nm<-df.nm[!is.na(df.nm[[dose]]) & !is.na(df.nm[[response]]),]
      if (nrow(df.nm)==0) next
      if (missing(vars)==TRUE) {
        num.dose<-as.numeric(as.character(df.nm[[dose]]))
        logdose<-log10(num.dose)
        cons<-NULL
        df.nm$cons<-switch(type, dose=num.dose,log=logdose)
        dose.lab<-switch(type, dose="dose", log="log(dose)")
        resp<-sym(response)
        if (nrow(df.nm)==0) next
        dose.plot[[nano.cat]]<-ggplot(df.nm,aes(x=cons,y=!!resp)) +
          geom_point(shape=16,size=3)+ylab(paste(response)) +
          xlab(paste(dose.lab)) +
          scale_shape_manual(values = c(1,2))+ggtitle(nano.cat) + theme_bw()
        dplot<-marrangeGrob(dose.plot, nrow=nrow, ncol=ncol)
      }
      else if (missing(vars)==FALSE) {
        if (length(vars)==1) {
          subdata<-split(df.nm, lapply(df.nm[vars], addNA),drop=T)
        } else {
          subdata<-split(df.nm,lapply(df.nm[,vars], addNA),drop=T)
        }

        sub.plot<-vector(mode = "list")
        subs<-subdata
        for (j in 1:length(subs)) {
          com.data<-subs[[j]]
          num.dose<-as.numeric(as.character(com.data[[dose]]))
          logdose<-log10(num.dose)
          cons<-NULL
          com.data$cons<-switch(type, dose=num.dose,log=logdose)
          dose.lab<-switch(type, dose="dose", log="log(dose)")
          resp<-sym(response)
          if (nrow(com.data)==0) next
          sub.plot[[j]]<-ggplot(com.data, aes(x=cons,y=!!resp)) +
            geom_point(shape=16,size=3)+ylab(paste(response)) +
            xlab(paste(dose.lab)) +
            scale_shape_manual(values = c(1,2))+ggtitle(paste(nano.cat,".",names(subs[j]))) + theme_bw()
        }
        dose.plot[[nano.cat]]<-marrangeGrob(sub.plot, nrow=nrow, ncol=ncol)
      }
    }
    if (missing(vars)==TRUE) {
      return(dplot)
    } else {
      return(dose.plot)
    }
  }
}

