#' Create plots of the dose and the response differentiated by specified
#' variables
#'
#' This function generates scatter plots of the dose and the response for every
#' unique value of a certain variable, with the colour of the data points
#' differentiated according to the value of another variable.
#'
#' @usage nmplot.ncat(data.nm, data.control, id, nano, response, dose, end,
#'   end.cat, unit, unit.cat,  cat, x.cat, type=c("dose","log"),
#'   control.opt=c("same","all"), nrow=1, ncol=1)
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
#' @param cat Plot is generated for every unique value of \code{cat}
#' @param x.cat Variable used to differentiate the colour of the data points in
#'   the plot(s)
#' @param type Type of the dose to be plotted ("\code{dose}" for dose and
#'   "\code{log}" for log(dose))
#' @param control.opt Option for the control doses. If only control doses with
#'   the same unit of measurement as the non-control ones are included, then
#'   specify "\code{same}" in the \code{control.opt}. If all control doses with
#'   any units of measurement are included, then specify "\code{all}".
#' @param nrow Number of rows in the plotting space (default is 1)
#' @param ncol Number of columns in the plotting space (default is 1)
#' @return This function produces dose-response plots for every unique value of
#'   a certain variable, with different colours of data points based on the
#'   value of another variable
#' @details
#' \itemize{
#' \item{This function generates plots for each nanomaterial in the dataset. The
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
#' controldata<-SubsetData(data=geninvitro, x="name", x.cat=c("control", "Control",
#'              "medium", "medium + BSA", "untreated"))
#'
#' # Exclude controls (which are named differently) from geninvitro dataset:
#' invitrodata<-SubsetData(data=geninvitro, x="name", x.cat=c("control", "Control",
#'              "medium", "medium + BSA", "untreated"), include=FALSE)
#'
#' # Generate dose-response plots for geninvitro, with DNA STRAND BREAKS as
#' # the endpoint and concentrations measured in "ug/cm2" (plot is generated for
#' # each study provider, with different colours represent different
#' # experiments):
#' #
#' \donttest{nmplot.ncat(data.nm=invitrodata, data.control=controldata,
#'             id="experimentID", nano="name", response="value",
#'             dose="concentration", end="endpoint", end.cat="DNA STRAND BREAKS",
#'             unit="concentration_unit", unit.cat="ug/cm2", cat="studyprovider",
#'             x.cat="experimentID", type="dose", control.opt="same", nrow=1,
#'             ncol=1)}
#'
#' @import ggplot2
#' @importFrom gridExtra marrangeGrob
#' @export
nmplot.ncat<-function(data.nm, data.control, id, nano, response, dose, end,
                          end.cat, unit, unit.cat, cat, x.cat,
                          type=c("dose","log"), control.opt=c("same","all"),
                          nrow=1, ncol=1) {
  if (missing(nano)==TRUE) {stop("argument \"nano\" is missing")
  } else if (missing(response)==TRUE) {stop("argument \"response\" is missing")
  } else if (missing(dose)==TRUE) {stop("argument \"dose\" is missing")
  } else if (missing(end)==TRUE) {stop("argument \"end\" is missing")
  } else if (missing(unit)==TRUE) {stop("argument \"unit\" is missing")
  } else if (missing(unit.cat)==TRUE) {stop("argument \"unit.cat\" is missing")
  } else if (missing(cat)==TRUE) {stop("argument \"cat\" is missing")
  } else
  {
    type<-match.arg(type)
    control.opt<-match.arg(control.opt)
    dose.plot<-vector(mode="list")
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
      df.nm<-df.nm[!is.na(df.nm[[dose]]) &
                     !is.na(df.nm[[response]]),]
      num.dose<-as.numeric(as.character(df.nm[[dose]]))
      logdose<-log10(num.dose)
      cons<-NULL
      df.nm$cons<-switch(type,dose=num.dose,log=logdose)
      dose.lab<-switch(type,dose="dose",log="log(dose)")
      resp<-sym(response)
      cats<-sym(x.cat)
      if (nrow(df.nm)==0) next
      plot.cat<-vector(mode="list")
      tab<-vector(mode = "list")
      for (u.cat in unique(df.nm[[cat]])) {
        if (is.na(u.cat)==T) {
          df.cat<-df.nm[is.na(df.nm[[cat]])==TRUE,]
        } else {
          df.cat<-df.nm[df.nm[[cat]]==u.cat & is.na(df.nm[[cat]])==FALSE,]
        }
        tab[[u.cat]]<-unique(df.cat[[x.cat]])
        if (length(tab[[u.cat]])==1 & any(is.na(tab[[u.cat]]))) {
          plot.cat[[u.cat]]<-ggplot(df.cat,
                                    aes(x=cons,y=!!resp,colour=!!cats)) +
            geom_point(shape=16,size=3) + ylab(paste(response)) +
            xlab(paste(dose.lab,"(",unit.cat,")")) +
            scale_shape_manual(values = c(1,2)) + labs(colour = paste(x.cat),
                                                       caption = paste(x.cat,": NA"))+
            scale_colour_brewer(palette = "Set1",
                                aesthetics="colour",
                                na.value="grey45")+ ggtitle(paste(nano.cat,u.cat,sep="\n")) + theme_bw()
        }else{
          plot.cat[[u.cat]]<-ggplot(df.cat,
                                    aes(x=cons,y=!!resp,colour=!!cats)) +
            geom_point(shape=16,size=3)+ylab(paste(response)) +
            xlab(paste(dose.lab,"(",unit.cat,")")) +
            scale_shape_manual(values = c(1,2))+
            labs(colour = paste(x.cat))+ggtitle(paste(nano.cat,u.cat,sep="\n")) + theme_bw()
        }
      }
      plots<-marrangeGrob(plot.cat, nrow=nrow, ncol=ncol)
      dose.plot[[nano.cat]]<-plots
    }
    return(dose.plot)
  }
}
