#' Create plot(s) of the observations, the sample means and the fitted isotonic
#' regression curve
#'
#' This function creates plot(s) of the observations, the sample means and the
#' fitted isotonic regression curve for one or more nanomaterials simultaneously
#'
#' @usage Isoplot(data.nm, data.control, id, nano, response, dose, end, end.cat,
#'   unit, unit.cat, dose.type=c("dose","log"), type = c("continuous",
#'   "ordinal"), control.opt=c("same","all"), add.curve = TRUE, vars, nrow=1,
#'   ncol=1, xlabel="Dose", ylabel="Response")
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
#' @param dose.type Type of the dose to be plotted ("dose" for dose and "log"
#'   for log(dose))
#' @param type Type of the dose (continuous or ordinal)
#' @param control.opt Option for the control doses if \code{unit} and
#'   \code{unit.cat} are specified. If only control doses with the same unit of
#'   measurement as the non-control ones are included, then specify
#'   "\code{same}" in the \code{control.opt}. If all control doses with any
#'   units of measurement are included, then specify "\code{all}".
#' @param add.curve Adding curve to the plot
#' @param vars Variable(s) used to subset the data
#' @param nrow Number of row in the plotting space
#' @param ncol Number of column in the plotting space
#' @param xlabel Label for x-axis
#' @param ylabel Label for y-axis
#' @return This function produces plot(s) consisting of observation data points,
#'   sample mean for each dose and fitted isotonic regression curve
#' @details
#' \itemize{
#' \item{This function performs data exploration for each nanomaterial in the
#' dataset (or for each subset of data). The different types of nanomaterials
#' are identified by their names. Therefore, if some control values are named
#' differently (see: \code{\link{geninvitro}} dataset and the \code{Examples}),
#' a separate dataset containing only these values first needs to be created.
#' Controls in the new dataset can be linked to the non-control observations
#' belonging to the same experiment through the identifier of the experiment
#' (the linking is performed inside this function). In this situation, it is
#' necessary to have an indicator that can identify different experiments (such
#' as experiment ID).}
#' \item{If all controls in the dataset are named according to the related
#' nanomaterial names, \code{data.control} and \code{id} do not need to be
#' specified.}
#' \item{If doses used in the experiment are all measured in the same unit of
#' measurement, then specify "\code{same}" in \code{control.opt}}.
#' \item{Dose-response plot can also be generated for subsets of data in each
#' nanomaterial by specifying the variables used to split the data in
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
#' # Generate dose-response plots for geninvitro, with DNA STRAND BREAKS
#' # as the endpoint, concentrations measured in "ug/cm2"
#' # and control doses measured in any units of measurement:
#' #
#' \donttest{Isoplot(data.nm=invitrodata, data.control=controldata, id="experimentID",
#'         nano="name", response="value", dose="concentration", end="endpoint",
#'         end.cat="DNA STRAND BREAKS", unit="concentration_unit", unit.cat="ug/cm2",
#'         dose.type="dose", control.opt="all")}
#'
#' # Example 2:
#' # Split geninvitro data according to the cell type, method, study provider and
#' # unit of the concentration and generate dose-response plot for each subset
#' # of data with DNA STRAND BREAKS as the endpoint:
#' #
#' \donttest{Isoplot(data.nm=invitrodata, data.control=controldata, id="experimentID",
#'         nano="name", response="value", dose="concentration", end="endpoint",
#'         end.cat="DNA STRAND BREAKS",  dose.type="dose",
#'         vars=c("celltype","method","studyprovider","concentration_unit"),
#'         nrow=2, ncol=2)}
#'
#' @references Lin D., Pramana, S., Verbeke, T., and Shkedy, Z. (2015). IsoGene:
#'   Order-Restricted Inference for Microarray Experiments. R package version
#'   1.0-24. \url{https://CRAN.R-project.org/package=IsoGene}
#'
#'   Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, L. (editors).
#'   (2012) Modeling Doseresponse Microarray Data in Early Drug Development
#'   Experiments Using R. Springer.
#' @importFrom dplyr filter
#' @import ggplot2
#' @export
Isoplot<-function(data.nm, data.control, id, nano, response, dose,
                      end, end.cat, unit, unit.cat, dose.type=c("dose","log"),
                      type = c("continuous", "ordinal"),
                      control.opt=c("same","all"), add.curve = TRUE,
                      vars, nrow=1, ncol=1, xlabel="Dose", ylabel="Response"){
  dose.type<-match.arg(dose.type)
  type<-match.arg(type)
  control.opt<-match.arg(control.opt)
  df.nm<-vector(mode="list")
  dose.plot<-vector(mode="list")
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
          c.exp<-data.control[data.control[[id]] %in%  df[[id]],]
          c.un<-c.exp[c.exp[[unit]]==unit.cat
                      & !is.na(c.exp[[unit]]),]
          df.same<-rbind(df.un,c.un)
          df.all<-rbind(df.com,c.exp)
        }
      }
      df.nm<-switch(control.opt,same=df.same,all=df.all)
    }
    else
    {
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
      df.nm$cons<-switch(dose.type,dose=num.dose,log=logdose)
      dose.lab<-switch(dose.type,dose="dose",log="log(dose)")
      if (nrow(df.nm)==0) next
      tab1 <- table(droplevels(as.factor(df.nm[[dose]])))
      if (length(tab1) == 1) next
      tab2 <- as.data.frame(tab1)
      if (!(0 %in% as.numeric(as.character(tab2$Var1)))) next
      com.data <- df.nm %>% group_by(df.nm[[dose]]) %>% filter(n() > 1)
      if (nrow(com.data) == 0) next
      dose.plot[[nano.cat]] <- IsoPlot.nm(data.nm = com.data, dose = "cons",
                                          response = response, type = type,
                                          add.curve = add.curve, nano.cat = nano.cat,
                                          xlabel, ylabel)
      plots <- marrangeGrob(dose.plot, nrow = nrow, ncol = ncol)
    }
    else if
    (missing(vars)==FALSE) {
      if (length(vars)==1) {
        subdata<-split(df.nm,lapply(df.nm[vars], addNA),drop=T)
      } else {
        subdata<-split(df.nm,lapply(df.nm[,vars], addNA),drop=T)
      }
      subs<-subdata
      sub.plot<-vector(mode="list")
      for (j in 1:length(subs)) {
        subs2<-subs[[j]]
        subs2 <-subs2[!is.na(subs2[[dose]]) & !is.na(subs2[[response]]), ]
        num.dose<-as.numeric(as.character(subs2[[dose]]))
        logdose<-log10(num.dose)
        subs2$cons<-switch(dose.type,dose=num.dose,log=logdose)
        dose.lab<-switch(dose.type,dose="dose",log="log(dose)")
        if (nrow(subs2)==0) next
        tab1 <- table(droplevels(as.factor(subs2[[dose]])))
        if (length(tab1) == 1) next
        tab2 <- as.data.frame(tab1)
        if (!(0 %in% as.numeric(as.character(tab2$Var1)))) next
        com.data <- subs2 %>% group_by(subs2[[dose]]) %>% filter(n() > 1)
        if (nrow(com.data) == 0) next
        sub.plot[[j]]<-IsoPlot.nm(data.nm = com.data, dose = "cons",
                                  response = response, type = type,
                                  add.curve = add.curve,
                                  nano.cat = paste(nano.cat,".",names(subs[j])),
                                  xlabel, ylabel)
      }
      sub.plot2 <- Filter(Negate(is.null), sub.plot)
      dose.plot[[nano.cat]]<-marrangeGrob(grobs=sub.plot2, nrow=nrow, ncol=ncol)
    }
  }
  if (missing(vars)==TRUE) {
    return(plots)
  } else {
    return(dose.plot)}
}

