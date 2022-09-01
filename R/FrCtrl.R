#' Create a frequency table of the control values
#'
#' This function provides a table containing the frequencies of the control
#' values and the number of available data in the specified dataset. Since it is
#' possible to have a dose with different units of measurement in one dataset,
#' this function helps in showing how many control values are available for each
#' nanomaterial for a specific unit of measurement of the dose.
#' @usage FrCtrl(data.nm, data.control, id, nano, dose, end, end.cat, unit,
#'   unit.cat, vars)
#' @param data.nm Data containing the result of toxicity study
#' @param data.control Data of control values
#' @param id Identifier of the experiment
#' @param nano Name of the nanomaterial
#' @param dose Dose or concentration
#' @param end Toxicity endpoint
#' @param end.cat Specific toxicity endpoint of interest
#' @param unit Unit of measurement of the dose
#' @param unit.cat Specific unit of measurement of the dose
#' @param vars Variables used to subset the data
#' @return The value returned from \code{FrCtrl} is a table containing:
#' \itemize{
#'   \item{\code{Freq(Dose=0.0).same}: frequency of control values with dose
#'   measured in \code{unit.cat}}
#'   \item{\code{Freq(obs).same}: number of observations with dose measured in
#'   \code{unit.cat}}
#'   \item{\code{Freq(Dose=0.0).all}: frequency of control values with dose
#'   measured \code{unit.cat} and in other units of measurement}
#'   \item{\code{Freq(obs).all}: number of observations with dose measured in
#'   \code{unit.cat} and control dose measured in any units},
#'   }
#'   for each nanomaterial or subset of data
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
#' \item{The exploration can also be performed on the subsets of data by
#' splitting the data of each nanomaterial according to the variable(s)
#' specified in \code{vars}}
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
#' # Frequency of controls for each nanomaterial in geninvitro dataset
#' # with DNA STRAND BREAKS as the toxicity endpoint and concentration
#' # measured in "ug/cm2":
#' #
#' FrCtrl(data.nm=invitrodata, data.control=controldata, nano="name",
#'        end="endpoint", end.cat="DNA STRAND BREAKS", id="experimentID",
#'        dose="concentration", unit="concentration_unit", unit.cat="ug/cm2")
#'
#' # Frequency of controls for each cell type in each nanomaterial
#' # (in geninvitro dataset) with DNA STRAND BREAKS as the toxicity
#' # endpoint and concentration measured in "ug/cm2":
#' #
#' FrCtrl(data.nm=invitrodata, data.control=controldata, nano="name",
#'        end="endpoint", end.cat="DNA STRAND BREAKS", id="experimentID",
#'        dose="concentration", unit="concentration_unit", unit.cat="ug/cm2",
#'        vars="celltype")
#' @import dplyr
#' @export
FrCtrl<-function(data.nm, data.control, id, nano, dose, end, end.cat,
                  unit, unit.cat, vars){
  if (missing(nano)==TRUE) {stop("argument \"nano\" is missing")
  } else if (missing(dose)==TRUE) {stop("argument \"dose\" is missing")
  } else if (missing(unit)==TRUE) {stop("argument \"unit\" is missing")
  } else
  {
    dose.fr<-vector(mode = "list")
    dose.frc<-vector(mode = "list")
    dose.fr2<-vector(mode = "list")
    dose.frc2<-vector(mode = "list")
    tabs<-vector(mode = "list")
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
      if (missing(vars)==TRUE) {
        if (nrow(df.same)==0 || nrow(df.all)==0) next
        sub.same<-df.same[as.numeric(as.character(df.same[[dose]]))==0.0
                          & !is.na(df.same[[dose]]),]
        dose.fr[[nano.cat]]<-table(droplevels(as.factor(as.character(sub.same[[dose]]))))
        dose.frc[[nano.cat]]<-nrow(df.same)
        sub.all<-df.all[as.numeric(as.character(df.all[[dose]]))==0.0
                        & !is.na(df.all[[dose]]),]
        dose.fr2[[nano.cat]]<-table(droplevels(as.factor(as.character(sub.all[[dose]]))))
        dose.frc2[[nano.cat]]<-nrow(df.all)
      } else {
        if (length(vars)==1) {
          subdata.same<-split(df.same, lapply(df.same[vars], addNA),drop=T)
          subdata.all<-split(df.all, lapply(df.all[vars], addNA),drop=T)
        } else {
          subdata.same<-split(df.same,lapply(df.same[,vars], addNA),drop=T)
          subdata.all<-split(df.all,lapply(df.all[,vars], addNA),drop=T)
        }
        subdose.fr<-vector(mode = "list")
        subdose.frc<-vector(mode = "list")
        for (j in 1:length(subdata.same)) {
          subs.same<-subdata.same[[j]]
          sub.same<-subs.same[as.numeric(as.character(subs.same[[dose]]))==0.0
                              & !is.na(subs.same[[dose]]),]
          subdose.fr[[j]]<-table(droplevels(as.factor(as.character(sub.same[[dose]]))))
          subdose.frc[[j]]<-nrow(subs.same)
        }
        subdose.fr2<-vector(mode = "list")
        subdose.frc2<-vector(mode = "list")
        for (j in 1:length(subdata.all)) {
          subs.all<-subdata.all[[j]]
          sub.all<-subs.all[as.numeric(as.character(subs.all[[dose]]))==0.0
                            & !is.na(subs.all[[dose]]),]
          subdose.fr2[[j]]<-table(droplevels(as.factor(as.character(sub.all[[dose]]))))
          subdose.frc2[[j]]<-nrow(subs.all)
        }
        sub<-data.frame()
        for(i in seq_along(subdose.fr)) {
          if(length(subdose.fr[[i]])==0){
            sub[i,1]<-nano.cat
            sub[i,2]<-paste(names(subdata.same)[i])
            sub[i,3]<-NA
          }else{
            sub[i,1]<-nano.cat
            sub[i,2]<-paste(names(subdata.same)[i])
            sub[i,3]<-subdose.fr[i][[1]]
          }
        }
        colnames(sub)<-c("Nanomaterial","Subset","Freq(Dose=0.0).same")
        sub2<-data.frame()
        for(i in seq_along(subdose.frc)) {
          sub2[i,1]<-nano.cat
          sub2[i,2]<-paste(names(subdata.same)[i])
          sub2[i,3]<-subdose.frc[i][[1]]
        }
        colnames(sub2)<-c("Nanomaterial","Subset","Freq(obs).same")
        sub3<-data.frame()
        for(i in seq_along(subdose.fr2)) {
          if(length(subdose.fr2[[i]])==0){
            sub3[i,1]<-nano.cat
            sub3[i,2]<-paste(names(subdata.all)[i])
            sub3[i,3]<-NA
          }else{
            sub3[i,1]<-nano.cat
            sub3[i,2]<-paste(names(subdata.all)[i])
            sub3[i,3]<-subdose.fr2[i][[1]]
          }
        }
        colnames(sub3)<-c("Nanomaterial","Subset","Freq(Dose=0.0).all")
        sub4<-data.frame()
        for(i in seq_along(subdose.frc2)) {
          sub4[i,1]<-nano.cat
          sub4[i,2]<-paste(names(subdata.all)[i])
          sub4[i,3]<-subdose.frc2[i][[1]]
        }
        colnames(sub4)<-c("Nanomaterial","Subset","Freq(obs).all")
        tab1<-merge(sub,sub2,by=c("Nanomaterial","Subset"))
        tab2<-merge(sub3,sub4,by=c("Nanomaterial","Subset"))
        tabs[[nano.cat]]<-merge(tab1,tab2,all=TRUE)
      }
    }
    if (missing(vars)==TRUE) {
      sub<-data.frame()
      for(i in seq_along(dose.fr)) {
        if(length(dose.fr[[i]])==0){
          sub[i,1]<-paste(names(dose.fr)[i])
          sub[i,2]<-NA
        }else{
          sub[i,1]<-paste(names(dose.fr)[i])
          sub[i,2]<-dose.fr[i][[1]]
        }
      }
      colnames(sub)<-c("Nanomaterial","Freq(Dose=0.0).same")
      sub2<-data.frame()
      for(i in seq_along(dose.frc)) {
        sub2[i,1]<-paste(names(dose.frc)[i])
        sub2[i,2]<-dose.frc[i][[1]]
      }
      colnames(sub2)<-c("Nanomaterial","Freq(obs).same")
      sub3<-data.frame()
      for(i in seq_along(dose.fr2)) {
        if(length(dose.fr2[[i]])==0){
          sub3[i,1]<-paste(names(dose.fr2)[i])
          sub3[i,2]<-NA
        }else{
          sub3[i,1]<-paste(names(dose.fr2)[i])
          sub3[i,2]<-dose.fr2[i][[1]]
        }
      }
      colnames(sub3)<-c("Nanomaterial","Freq(Dose=0.0).all")
      sub4<-data.frame()
      for(i in seq_along(dose.frc2)) {
        sub4[i,1]<-paste(names(dose.frc2)[i])
        sub4[i,2]<-dose.frc2[i][[1]]
      }
      colnames(sub4)<-c("Nanomaterial","Freq(obs).all")
      tab1<-merge(sub,sub2,by="Nanomaterial")
      tab2<-merge(sub3,sub4,by="Nanomaterial")
      tab<-merge(tab1,tab2,by="Nanomaterial")
      return(tab)
    } else {
      tabs2<- bind_rows(tabs)
      return(tabs2)
    }
  }
}
