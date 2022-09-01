#' Calculate the test statistics and p-values for testing monotonic trend
#'
#' This function provides the value of the test statistics, p-values, plot of
#' the null distribution and the observed test statistics, and adjusted p-values
#' according to the method chosen.
#' @usage Isotest(data.nm, data.control, id, nano, response, dose, end, end.cat,
#'   unit, unit.cat, stat=c("E2","Williams","Marcus","M","ModM"), niter,
#'   method=p.adjust.methods, control.opt=c("same","all"), set.seed, vars)
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
#' @param stat Test statistics ("\code{E2}" for the global likelihood test,
#'   "\code{Williams}" for Williams test, "\code{Marcus}" for Marcus test,
#'   "\code{M}" for M test or "\code{ModM}" for modified M test)
#' @param niter Number of permutations
#' @param method Method used to adjust for the multiplicity (see \code{p.adjust.methods})
#' @param control.opt Option for the control doses if \code{unit} and
#'   \code{unit.cat} are specified. If only control doses with the same unit of
#'   measurement as the non-control ones are included, then specify
#'   "\code{same}" in the \code{control.opt}. If all control doses with any
#'   units of measurement are included, then specify "\code{all}".
#' @param set.seed Specify seed
#' @param vars Variable(s) used to subset the data
#' @return This function calculates the value for the test statistics, the
#'   p-value (and its plot) based on the permutation, the adjusted p-value and
#'   the more likely direction of the monotonic trend.
#' @details
#' \itemize{
#' \item{This function performs trend testing for each nanomaterial in the
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
#' \item{Trend testing can also be performed for subsets of data in each
#' nanomaterial, by specifying the variables used to split the data in
#' \code{vars}}.
#' }
#' @references Lin D., Pramana, S., Verbeke, T., and Shkedy, Z. (2015). IsoGene:
#'   Order-Restricted Inference for Microarray Experiments. R package version
#'   1.0-24. \url{https://CRAN.R-project.org/package=IsoGene}
#'
#'   Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, L. (editors).
#'   (2012) Modeling Doseresponse Microarray Data in Early Drug Development
#'   Experiments Using R. Springer.
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
#' # Perform test for the monotonic trend in geninvitro, with DNA STRAND BREAKS
#' # as the endpoint, concentrations measured in "ug/cm2" and control doses
#' # measured in any units of measurement:
#' #
#' \donttest{Isotest(data.nm=invitrodata, data.control=controldata, id="experimentID",
#'         nano="name", response="value", dose="concentration", end="endpoint",
#'         end.cat="DNA STRAND BREAKS", unit="concentration_unit", unit.cat="ug/cm2",
#'         stat="E2", niter=1000, method="BH", control.opt="all", set.seed=1234)}
#'
#' #' # Example 2:
#' # Split geninvitro according to the endpoint and unit of the concentration
#' # and then perform test for the monotonic trend for each subset of data:
#' #
#' \donttest{Isotest(data.nm=invitrodata, data.control=controldata,
#'         id="experimentID", nano="name", dose="concentration",
#'         response="value", stat="E2", niter=1000, method="BH", set.seed=1234,
#'         vars=c("endpoint","concentration_unit"))}

#' @importFrom dplyr filter
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @export
Isotest<-function(data.nm, data.control, id, nano, response, dose, end, end.cat,
                      unit, unit.cat, stat=c("E2","Williams","Marcus","M","ModM"),
                      niter, method=p.adjust.methods, control.opt=c("same","all"),
                      set.seed, vars){
  stat<-match.arg(stat)
  control.opt<-match.arg(control.opt)
  pval<-vector(mode="list")
  test<-vector(mode="list")
  mod<-vector(mode="list")
  snames<-vector(mode="list")
  for (nano.cat in unique(data.nm[[nano]])) {
    set.seed(set.seed)
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
      }
      else
      {
        if (missing(end)==FALSE && missing(end.cat)==FALSE) {
          c.end<-data.control[data.control[[end]]==end.cat &
                                !is.na(data.control[[end]]) &
                                !is.na(data.control[[dose]]),]
          c.exp<-c.end[c.end[[id]] %in%  df[[id]],]
          c.un<-c.exp[c.exp[[unit]]==unit.cat
                      & !is.na(c.exp[[unit]]),]
          df.same<-rbind(df.un,c.un)
          df.all<-rbind(df.com,c.exp)
        }
        else
        {
          c.exp<-data.control[data.control[[id]] %in%  df[[id]] &
                                !is.na(data.control[[dose]]),]
          #### above
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
      }
      else
      {
        if (missing(end)==FALSE && missing(end.cat)==FALSE) {
          c.end<-data.control[data.control[[end]]==end.cat &
                                !is.na(data.control[[end]]) &
                                !is.na(data.control[[dose]]),]
          c.exp<-c.end[c.end[[id]] %in%  df[[id]],]
        }
        else
        {
          c.exp<-data.control[data.control[[id]] %in%  df[[id]] &
                                !is.na(data.control[[dose]]),]
          ##### above
        }
        df.nm<-rbind(df,c.exp)
      }
    }

    df.nm<-df.nm[!is.na(df.nm[[dose]]) & !is.na(df.nm[[response]]),]
    if (nrow(df.nm)==0) next

    if (missing(vars)==TRUE) {
      tab1 <- table(droplevels(as.factor(df.nm[[dose]])))
      if (length(tab1)==1) next
      tab2 <- as.data.frame(tab1)
      if (!(0 %in% as.numeric(as.character(tab2$Var1)))) next
      com.data<-df.nm %>% group_by(df.nm[[dose]]) %>% filter(n() > 1)
      if (length(unique(com.data[[response]]))==1) next
      if (nrow(com.data)==0) next
      test[[nano.cat]]<-IsoTest.nm(data.nm=com.data,
                                   dose=dose,
                                   response=response,
                                   stat=stat)
      dev.new()
      pval[[nano.cat]]<-IsoPval.nm(data.nm=com.data,
                                   niter=niter,stat=stat,
                                   dose=dose,response=response,
                                   nano.cat=nano.cat)
    }
    else if (missing(vars)==FALSE) {
      df.nm<-df.nm[!is.na(df.nm[[dose]]) & !is.na(df.nm[[response]]),]
      if (length(vars)==1) {
        subdata<-split(df.nm,lapply(df.nm[vars], addNA),drop=T)
      } else {
        subdata<-split(df.nm,lapply(df.nm[,vars], addNA),drop=T)
      }
      test.sub<-vector(mode="list")
      pval.sub<-vector(mode="list")
      subnames<-vector(mode="list")
      for (j in 1:length(subdata)) {
        tab1<-table(droplevels(as.factor(subdata[[j]][[dose]])))
        if (length(tab1)==1) next
        tab2<-as.data.frame(tab1)
        if (!(0 %in% as.numeric(as.character(tab2$Var1)))) next
        com.data<-subdata[[j]] %>% group_by(subdata[[j]][[dose]]) %>% filter(n() > 1)
        if (length(unique(com.data[[response]]))==1) next
        if (nrow(com.data)==0) next
        test.sub[[j]]<-IsoTest.nm(data.nm=com.data,
                                  dose=dose,
                                  response=response,
                                  stat=stat)
        dev.new()
        pval.sub[[j]]<-IsoPval.nm(data.nm=com.data,
                                  niter=niter,stat=stat,
                                  dose=dose,response=response,
                                  nano.cat=paste(nano.cat,".",names(subdata[j])))
        subnames[[j]]<-paste(nano.cat,";",names(subdata[j]))
      }
      pval[[nano.cat]]<-pval.sub
      test[[nano.cat]]<-test.sub
      snames[[nano.cat]]<-subnames
    }

  }


  if (missing(vars)==TRUE) {
    pval.up<-lapply(pval,function(x) x[,1])
    pval.down<-lapply(pval,function(x) x[,2])
    t.stat.up<-lapply(test,function(x) x[,1])
    t.stat.down<-lapply(test,function(x) x[,2])
    tst.up<-lapply(t.stat.up,round,3)
    tst.dn<-lapply(t.stat.down,round,3)
    drc<-lapply(test,function(x) x[,3])
    if (length(pval) > 1) {
      adj.p.up <- p.adjust(pval.up, method, n = length(pval))
      adj.p.down <- p.adjust(pval.down, method, n = length(pval))
      adj.u <- lapply(adj.p.up, round, 3)
      adj.dn <- lapply(adj.p.down, round, 3)
      all <- cbind(tst.up, tst.dn, pval.up, pval.down, adj.u,
                   adj.dn, drc)
      colnames(all) <- c(paste("tstat.up:", stat),
                         paste("tstat.dn:", stat),
                         "pvalue.up", "pvalue.dn", paste("adj.p.up:",method),
                         paste("adj.p.dn:", method), "direction")
    }
    else {
      pvalue.up <- lapply(pval.up, round, 3)
      pvalue.dn <- lapply(pval.down, round, 3)
      all <- cbind(tst.up, tst.dn, pvalue.up, pvalue.dn, drc)
      colnames(all) <- c(paste("tstat.up:", stat),
                         paste("tstat.dn:", stat),
                         "pvalue.up", "pvalue.dn", "direction")
    }
    return(all)
  }
  else if (missing(vars)==FALSE) {
    pval.u<-lapply(unlist(pval, recursive = FALSE), `[`, 1)
    pval.d<-lapply(unlist(pval, recursive = FALSE), `[`, 2)
    t.stat.u<-lapply(unlist(test, recursive = FALSE), `[`, 1)
    t.stat.d<-lapply(unlist(test, recursive = FALSE), `[`, 2)
    pval.up<-as.vector(unlist(pval.u))
    pval.down<-as.vector(unlist(pval.d))
    t.stat.up<-as.vector(unlist(t.stat.u))
    t.stat.down<-as.vector(unlist(t.stat.d))
    tst.up<-lapply(t.stat.up,round,3)
    tst.dn<-lapply(t.stat.down,round,3)
    dr<-lapply(unlist(test, recursive = FALSE), `[`, 3)
    drc<-as.vector(unlist(dr))
    sub.names<-(unlist(lapply(unlist(snames, recursive = FALSE), `[`, 1)))
    sub.names2<-data.frame(sub.names)
    rownames(sub.names2)=NULL
    if (length(pval.up)>1) {
      adj.p.up<-p.adjust(pval.up, method, n = length(pval.up))
      adj.p.down<-p.adjust(pval.down, method, n = length(pval.down))
      adj.u<-lapply(adj.p.up,round,3)
      adj.dn<-lapply(adj.p.down,round,3)
      allvar=cbind(tst.up,tst.dn,pval.up,pval.down,adj.u,adj.dn,drc)
      rownames(allvar)<-NULL
      all<-cbind(sub.names2,allvar)
      colnames(all)<-c("name",paste("tstat.up:",stat),
                       paste("tstat.dn:",stat),"pvalue.up",
                       "pvalue.dn",paste("adj.p.up:",method),
                       paste("adj.p.dn:",method),"direction")
    } else {
      pvalue.up<-lapply(pval.up,round,3)
      pvalue.dn<-lapply(pval.down,round,3)
      allvar<-cbind(tst.up,tst.dn,pvalue.up,pvalue.dn,drc)
      rownames(allvar)<-NULL
      all<-cbind(sub.names2,allvar)
      colnames(all)<-c("name",paste("tstat.up:",stat),
                       paste("tstat.dn:",stat),"pvalue.up",
                       "pvalue.dn","direction")
    }
    return(all)
  }
}
