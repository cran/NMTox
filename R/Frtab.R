#' Generate the number of unique values, the list of unique values, or the number
#' of observations of specified variable(s)
#'
#' This function produces the number of unique values, the list of unique
#' values, or the number of observations of certain variable(s) in the dataset
#'
#' @usage Frtab(data, x, cat, val, opt=c("un", "list", "obs", "ls.obs"),
#'   na.rm=FALSE)
#' @param data Dataset
#' @param x Variable(s) to be explored
#' @param cat Variable used to group or subset the data
#' @param val Specific value of variable \code{cat}. If \code{val} is specified,
#'   then only result for this value will be shown.
#' @param opt Options for the result: If \code{opt} is \code{un}, then the
#'   number of unique values of \code{x} will be counted. If \code{opt} is
#'   \code{list}, then the list of unique values of \code{x} will be
#'   generated. If \code{opt} is \code{obs}, then the number of observations of
#'   \code{x} will be counted. If \code{opt} is \code{ls.obs}, then the number
#'   of observations for each unique value will be generated.
#' @param na.rm If \code{na.rm} is \code{FALSE} (the default), then missing
#'   values are not removed from the dataset. If \code{na.rm}=\code{TRUE}, then
#'   missing values are removed from the dataset
#' @return This function generates the number of unique values, the list of
#'   unique values, or the number of observations of specified variable(s).
#' @examples
#'
#' # List of nanomaterial in geninvitro dataset:
#' Frtab(data=geninvitro, x="name", opt="list")
#'
#' # How many types of endpoint are measured in geninvitro data?
#' Frtab(data=geninvitro, x="endpoint", opt="un")
#'
#' # How many observations are available for variable "unit"?
#' Frtab(data=geninvitro, x="unit", opt="obs")
#'
#' # How many types of endpoint and nanomaterial are available in geninvitro data?
#' Frtab(data=geninvitro, x=c("name","endpoint"), opt="un")
#'
#' # How many types of endpoint are available for each nanomaterial
#' # in geninvitro data?
#' Frtab(data=geninvitro, x="endpoint", cat="name", opt="un")
#'
#' # How many observations are available for each endpoint in each nanomaterial?
#' Frtab(data=geninvitro, x=c("name","endpoint"), opt="ls.obs")
#'
#' # How many observations with "DNA STRAND BREAKS" as the endpoint are available
#' # for each nanomaterial in geninvitro data?
#' Frtab(data=geninvitro, x="name", cat="endpoint", val="DNA STRAND BREAKS",
#' opt="ls.obs")
#'
#'
#' @import dplyr
#' @importFrom tidyr drop_na
#' @importFrom stats na.omit
#' @export
Frtab<-function(data, x, cat, val, opt=c("un", "list", "obs", "ls.obs"), na.rm=FALSE) {
  opt<-match.arg(opt)
  if(missing(cat)==F && missing(val)==F) {
    data2<-data[data[[cat]]==val & !is.na(data[[cat]]),]
    un<-data2 %>% summarize_at(x, n_distinct, na.rm=na.rm)
    un.df<-as.data.frame(un)
    obs<-data2 %>% summarise_at(x, ~sum(!is.na(.)))
    obs.df<-as.data.frame(obs)

    if (na.rm==T){
      lis<-unique(na.omit(data2[,x]))
      lis.ob<- data2 %>% drop_na(all_of(x)) %>% count(across(all_of(x)))
    } else {
      lis<-unique(data2[,x])
      lis.ob<- data2 %>% count(across(all_of(x)))
    }

  }
  else if (missing(cat)==T && missing(val)==T) {
    un<-data %>% summarize_at(x, n_distinct, na.rm=na.rm)
    un.df<-as.data.frame(un)

    obs<-data %>% summarise_at(x, ~sum(!is.na(.)))
    obs.df<-as.data.frame(obs)

    if(na.rm==T){
      lis<-unique(na.omit(data[,x]))
      lis.ob<- data %>% drop_na(all_of(x)) %>% count(across(all_of(x)))
    } else {
      lis<-unique(data[,x])
      lis.ob<- data %>% count(across(all_of(x)))
    }
  }
  else if (missing(cat)==F && missing(val)==T) {
    cats<-sym(cat)
    if (na.rm==T) {
      un<-data %>%
        filter(!is.na(!!cats)) %>%
        group_by(!!cats) %>%
        summarize_at(x, n_distinct, na.rm=na.rm)

      obs<-data %>%
        filter(!is.na(!!cats)) %>%
        group_by(!!cats) %>%
        summarise_at(x, ~sum(!is.na(.)))

      lis<-unique(na.omit(data[,c(cat,x)]))

      lis.ob<- data %>% filter_at(c(cat,x),~!is.na(.)) %>%
        count(across(all_of(c(cat,x))))

    } else {
      un<-data %>%
        group_by(!!cats) %>%
        summarize_at(x, n_distinct, na.rm=na.rm)

      obs<-data %>%
        group_by(!!cats) %>%
        summarise_at(x, ~sum(!is.na(.)))

      lis<-unique(data[,c(cat,x)])
      lis.ob<- data %>% count(across(all_of(c(cat,x))))
    }
    un.df<-as.data.frame(un)
    obs.df<-as.data.frame(obs)
  }

  out<-switch (opt, "un" = un.df, "obs" = obs.df, "list" = lis, "ls.obs" =lis.ob)
  return(out)
}
