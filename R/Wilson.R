

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Calculate Qhats (mean of response for each imputed dataset)
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of binary response variable
#'
#' @return Qhats: vector of mean response for each imputed dataset
#' @export
#' @importFrom dplyr all_of
#' @importFrom dplyr select
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Qhats(imp, "hyp")
#'
Qhats <- function(mids_obj,response) {
  if ("constant" %in% mids_obj$loggedEvents$meth) {
    stop(paste("MICE unable to impute",response," due to constant observed response values."))
  }

  qhats = lapply(1:mids_obj$m,
                 function(i) mice::complete(mids_obj,i) %>%
                   dplyr::select(all_of(response)) %>%
                   colMeans()) %>% unlist()

  return(qhats)
}

#' Calculate Qbar (average response over MICE datasets)
#'
#' @param response string name of binary response variable
#' @param mids_obj mids object created by mice package
#'
#' @return Qbar: the average response over MICEd datasets.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Qbar(imp, "hyp")
#'
Qbar <- function(mids_obj, response) {
  qbar = Qhats(mids_obj, response) %>% mean()

  return(qbar)
}

#' Calculate Uhats (variance for each imputed dataset)
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of binary response variable
#'
#' @return Uhats: vector of response variances for each imputed dataset
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Uhats(imp,"hyp")
#'
Uhats <- function(mids_obj, response) {
  qhats = Qhats(mids_obj, response)
  return(qhats*(1-qhats)/(mids_obj$data %>% nrow()))
}

#' Calculate Ubar (average response variance over MICE datasets)
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#'
#' @return Ubar: average response variance over MICE datasets
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Ubar(imp, "hyp")
#'
Ubar <- function(mids_obj, response) {
  uhats = Uhats(mids_obj, response)

  return(sum(uhats)/mids_obj$m)
}


#' Calculate between-imputation variance of the response mean
#' \deqn{\frac{\sum (\hat{Q}_l-\bar{Q})}{m-1}}
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#'
#' @return Bm: the between-dataset variance of the response mean
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Bm(imp, "hyp")
#'
Bm <- function(mids_obj, response) {
  qbar = Qbar(mids_obj, response)
  qhats = Qhats(mids_obj, response)

  return(sum((qhats-qbar)^2)/(mids_obj$m-1))
}

#' Estimate variance of proportion point estimate \eqn{\bar{Q}_m}
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#'
#' @return variance of proportion point estimate
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Tm(imp, "hyp")
#'
Tm <- function(mids_obj, response) {
  ubar = Ubar(mids_obj, response)
  bm = Bm(mids_obj, response)

  return((1 + 1/mids_obj$m)*bm + ubar)
}

#' Helper function for getting rm, a key component for
#' calculating degrees of freedom and the wilson CI directly
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#'
#' @return rm
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Rm(imp, "hyp")
#'
Rm <- function(mids_obj, response) {
  ubar = Ubar(mids_obj, response)
  bm = Bm(mids_obj, response)
  m = mids_obj$m

  if(bm==0) {
    rm = 0
    message("Bm is 0, imputed and observed values of bin. variable are the same.")
  }
  else {
    rm = (1 + 1/m)*bm/ubar
  }

  return(rm)
}

#' Calculate degrees of freedom used in calculating
#' confidence intervals of t-distributed proportion point
#' estimate \eqn{\bar{Q}_m}
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#'
#' @return degrees of freedom
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' dof(imp, "hyp")
#'
dof <- function(mids_obj, response) {
  m = mids_obj$m
  rm = Rm(mids_obj, response)

  return((m - 1)*(1 + 1/rm)^2)
}

#' Calculates the specified Wilson CI of a binomial proportion
#' variable, given imputed data sets.
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable (must be 0-1 valued)
#' @param ci_level desired confidence interval level (defaults to 95%)
#'
#' @return two-length vector of Wilson lower CI and upper CI
#' @export
#' @importFrom dplyr mutate
#'
#' @examples
#' imp = mice::mice(mice::nhanes %>% dplyr::mutate(hyp = hyp-1))
#' mi_wilson(imp, "hyp", 0.95)
#'
mi_wilson <- function(mids_obj, response, ci_level=0.95) {

  #if confidence interval is invalid
  if(ci_level<=0 | ci_level>= 1) {
    stop("CI level must be between 0 and 1.")
  }

  #if response is not 0-1 valued
  if(all(lapply(mids_obj$data %>% select(all_of(response)),
                function(x) x %in% c(NA,0,1)) %>% unlist())==FALSE) {
    stop(paste(response,"must be 0-1 binary encoded."))
  }

  qbar = Qbar(mids_obj, response)
  print(paste("Qbar: ",qbar))
  rm = Rm(mids_obj, response)
  print(paste("RM: ",rm))

  #if response variable has one value only
  if(rm==0) {
    df = Inf
  }
  else {
    df = dof(mids_obj, response)
  }

  t_score = stats::qt(ci_level, df)
  n = mids_obj$data %>% nrow()

  tquad = t_score^2 / n + t_score^2 * rm /n
  center = (2 * qbar + tquad)/(2*(1 + tquad))
  half_width = sqrt( (2*qbar + tquad)^2 / (4*(1+tquad)^2) -
                       qbar^2 / (1+tquad) )

  return(c(center - half_width, center + half_width))
}

#' Calculates the specified Wald CI of a binomial proportion
#' variable, given imputed data sets. TODO: create edited nhanes dataset
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable (must be 0-1 valued)
#' @param ci_level desired confidence interval level (defaults to 95%)
#'
#' @return two-length vector of Wald lower CI and upper CI
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' imp = mice::mice(mice::nhanes %>% dplyr::mutate(hyp = hyp-1))
#' mi_wald(imp, "hyp", 0.95)
#'
mi_wald <- function(mids_obj, response, ci_level=0.95) {

  if(ci_level<=0 | ci_level>= 1) {
    stop("CI level must be between 0 and 1.")
  }

  #if response is not 0-1 valued
  if(all(lapply(mids_obj$data %>% select(all_of(response)),
                function(x) x %in% c(NA,0,1)) %>% unlist())==FALSE) {
    stop(paste(response,"must be 0-1 binary encoded."))
  }

  qbar = Qbar(mids_obj, response)
  tm = Tm(mids_obj, response)
  df = dof(mids_obj, response)
  t_score = stats::qt(ci_level, df)

  return(c(qbar - sqrt(tm)*t_score, qbar + sqrt(tm)*t_score))
}
