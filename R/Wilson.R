

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Calculate Qhats (means of response for each imputed dataset)
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of binary response variable
#'
#' @return Qhats: vector of response means for each imputed dataset
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
#' @param qhats vector of Qhats(response means for each imputed dataset)
#'
#' @return Qbar: the average response over MICEd datasets.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' Qbar(qhats)
#'
Qbar <- function(qhats) {
  return(qhats %>% mean())
}

#' Calculate Uhats (variance for each imputed dataset)
#'
#' @param qhats vector of Qhats(means of response for each imputed dataset)
#' @param nrow number of observations in the imputed dataset
#'
#' @return Uhats: vector of response variances for each imputed dataset
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' nrow = imp$data %>% nrow()
#' Uhats(qhats, nrow)
#'
Uhats <- function(qhats, nrow) {
  return(qhats*(1-qhats)/nrow)
}

#' Calculate Ubar (average response variance over MICE datasets)
#'
#' @param qhats vector of Qhats(means of response for each imputed dataset)
#' @param m number of imputed datasets
#' @param nrow number of observations in the imputed dataset
#'
#' @return Ubar: average response variance over MICE datasets
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' m = imp$m
#' nrow = imp$data %>% nrow()
#' Ubar(qhats, m, nrow)
#'
Ubar <- function(qhats, m, nrow) {
  uhats = Uhats(qhats, nrow)

  return(sum(uhats)/m)
}


#' Calculate between-imputation variance of the response mean
#' \deqn{\frac{\sum (\hat{Q}_l-\bar{Q})}{m-1}}
#'
#' @param qhats vector of Qhats(means of response for each imputed dataset)
#' @param m number of imputed datasets
#'
#' @return Bm: the between-dataset variance of the response mean
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' m = imp$m
#' Bm(qhats, m)
#'
Bm <- function(qhats, m) {
  qbar = Qbar(qhats)

  return(sum((qhats-qbar)^2)/(m-1))
}

#' Estimate variance of proportion point estimate \eqn{\bar{Q}_m}
#'
#' @param qhats vector of Qhats(means of response for each imputed dataset)
#' @param m number of imputed datasets
#' @param nrow number of observations in the imputed dataset
#'
#' @return variance of proportion point estimate
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' m = imp$m
#' nrow = imp$data %>% nrow()
#' Tm(qhats, m, nrow)
#'
Tm <- function(qhats, m, nrow) {
  ubar = Ubar(qhats, m, nrow)
  bm = Bm(qhats, m)

  return((1 + 1/m)*bm + ubar)
}

#' Helper function for getting rm, a key component for
#' calculating degrees of freedom and the wilson CI directly
#'
#' @param qhats vector of Qhats(means of response for each imputed dataset)
#' @param m number of imputed datasets
#' @param nrow number of observations in the imputed dataset
#'
#' @return rm
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' m = imp$m
#' nrow = imp$data %>% nrow()
#' Rm(qhats, m, nrow)
#'
Rm <- function(qhats, m, nrow) {
  ubar = Ubar(qhats, m, nrow)
  bm = Bm(qhats, m)

  if(bm==0) {
    rm = 0
    #message("Bm is 0, imputed and observed values of bin. variable are the same.")
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
#' @param qhats vector of Qhats(means of response for each imputed dataset)
#' @param m number of imputed datasets
#' @param nrow number of observations in the imputed dataset
#'
#' @return degrees of freedom
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' qhats = Qhats(imp, "hyp")
#' m = imp$m
#' nrow = imp$data %>% nrow()
#' dof(qhats, m, nrow)
#'
dof <- function(qhats, m, nrow) {
  rm = Rm(qhats, m, nrow)

  return((m - 1)*(1 + 1/rm)^2)
}

#' Calculates the specified Wilson CI of a binomial proportion
#' variable, given imputed data sets.
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable (must be 0-1 valued)
#' @param ci_level desired confidence interval level (defaults to 95%)
#' @param summaries boolean: should summary helper values be printed (default TRUE)
#'
#'
#' @return two-length vector of Wilson lower CI and upper CI
#' @export
#' @importFrom dplyr mutate
#'
#' @examples
#' imp = mice::mice(mice::nhanes %>% dplyr::mutate(hyp = hyp-1))
#' mi_wilson(imp, "hyp", 0.95)
#'
mi_wilson <- function(mids_obj=NA, response=NA, ci_level=0.95, summaries=TRUE) {

  #if confidence interval is invalid
  if(ci_level<=0 | ci_level>= 1) {
    stop("CI level must be between 0 and 1.")
  }

  #test if provided mids object and response variable is valid
  resp = tryCatch(mids_obj$data %>% select(all_of(response)),
                  error = function(e)
                    stop("Invalid mids object and/or response variable."))

  #if response is not 0-1 valued
  if(all(lapply(resp,
                function(x) x %in% c(NA,0,1)) %>% unlist())==FALSE) {
    stop(paste(response,"must be 0-1 binary encoded."))
  }

  qhats = Qhats(mids_obj, response)
  m = mids_obj$m
  nrow = mids_obj$data %>% nrow()


  qbar = Qbar(qhats)
  rm = Rm(qhats, m, nrow)

  #if imputed binomial proportions are the same
  if(rm==0) {
    df = Inf
    warning("Imputed binomial proportions are identical; degrees of freedom set to infinity.")
  }
  else {
    df = dof(qhats, m, nrow)
  }


  #print summaries if desired
  if(summaries) {
    print(paste("Qbar: ", qbar))
    print(paste("Rm: ", rm))
    print(paste("dof: ",df))
  }

  t_score = stats::qt(ci_level, df)

  tquad = t_score^2 / nrow + t_score^2 * rm /nrow
  center = (2 * qbar + tquad)/(2*(1 + tquad))
  half_width = sqrt( (2*qbar + tquad)^2 / (4*(1+tquad)^2) -
                       qbar^2 / (1+tquad) )

  return(c(center - half_width, center + half_width))
}

#' Calculates the specified Wald CI of a binomial proportion
#' variable, given imputed data sets.
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable (must be 0-1 valued)
#' @param ci_level desired confidence interval level (defaults to 95%)
#' @param summaries boolean: should summary helper values be printed (default TRUE)
#'
#' @return two-length vector of Wald lower CI and upper CI
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' imp = mice::mice(mice::nhanes %>% dplyr::mutate(hyp = hyp-1))
#' mi_wald(imp, "hyp", 0.95)
#'
mi_wald <- function(mids_obj=NA, response=NA, ci_level=0.95, summaries=TRUE) {

  #test if provided mids object and response variable is valid
  resp = tryCatch(mids_obj$data %>% select(all_of(response)),
                  error = function(e)
                    stop("Invalid mids object and/or response variable."))


  #if CI not between 0 and 1
  if(ci_level<=0 | ci_level>= 1) {
    stop("CI level must be between 0 and 1.")
  }

  #if response is not 0-1 valued
  if(all(lapply(resp,
                function(x) x %in% c(NA,0,1)) %>% unlist())==FALSE) {
    stop(paste(response,"must be 0-1 binary encoded."))
  }

  qhats = Qhats(mids_obj, response)
  m = mids_obj$m
  nrow = mids_obj$data %>% nrow()

  qbar = Qbar(qhats)
  tm = Tm(qhats, m, nrow)
  df = dof(qhats, m, nrow)
  t_score = stats::qt(ci_level, df)

  #print summaries if desired
  if(summaries) {
    print(paste("Qbar: ", qbar))
    print(paste("Tm: ", tm))
    print(paste("dof: ",df))
  }

  return(c(qbar - sqrt(tm)*t_score, qbar + sqrt(tm)*t_score))
}


#' Calculates the MI-Wilson interval if given a vector of observed
#' binomial proportions (one for each imputed data frame)
#'
#' @param phats vector of binomial proportions (one for each imputation)
#' @param n the common number of observations over the imputed dataframes
#' @param ci_level desired confidence interval level (default 95%)
#' @param summaries boolean: should summary helper values be printed (default TRUE)
#'
#' @return two-length vector of Wilson lower CI and upper CI
#' @export
#'
#' @examples
#' phats = c(0.2, 0.23, 0.25)
#' mi_wilson_phat(phats, 100, 0.99, TRUE)
#'
mi_wilson_phat <- function(phats, n, ci_level =0.95, summaries = TRUE) {

  #if CI not between 0 and 1
  if(ci_level<=0 | ci_level>= 1) {
    stop("CI level must be between 0 and 1.")
  }

  #if n not specified
  if (missing(n)) stop("Must specify number of observations as function parameter.")

  qhats = phats
  m = length(phats)
  nrow = n

  qbar = Qbar(qhats)
  rm = Rm(qhats, m, nrow)

  #if imputed binomial proportions are the same
  if(rm==0) {
    df = Inf
    warning("Imputed binomial proportions are identical; degrees of freedom set to infinity.")
  }
  else {
    df = dof(qhats, m, nrow)
  }


  #print summaries if desired
  if(summaries) {
    print(paste("Qbar: ", qbar))
    print(paste("Rm: ", rm))
    print(paste("dof: ",df))
  }

  t_score = stats::qt(ci_level, df)

  tquad = t_score^2 / nrow + t_score^2 * rm /nrow
  center = (2 * qbar + tquad)/(2*(1 + tquad))
  half_width = sqrt( (2*qbar + tquad)^2 / (4*(1+tquad)^2) -
                       qbar^2 / (1+tquad) )

  return(c(center - half_width, center + half_width))
}

#' Calculates the MI-Wald interval if given a vector of observed
#' binomial proportions (one for each imputed data frame)
#'
#' @param phats vector of binomial proportions (one for each imputation)
#' @param n the common number of observations over the imputed dataframes
#' @param ci_level desired confidence interval level (default 95%)
#' @param summaries boolean: should summary helper values be printed (default TRUE)
#'
#' @return two-length vector of Wilson lower CI and upper CI
#' @export
#'
#' @examples
#' phats = c(0.2, 0.23, 0.25)
#' mi_wald_phat(phats, 100, 0.99, TRUE)
#'
mi_wald_phat <- function(phats, n, ci_level = 0.95, summaries = TRUE) {

  #if CI not between 0 and 1
  if(ci_level<=0 | ci_level>= 1) {
    stop("CI level must be between 0 and 1.")
  }

  #if n not specified
  if (!is.numeric(n)) stop("Must specify number of observations as function parameter.")

  qhats = phats
  m = length(phats)
  nrow = n

  qbar = Qbar(qhats)
  tm = Tm(qhats, m, nrow)
  df = dof(qhats, m, nrow)
  t_score = stats::qt(ci_level, df)

  #print summaries if desired
  if(summaries) {
    print(paste("Qbar: ", qbar))
    print(paste("Tm: ", tm))
    print(paste("dof: ",df))
  }

  return(c(qbar - sqrt(tm)*t_score, qbar + sqrt(tm)*t_score))
}
