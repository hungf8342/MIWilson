#' Calculate Qbar (average response over MICE datasets)
#'
#' @param response string name of endogenous variable
#' @param mids_obj mids object created by mice package
#'
#' @return Qbar: the average response over MICEd datasets.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Qbar(imp, "bmi")
#'
Qbar <- function(mids_obj,response) {
  qbar = lapply(1:mids_obj$m,
             function(i) mice::complete(mids_obj, i) %>%
               dplyr::select(response) %>%
               colMeans()) %>%
    unlist() %>% mean()

  return(qbar)
}

#' Calculate Ubar (average response variance over MICE datasets)
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#'
#' @return Ubar: average response variance over MICE datasets
#' @export
#'
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Ubar(imp, "bmi")
#'
Ubar <- function(mids_obj, response) {
  ubar = lapply(1:mids_obj$m,
                function(i) mice::complete(mids_obj, i) %>%
                  dplyr::select(response) %>%
                  unlist() %>% stats::var()) %>%
    unlist() %>% mean()

  return(ubar)
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
#' Bm(imp, "bmi")
#'
Bm <- function(mids_obj, response) {
  qbar = Qbar(mids_obj, response)
  bm = lapply(1:mids_obj$m,
              function(i) ((mice::complete(mids_obj,i) %>%
                              dplyr::select(response) %>%
                              unlist() %>% mean())-qbar)^2) %>%
    unlist() %>% sum()

  return(bm/(mids_obj$m-1))
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
#' Tm(imp, "bmi")
#'
Tm <- function(mids_obj, response) {
  ubar = Ubar(mids_obj, response)
  bm = Bm(mids_obj, response)

  return((1 + 1/mids_obj$m)*bm + ubar)
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
#' dof(imp, "bmi")
#'
dof <- function(mids_obj, response) {
  ubar = Ubar(mids_obj, response)
  bm = Bm(mids_obj, response)
  m = mids_obj$m
  rm = (1 + 1/m)*bm/ubar

  return((m - 1)*(1 + 1/rm)^2)
}

#' Calculates the specified Wilson CI of a binomial proportion
#' variable, given imputed data sets.
#'
#' @param mids_obj mids object created by mice package
#' @param response string name of response variable
#' @param ci_level desired confidence interval level
#'
#' @return two-length vector of lower CI and upper CI
#' @export
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' mi_wilson(imp, "bmi", 0.95)
#'
mi_wilson <- function(mids_obj, response, ci_level) {
  qbar = Qbar(mids_obj, response)
  tm = Tm(mids_obj, response)
  df = dof(mids_obj, response)
  t_score = stats::qt(ci_level, df)

  return(c(qbar - sqrt(tm)*t_score, qbar + sqrt(tm)*t_score))
}
