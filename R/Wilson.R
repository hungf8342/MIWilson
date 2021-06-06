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
#' @param response string name of endogenous variable
#'
#' @return Ubar: average response variance over MICE datasets
#' @export
#'
#'
#' @examples
#' imp = mice::mice(mice::nhanes)
#' Ubar(imp, "bmi")
Ubar <- function(mids_obj, response) {
  ubar = lapply(1:mids_obj$m,
                function(i) mice::complete(mids_obj, i) %>%
                  dplyr::select(response) %>%
                  unlist() %>% stats::var()) %>%
    unlist() %>% mean()

  return(ubar)
}
