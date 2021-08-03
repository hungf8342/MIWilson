library(mice)
library(dplyr)

test_that("MI-Wilson works", {
  nhanes = mice::nhanes %>%
    mutate(hyp = hyp-1)
  set.seed(47)
  imp = mice(nhanes)
  qhats = Qhats(imp, "hyp")
  m = imp$m

  expect_equal(Qbar(qhats), 0.24)
  expect_equal(Bm(qhats, m), 0.0016)
  expect_equal(mi_wilson(imp, "hyp", 0.95)[1]>0, TRUE)
})

test_that("Correct errors thrown", {
  nhanes_xtreme = mice::nhanes %>% mutate(hyp=ifelse(is.na(hyp),NA,0))
  imp_x = mice::mice(nhanes_xtreme)

  imp = mice::mice(mice::nhanes)

  expect_error(mi_wald(imp_x, "hyp", 1.99),"^CI.*")
  expect_error(mi_wilson(imp_x,"hyp"),".*unable to impute.*")
  expect_error(mi_wilson(imp, "hyp"), ".*binary encoded.*")
  expect_warning(mi_wilson_phat(phats = rep(0.3,3), n = 10), ".*degrees of freedom.*")
})
