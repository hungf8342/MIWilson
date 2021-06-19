context("Tests whether the MI-Wilson function, given a mice object,
        outputs correct CI")
library(mice)
library(dplyr)

test_that("MI-Wilson works", {
  nhanes = mice::nhanes %>%
    mutate(hyp = hyp-1)
  set.seed(47)
  imp = mice(nhanes)

  expect_equal(Qbar(imp,"hyp"), 0.24)
  expect_equal(Bm(imp, "hyp"), 0.0016)
  expect_equal(mi_wilson(imp, "hyp", 0.95)[1]>0, TRUE)
})


