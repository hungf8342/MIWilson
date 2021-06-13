context("Tests whether the MI-Wilson function, given a mice object,
        outputs correct CI")

test_that("multiplication works", {
  nhanes = mice::nhanes %>%
    dplyr::mutate(hyp = hyp-1)
  set.seed(47)
  imp = mice::mice(nhanes)

  expect_equal(Qbar(imp,"hyp"), 0.24)
})
