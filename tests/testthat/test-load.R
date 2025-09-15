test_that("package namespace is loadable", {
  expect_true(requireNamespace("penalizedSVM", quietly = TRUE))
  expect_true(isNamespaceLoaded("penalizedSVM"))
})
