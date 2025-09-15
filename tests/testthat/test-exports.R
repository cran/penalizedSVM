test_that("exported objects exist (if any)", {
  ns <- asNamespace("penalizedSVM")
  exports <- getNamespaceExports(ns)
  if (length(exports) == 0L) skip("No exports in NAMESPACE.")
  for (x in exports) {
    ok <- exists(x, envir = ns, inherits = FALSE) || methods::isGeneric(x)
    expect_true(ok, info = x)
  }
})
