test_that("svmfs(scad) returns a valid penSVM model and predict() works (tiny synthetic)", {
  skip_on_cran()
  set.seed(42)
  n <- 40; p <- 12
  # strong signal on first 2 features
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", seq_len(p))
  y <- rep(c(-1, 1), each = n/2)
  x[y == -1, 1] <- x[y == -1, 1] - 2.5
  x[y ==  1, 1] <- x[y ==  1, 1] + 2.5
  x[y == -1, 2] <- x[y == -1, 2] - 2
  x[y ==  1, 2] <- x[y ==  1, 2] + 2

  fit <- suppressWarnings(svmfs(
    x, y,
    fs.method = "scad",
    grid.search = "interval",
    inner.val.method = "cv",
    cross.inner = 3,
    maxevals = 12,
    show = "none",
    seed = 1
  ))
  expect_true(inherits(fit, "penSVM"))
  expect_true(is.list(fit$model))
  # model internals present
  f <- fit$model
  expect_true(is.numeric(f$w))
  expect_true(is.numeric(f$b) && length(f$b) == 1)
  expect_true(is.integer(f$xind) || is.numeric(f$xind))
  expect_true(length(f$w) == length(f$xind))
  expect_true(length(f$w) > 0)
  # indices within bounds and names consistent
  expect_true(all(f$xind >= 1 & f$xind <= ncol(x)))
  expect_true(all(names(f$w) %in% colnames(x)))

  # predict() on training data
  pr <- predict(fit, x, newdata.labels = y)
  expect_type(pr, "list")
  expect_equal(length(pr$pred.class), n)
  expect_true(all(levels(pr$pred.class) == c("-1", "1")))
  # sep (fitted) length
  expect_equal(length(as.vector(pr$fitted)), n)
  # not asserting tight accuracy to keep runtime small / stable, but sanity check:
  expect_true(is.numeric(pr$error) || is.na(pr$error))
})

test_that("svmfs(1norm) produces a model with coherent shapes and print() works", {
  skip_on_cran()
  set.seed(7)
  n <- 30; p <- 10
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("X", seq_len(p))
  y <- rep(c(-1, 1), each = n/2)
  x[y == -1, 1] <- x[y == -1, 1] - 3
  x[y ==  1, 1] <- x[y ==  1, 1] + 3

  fit <- svmfs(
    x, y,
    fs.method = "1norm",
    grid.search = "interval",
    inner.val.method = "cv",
    cross.inner = 3,
    maxevals = 10,
    show = "none",
    seed = 2
  )
  expect_true(inherits(fit, "penSVM"))
  f <- fit$model
  expect_true(is.numeric(f$w) && length(f$w) > 0)
  expect_true(is.numeric(f$b) && length(f$b) == 1)
  expect_true(all(f$xind >= 1 & f$xind <= ncol(x)))
  expect_equal(length(f$w), length(f$xind))

  # S3 print should not error
  expect_output(print(fit))

  # predict on train
  pr <- predict(fit, x, newdata.labels = y)
  expect_equal(length(pr$pred.class), n)
  expect_true(all(levels(pr$pred.class) == c("-1", "1")))
})
