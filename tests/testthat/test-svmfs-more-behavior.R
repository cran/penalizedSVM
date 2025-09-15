test_that("svmfs(DrHSVM) returns coherent model and predictions", {
  skip_on_cran()
  seed <- 123
  n=200
  train<-sim.data(n = n, ng = 100, nsg = 10, corr=FALSE, seed=seed )
  print(str(train)) 
  test<-sim.data(n = n, ng = 100, nsg = 10, corr=FALSE, seed=seed+1 )
  print(str(test)) 
  bounds=t(data.frame(log2lambda1=c(-10, 10)))
  colnames(bounds)<-c("lower", "upper")	
  
  # computation intensive; for demostration reasons only for the first 100 features 
  # and only for 10 Iterations maxIter=10, default maxIter=700
  print("start interval search")
  suppressWarnings( fit<- svmfs(t(train$x)[,1:100], y=train$y,
                            fs.method="scad", bounds=bounds, 
                            cross.outer= 0, grid.search = "interval",  maxIter = 10, 
                            inner.val.method = "cv", cross.inner= 5, maxevals=500,
                            seed=seed, parms.coding = "log2", show="none", verbose=FALSE ) )
  
  expect_true(inherits(fit, "penSVM"))
  f <- fit$model
  expect_true(is.numeric(f$w))
  expect_true(is.numeric(f$b) && length(f$b) == 1)
  expect_true(is.numeric(f$xind) || is.integer(f$xind))
  expect_true(length(f$w) == length(f$xind))
  expect_true(all(f$xind >= 1 & f$xind <= ncol(test$x)))
  
  pr<-predict.penSVM(fit, t(test$x)[,1:100], newdata.labels=test$y)

  expect_type(pr, "list")
  expect_equal(length(pr$pred.class), n)
  expect_equal(length(as.vector(pr$fitted)), n)
})

test_that("svmfs(scad+L2) returns coherent model and predictions", {
  skip_on_cran()
  set.seed(33)
  n <- 40; p <- 12
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("G", seq_len(p))
  y <- rep(c(-1, 1), each = n/2)
  x[y == -1, 1] <- x[y == -1, 1] - 3.0
  x[y ==  1, 1] <- x[y ==  1, 1] + 3.0
  x[y == -1, 3] <- x[y == -1, 3] - 1.2
  x[y ==  1, 3] <- x[y ==  1, 3] + 1.2

  suppressWarnings(fit <- svmfs(
    x, y,
    fs.method = "scad+L2",
    grid.search = "interval",
    inner.val.method = "cv",
    cross.inner = 3,
    maxevals = 12,
    show = "none",
    seed = 33
  ))
  expect_true(inherits(fit, "penSVM"))
  f <- fit$model
  expect_true(is.numeric(f$w))
  expect_true(is.numeric(f$b) && length(f$b) == 1)
  expect_true(is.numeric(f$xind) || is.integer(f$xind))
  expect_true(length(f$w) == length(f$xind))
  expect_true(all(f$xind >= 1 & f$xind <= ncol(x)))

  pr <- predict(fit, x, newdata.labels = y)
  expect_type(pr, "list")
  expect_equal(length(pr$pred.class), n)
  expect_equal(length(as.vector(pr$fitted)), n)
})
