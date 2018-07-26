`print.penSVM` <-
  function(x,...){
    cat("\nBias = ", x$b)
    cat("\nSelected Variables= ", names(x$w))
    cat("\nCoefficients:\n  ")
    print(x$w)
    cat("\n")
}

