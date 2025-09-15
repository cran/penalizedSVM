---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# penalizedSVM

<https://doi.org/10.32614/CRAN.package.penalizedSVM>

<!-- badges: start -->
[![DOI](https://img.shields.io/badge/doi-10.32614/CRAN.package.penalizedSVM-blue.svg)](https://doi.org/10.32614/CRAN.package.penalizedSVM)
[![CRAN status](https://www.r-pkg.org/badges/version/penalizedSVM)](https://CRAN.R-project.org/package=penalizedSVM)
<!-- badges: end -->

The goal of penalizedSVM is to support Vector Machine (SVM) classification with simultaneous feature selection using penalty functions is implemented. The smoothly clipped absolute deviation (SCAD), 'L1-norm', 'Elastic Net' ('L1-norm' and 'L2-norm') and 'Elastic SCAD' (SCAD and 'L2-norm') penalties are available. The tuning parameters can be found using either a fixed grid or a interval search.


## Installation

You can install the released version of peperr from [CRAN](https://CRAN.R-project.org) with:


``` r
install.packages("penalizedSVM")
```

And the development version from [GitHub](https://github.com/) with:


``` r
install.packages("devtools")
devtools::install_github("fbertran/penalizedSVM")
```

## Example

This is a basic example.


``` r
library(penalizedSVM)
```


``` r
seed<- 123
train<-sim.data(n = 200, ng = 100, nsg = 10, corr=FALSE, seed=seed )
print(str(train))
#> List of 3
#>  $ x   : num [1:100, 1:200] 0.547 0.635 -0.894 0.786 2.028 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:100] "pos1" "pos2" "pos3" "pos4" ...
#>   .. ..$ : chr [1:200] "1" "2" "3" "4" ...
#>  $ y   : Named num [1:200] 1 -1 -1 1 -1 -1 1 -1 -1 -1 ...
#>   ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
#>  $ seed: num 123
#> NULL
```

Train standard svm

``` r
my.svm<-svm(x=t(train$x), y=train$y, kernel="linear")
```

Test with other data

``` r
test<- sim.data(n = 200, ng = 100, nsg = 10, seed=(seed+1) )
```

Check accuracy standard SVM

``` r
my.pred <-ifelse( predict(my.svm, t(test$x)) >0,1,-1)
```

Check accuracy:

``` r
table(my.pred, test$y)
#>        
#> my.pred -1  1
#>      -1 81 23
#>      1  30 66
```
