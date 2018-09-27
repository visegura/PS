# R functions used in Rincent *et al.* for evaluating the predictive abilities of SNP and NIRS through cross-validations

## Setup

### PS functions


```r
source("https://raw.githubusercontent.com/visegura/PS/master/rfuncPS.r")
```

### Dependencies

Make use of the very cool `anyLib` package to check and install all dependencies but `emma`:


```r
install.packages("anyLib")
```

```
## Installing package into '/usr/local/lib/R/site-library'
## (as 'lib' is unspecified)
```

```r
library(anyLib)
anyLib(c("doParallel", "MASS", "corpcor", "BGLR"))
```

```
## Loading required package: doParallel
```

```
## Loading required package: foreach
```

```
## Loading required package: iterators
```

```
## Loading required package: parallel
```

```
## Loading required package: MASS
```

```
## Loading required package: corpcor
```

```
## Loading required package: BGLR
```

```
## doParallel       MASS    corpcor       BGLR 
##       TRUE       TRUE       TRUE       TRUE
```

`emma` package is different from the one available on CRAN, it must be installed from github as follows:


```r
install.packages("https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz",
                 repos = NULL)
```

```
## Installing package into '/usr/local/lib/R/site-library'
## (as 'lib' is unspecified)
```

```r
anyLib("emma")
```

```
## Loading required package: emma
```

```
## emma 
## TRUE
```

## Datasets

The example is carried out on the poplar dataset used in the paper for the trait "Bud Set" evaluated in Orléans trial.


```r
print('test')
```

```
## [1] "test"
```


