# R functions used in Rincent *et al.* 2018

This repository includes R functions for evaluating the predictive abilities of SNP and NIRS through cross-validations as reported by Rincent *et al.* 2018

The functions can be accessed within R as follows:

```
source("https://raw.githubusercontent.com/visegura/PS/master/rfuncPS.r")
```

**Dependencies**

The following r packages are required: `doParallel`, `MASS`, `emma`, `corpcor` and `BGLR`.

Please note that emma is not the one available from CRAN, it should be downloaded from [UCLA](http://mouse.cs.ucla.edu/emma/) or for Windows users from [here](https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz).

**Use**

Please see the vignette [here](https://github.com/visegura/PS/blob/master/PS_R_functions_github.md) or [there in pdf](https://github.com/visegura/PS/blob/master/PS_R_functions.pdf)
