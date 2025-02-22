# varGuid: an R package Implementing Variance Guided Regression Models for Heteroscedastic Data

## Authors
Sibei Liu (sxl4188@miami.edu) and Min Lu (m.lu6@umiami.edu)

## Reference
Liu S. and Lu M. Variance guided regression models for heteroscedastic data (under review)

#### Description
Advanced regression techniques to address heteroscedasticity in linear models. It features two key algorithms: an iteratively reweighted least squares (IRLS) method for robust coefficient estimation under linear variance patterns, and a biconvex algorithm that creates artificial grouping effects to capture nonlinear relationships. 

### Installation
#### Step 1: download and install the "cvxclustr" package using "cvxclustr_1.1.0.tar.gz" from the root folder

```
install.packages("/YourFolder/cvxclustr_1.1.0.tar.gz", repos = NULL, type = "source")
```
#### Step 2: install and library the "varGuid" package
```
## install.packages("devtools") ## install devtools if not already installed
devtools::install_github("Sibeiliunew/varGuid")
library(varGuid)
```
### Examples

* Step 1:  Obtain linear coefficients using the IRLS algorithm:
```
data(tnbc)
obj <- rforest(subtype~., data = tnbc[1:100,c(1:5,337)])
importance(obj)
predict(obj)$label
predict(obj, tnbc[101:110,1:5])$label

### pair() to convert continuous variables to binary ranked pairs
tnbc[101:110,1:5]
datp <- pair(tnbc[101:110,1:5])
datp
predict(obj, datp, newdata.pair = TRUE)$label
```

* Step 2: Create artificial grouping effect for nonlinear prediction:
```
objr <- extract.rules(obj)
objr$rule[1:5,]
predict(objr)$label[1:5]

objrs <- select.rules(objr,tnbc[110:130,c(1:5,337)])
predict(objrs, tnbc[111:120,1:5])$label
objrs$rule[1:5,]
```

* Note that Step 2 does not change the coefficient estimates and may be skipped unless outcome prediction is required.

  
