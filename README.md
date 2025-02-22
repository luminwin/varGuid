# varGuid: an R package Implementing Variance Guided Regression Models for Heteroscedastic Data

### Authors
Sibei Liu (sxl4188@miami.edu) and Min Lu (m.lu6@umiami.edu)

### Reference
Liu S. and Lu M. (2025) Variance guided regression models for heteroscedastic data (under review)

### Description
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
data(cobra2d, package = "varGuid")
dat <- cobra2d
tid <- sample(1:nrow(dat), 200)
train <- dat[-tid,]
test <- dat[tid,]
yid <- which(colnames(dat) == "y")

o <- lmv(X = train[,-yid] , Y = train[,yid], lasso = FALSE) 
o$obj.varGuid.coef$HC3 ## coefficient estimator from VarGuid regression
summary(o$obj.OLS) ## coefficient estimator from OLS regression

o2 <- lmv(X = train[,-yid] , Y = train[,yid], lasso = TRUE) 
o2$beta ## coefficient estimator from VarGuid-Lasso regression
o2$obj.lasso$beta ## coefficient estimator from Lasso regression
```

* Step 2: Create artificial grouping effect for nonlinear prediction:
```
# create artificial grouping effects
y.obj <- ymodv(obj = o) 

# outcome prediction on new data
pred <- predict.varGuid(mod=y.obj,lmvo = o,newdata = test[,-yid]) 

# RMSE
sqrt(colMeans((  matrix(replicate(ncol(pred),test[,yid]),ncol=ncol(pred))-pred)^2, na.rm = TRUE)) 

```

* Note that Step 2 does not change the coefficient estimates and may be skipped unless outcome prediction is required.

  
