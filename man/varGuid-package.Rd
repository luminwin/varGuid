\name{varGuid}
\alias{varGuid}
\title{
  Variance Guided Regression Models for Heteroscedastic Data
}
\description{
The package \code{varGuid} implements advanced regression techniques to address heteroscedasticity in linear models. It features two key algorithms: an iteratively reweighted least squares (IRLS) method for robust coefficient estimation under linear variance patterns, and a biconvex algorithm that creates artificial grouping effects to capture nonlinear relationships. 
}
\author{
Sibei Liu (Maintainer,<sxl4188@miami.edu>) and Min Lu
}
\references{
Liu S. and Lu M. (2025) Variance guided regression models for heteroscedastic data (under review)
}
\examples{
\donttest{
library(varGuid)
data(cobra2d, package = "varGuid")
dat <- cobra2d
tid <- sample(1:nrow(dat), 200)
train <- dat[-tid,]
test <- dat[tid,]
yid <- which(colnames(dat) == "y")

########### Step 1:  Obtain linear coefficients using the IRLS algorithm:

o <- lmv(X = train[,-yid] , Y = train[,yid], lasso = FALSE) 
o$obj.varGuid.coef$HC3 ## coefficient estimator from VarGuid regression
summary(o$obj.OLS) ## coefficient estimator from OLS regression

# o <- lmv(X = train[,-yid] , Y = train[,yid], lasso = TRUE) 
# o$beta ## coefficient estimator from VarGuid-Lasso regression
# o$obj.lasso$beta ## coefficient estimator from Lasso regression

########### Step 2: Create artificial grouping effect for nonlinear prediction:

# create artificial grouping effects
y.obj <- ymodv(obj = o) 

# outcome prediction on new data
pred <- predict.varGuid(mod=y.obj,lmvo = o,newdata = test[,-yid]) 

# RMSE
sqrt(colMeans((  matrix(replicate(ncol(pred),test[,yid]),ncol=ncol(pred))-pred)^2, na.rm = TRUE)) 

########### Note that Step 2 does not change the coefficient estimates and 
########### may be skipped unless outcome prediction is required.
}
}
