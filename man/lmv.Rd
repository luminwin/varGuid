\name{lmv}
\alias{lmv}
\alias{lmv}
\title{
Coefficient Estimation from the VarGuid Reweighted Least Squares or Lasso Method
}
\description{
The function implements the iteratively reweighted least squares (IRLS) or iteratively reweighted Lasso method for robust coefficient estimation under linear variance patterns.   The inference of coefficient estimator was directly imported from the \pkg{sandwich}
  package.  Use the command \code{package?sandwich} for more information.
}
\usage{
lmv(X, Y, M =  10, step = 1, tol = exp(-10), lasso = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{X}{input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format (inherit from class \code{"sparseMatrix"} as in package \code{Matrix}). Requirement: nvars >1; in other words, x should have 2 or more columns. }

  \item{Y}{response variable. Currently only for \code{family="gaussian"}.}

  \item{M}{maximum number of iterations.}

  \item{step}{ scale parameter for the data weight. }

  \item{tol}{tolerance parameter for convergence.}

  \item{lasso}{should the iteratively reweighted Lasso method be used? IRLS will be used if set to \code{FALSE}.}
}
\value{
An object of class \code{(varGuid)} with the following
  components:
  \item{beta}{Coefficient estimator from the VarGuid method.}
  \item{obj.OLS}{Object from the OLS method using the \code{lm()} without weights. This object is used to initialize  the VarGuid method with IRLS.}
  \item{obj.lasso}{Object from the Lasso method using the \code{glmnet()} without weights. This object is used to initialize  the VarGuid method with iteratively reweighted Lasso.}
  \item{obj.varGuid}{Object from the last iteration of the VarGuid model from either \code{lm()} or \code{glmnet()}, depending if argument \code{lasso = FALSE} or \code{lasso = TRUE}.}
  \item{res}{Residual from the VarGuid model.}
  \item{obj.varGuid.coef}{Summary tables for the coefficient estimators. Inferences are obtained from the \code{vcovHC()} function from the \code{sandwich} package. }
  \item{X}{The same as argument \code{X}.}

}
\references{
Liu S. and Lu M. (2025) Variance guided regression models for heteroscedastic data (under review)
}
\author{
Sibei Liu (Maintainer,<sxl4188@miami.edu>) and Min Lu
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

########### Obtain linear coefficients using the IRLS algorithm:

o <- lmv(X = train[,-yid] , Y = train[,yid], lasso = FALSE) 
o$obj.varGuid.coef$HC3 ## coefficient estimator from VarGuid regression
summary(o$obj.OLS) ## coefficient estimator from OLS regression

 o <- lmv(X = train[,-yid] , Y = train[,yid], lasso = TRUE) 
 o$beta ## coefficient estimator from VarGuid-Lasso regression
 o$obj.lasso$beta ## coefficient estimator from Lasso regression
}
}
