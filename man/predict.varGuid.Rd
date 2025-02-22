\name{predict.varGuid}
\alias{predict.varGuid}
\alias{predict}
\title{
Variance Guided Prediction from Artificial Grouping Effects to Capture Nonlinear Relationships
}
\description{
The function makes predictions for new data for the objects from the \code{ymodv()} function and the \code{lmv()} function.
}
\usage{
\method{predict}{varGuid}(mod,lmvo,newdata, redDim = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 \item{mod}{Object from the \code{ymodv()} function. }
  \item{lmvo}{Object from the \code{lmv()} function. }

  \item{newdata}{ 	A data frame or matrix in which to look for variables with which to predict. }

   \item{redDim}{Set \code{redDim = TRUE} if only variables with non-zero coefficients from the \code{lmv()} function were used for artificial grouping effects. }
}

\value{
An object of class \code{(varGuid)} with the following
  components:
  \item{VarGuidOriginal}{Predicted outcome from the \code{lmv()} object without adding rtificial grouping effects.}
  \item{VarGuid}{Predicted outcome with artificial grouping effects.}
  \item{OLS}{Predicted outcome from the OLS method using the \code{lm()} function without weights. This model is used to initialize the VarGuid method with IRLS and corresponds to the \code{obj.OLS} output in the \code{lmv()} object. This is \code{NULL} if \code{lasso = TRUE}.}
    \item{Lasso}{Predicted outcome from the Lasso method using the \code{glmnet()} function without weights. This model is used to initialize the VarGuid method with iteratively reweighted Lasso and corresponds to the \code{obj.lasso} output in the \code{lmv()} object. This is \code{NULL} if \code{lasso = FALSE}.}

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

########### Step 2: Create artificial grouping effect for nonlinear prediction:

# create artificial grouping effects
y.obj <- ymodv(obj = o) 

# outcome prediction on new data
pred <- predict.varGuid(mod=y.obj,lmvo = o,newdata = test[,-yid]) 

# RMSE
sqrt(colMeans((  matrix(replicate(ncol(pred),test[,yid]),ncol=ncol(pred))-pred)^2, na.rm = TRUE)) 
}
}
