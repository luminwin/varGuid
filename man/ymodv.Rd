\name{ymodv}
\alias{ymodv}
\title{
Convex Clustering for Creating Artificial Grouping Effects to Capture Nonlinear Relationships
}
\description{
The function implements convex clustering for creating artificial grouping effects for the object from the \code{lmv()} function.   The convex clustering algorithm was directly imported from the \pkg{cvxclustr}
  package.  Use the command \code{package?cvxclustr} for more information.
}
\usage{
ymodv(obj, nu = c(seq(0,9, length.out=5)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{obj}{Object from the \code{lmv()} function. }

  \item{nu}{Parameter to control the number of artificial groups. This is the same as the \code{gamma} argument from the \code{cvxclust()} function from the \pkg{cvxclustr}
  package. }

 
}

\value{
An object of class \code{(varGuid)} with the following
  components:
  \item{clust}{A list in which \code{sol} contains the output from the \code{cvxclust()} function, along with other parameters used as arguments for \code{cvxclust()}.}
  \item{data}{A dataframe from the clustering algorithm, where \code{z} contains the grouping factor.}
  \item{cent}{Centers of the artificial groups.}

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
y.obj <- ymodv(obj = o, nu = c(seq(90,190, length.out=5))) 

# outcome prediction on new data
pred <- predict.varGuid(mod=y.obj,lmvo = o,newdata = test[,-yid]) 

# RMSE
sqrt(colMeans((  matrix(replicate(ncol(pred),test[,yid]),ncol=ncol(pred))-pred)^2, na.rm = TRUE)) 
}
}
