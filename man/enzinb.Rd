\name{ezinb}
\alias{ezinb}
\title{The expected value of the censored zero-inflated negative binomial model}
\description{
   'ezinb' computes the expected value of the censored zero-inflated 
    negative binomial regression model
}
\usage{
ezinb((coeff = NULL, resp = y ~ ., upper = r ~ ., x =  ~ 1, z =  ~ 1, data=NULL,
	sub = NULL)
}

\arguments{
   \item{coeff}{Regression coefficients from 'zicensor' 
	a.k.a censored ZINB regression model.}
   \item{resp}{a symbolic description of the response, of the form
   	resp=y~., where y is the response variable.}
   \item{x}{a formula to describe the covariates to the Poisson (or Negative binomial)
    part of the model, of the form x=~x1 + x2 + ..., where x1,x2,... 
    are explanatory variables.}  
   \item{z}{a formula to describe the covariates to zero-indlated part, of the form
   	z=~z1 + z2 + ..., where z1, z2,... are explanatory variables.}  
   \item{data}{a frame containing the variables in the model.}
   \item{sub}{an optional logical expression specifying a subset of observations to be
           used in the fitting process. This should be specified in the form: e.g. 
           subsets =expression(w==0).}
 }
\value{
 returns a vector of the censored zinb expected values. 
 }
  

\seealso{
\code{\link{zicounts},\link{zicensor}}.
}

\examples{
## load the artificial 'teeth' data
 data(teeth)
 names(teeth)
 
##  fit zero-inflated negative binomial regression model
zinb.cens <- zicensor(resp = dmft~., upper = decid~.,x =~gender + age,z =~gender + age, data=teeth)
zinb.cens 

## expected value of the censored zinb model
## uncomment to run
#zinb.e <- ezinb(coeff=zinb.cens$coeff,resp = dmft~., upper = decid~.,
#	x =~gender + age,z =~gender + age, data=teeth)
#zinb.e
}

\keyword{models}
\author{S. M. Mwalili <samuelmwalili@yahoo.com>}
