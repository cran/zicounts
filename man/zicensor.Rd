\name{zicensor}
\alias{zicensor}
\title{Fitting classical and zero-inflated count regression models}
\description{
   'zicensor' is used to fit the classical and the zero-inflated
    interval censored count data regression model. In particular, the censored Poisson,
    censored negative binomial (NB),censored  zero-inflated Poisson (ZIP) and 
    censored zero-inflated negative binomial (ZINB) are fitted.
    Type help(zicounts) for count data models.
}
\usage{
zicensor(parm = NULL, resp = y ~ ., upper = r ~ ., x = ~1, z = ~1, data,
	distrname = "ZINB", sub = NULL, ntries = 5, method = "BFGS", ...)
}

\arguments{
  \item{resp}{a symbolic description of the response, of the form
   	resp=y~., where y is the response variable.}
   \item{upper}{a symbolic description of the response's upper bound, of the form
      	upper=r~., where r is the upper bound variable.}
   \item{x}{a formula to describe the covariates to the Poisson (or Negative binomial)
    part of the model, of the form x=~x1 + x2 + ..., where x1,x2,... 
    are explanatory variables.}  
   \item{z}{a formula to describe the covariates to zero-indlated part, of the form
   	z=~z1 + z2 + ..., where z1, z2,... are explanatory variables.}  
   \item{distrname}{a description of the error distribution to
	be used in the model. The default is "ZINB" (zero inflated negative binomial).
	Other choices include "Poisson", "ZIP" (zero-inflated Poisson) and 
	"NB" (negative binomial).}
   \item{parm}{starting values for all the parameters in the model.
    Number of parameters: Poisson (number of x variables), negative binomial 
    (number of x variables + 1), ZIP (number of x  and z variables) and
    ZINB (number of x  and z variables +1). If not specificied then the starting
    values are estimated from the data. It is advisable to use the default values.}
   \item{data}{a frame containing the variables in the model.}
   \item{sub}{an optional logical expression specifying a subset of observations to be
           used in the fitting process. This should be specified in the form: e.g. 
           subsets =expression(w==0).}
   \item{method}{the method to be used in fitting the model. The default
           method  is'"BFGS"', which  is a quasi-Newton method.
           For other method see optim() function in R.}
    \item{ntries}{The number of times the optimization may be repeated till convergence.
	The default value is 5. Normally this is performed once.} 
    \item{...}{other options to optim function.} 
}
\value{
 \item{coefficients}{The best set of parameters found. Subscript x and z respectively
 	indicates the estimates from the regression and the zero-inflated part.}
 \item{se}{Standard error of the best set of parameters found.}
 \item{maxlike}{The maximum likelihood value.}
  \item{corr}{The correlation matrix of the regression coefficients.}
  \item{counts}{A two-element integer vector giving the number of calls to
	the likelihood and the first derivative of the likelihood respectively.}
 \item{convergence}{ An integer code. '0' indicates successful convergence, otherwise no convergence.}
 \item{message}{A character string giving any additional information returned
	by the optimizer, or 'NULL'.}
 \item{hessian}{ A symmetric matrix giving an estimate of the Hessian at the solution found.}
 \item{data}{A vector of the response variable.} 
 \item{ntries}{The number of times the optimization is performed.} 
 \item{aic}{Akaike information criteria.}
 \item{bic}{Bayesian information criteria.} 
 }
  
\references{
Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U. Kirchner (1999). 
	The zero-inflated Poisson model and the decayed, missing and filled 
	teeth index in dental epidemiology. Journal of the Royal Statistical Society,
	Series A 162, 195--209.
	
Carrivick, P. J. W., A. H. Lee, and K. K. W. Yau (2003). Zero-inflated Poisson 
	to evaluate occupational safety interventions.
	Safety Science 41, 53--63.

Lambert, D. (1992). Zero-inflated Poisson regression, with an application to defects in
	manufacturing. Technometrics 34, 1--14.

Lewsey, J. D. and W. M. Thomson (2004). The utility of the zero-inflated Poisson and
	zero-inflated negative binomial models: a case study of cross-sectional and
	longitudinal DMF data examining the effect of socio-economic status. 
	Community Dentistry and Oral Epidemiology 32, 183--189.
	
Mwalili, S. M., Lesaffre, E. and Declerck, D. (2004). The interval censored zero-inflated
	negative binomial regression model: an application in caries research. submitted.

}

\seealso{
\code{\link{zicounts}, \link{ezinb}}.
}

\examples{
## load the artificial 'teeth' data
 data(teeth)
 names(teeth)
  
## a) fit a Poisson regression model
pois.cens1 <- zicensor(parm=c(-2,-2,1),resp = dmft~., upper = decid~.,x =~gender + age, 
	data=teeth, distr = "Poisson")
pois.cens1 
# fit for  boys only -- using 'sub'
pois.boys <- zicensor(resp = dmft~., upper = decid~.,x =~age, data=teeth, distr = "Poisson",
	sub=expression(gender==0))
pois.boys  

## b) fit zero-inflated  Poisson regression model
zip.cens <- zicensor(resp = dmft~., upper = decid~.,x =~gender + age,z =~gender + age, 
	data=teeth, distr = "ZIP")
zip.cens
#summary(zip.cens)

## c) fit negative binomial regression model
nb.cens <- zicensor(resp = dmft~., upper = decid~.,x =~gender + age,data=teeth, distr = "NB")
nb.cens
#summary(nb.cens)

## d) fit zero-inflated negative binomial regression model
zinb.cens <- zicensor(resp = dmft~., upper = decid~.,x =~gender + age,z =~gender + age, data=teeth)
zinb.cens
#summary(zinb.cens)
}

\keyword{models}
\author{S. M. Mwalili <samuelmwalili@yahoo.com>}
