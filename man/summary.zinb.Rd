\name{summary.zinb}
\alias{summary.zinb}
\alias{summary.czinb}
\alias{print.zinb}
\alias{print.czinb}
\alias{prnt.zinb.int}
\title{Printing summary of the (censored) zinb regression model}
\description{
    'print.zinb' and  'summary.zinb', or 'print.czinb' and  'summary.czinb' 
	provides an efficient way of summarizing the 'zicounts' or 'zicensor'
	object. NOTE 'prnt.zinb.int' is used ONLY for model without covariates.
}
\usage{
 print.zinb(x, digits = max(3, getOption("digits") - 3),...)

 summary.zinb(object, digits = max(3, getOption("digits") - 3),...)

 print.czinb(x, digits = max(3, getOption("digits") - 3),...)

 summary.czinb(object, digits = max(3, getOption("digits") - 3),...)

 prnt.zinb.int(object, digits = 1)

}


\arguments{
  \item{x}{an object from the 'zicounts' or 'zicensor' fit}
  \item{object}{an object from the 'zicounts' or 'zicensor' fit}
  \item{digits}{number of significant decimal points.}  
  \item{...}{other options 'summary' arguments.}  

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
\code{\link{zicounts}, \link{zicensor}}.
}

\examples{
## See help for ?zicounts or ?zicensor
}

\keyword{models}
\author{S. M. Mwalili <samuelmwalili@yahoo.com>}
