.packageName <- "zicounts"
#------------------------------------------------------------------------------------
#        Classical and Zero-inflated Models for Count Data  Regression model       
#------------------------------------------------------------------------------------
zicounts <- function(parm = NULL, resp = y ~ ., x =  ~ 1, z =  ~ 1, data, distrname = 
	"ZINB", offset=NULL,sub = NULL, ntries = 5, method = "BFGS", ...)
{
	# -log(likelihood)               
	if(distrname == "Poisson") {
		neg.like <- function(p, y,ee,xx, ...)
		{
			.C("PoisNLL",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZIP") {
		neg.like <- function(p, y,ee,xx, zz, ...)
		{
			.C("ZipNLL",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "NB") {
		neg.like <- function(p, y,ee,xx, ...)
		{
			.C("NegbinNLL",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZINB") {
		neg.like <- function(p, y,ee,xx, zz, ...)
		{
			.C("ZinbNLL",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	if(distrname == "Poisson") {
		neg.grad <- function(p, y,ee,xx, ...)
		{
			.C("PoisGrad",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(ncol(xx)),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZIP") {
		neg.grad <- function(p, y,ee,xx, zz, ...)
		{
			.C("ZipGrad",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(ncol(xx) + ncol(zz)),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "NB") {
		neg.grad <- function(p, y,ee,xx, ...)
		{
			.C("NegbinGrad",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(ncol(xx) + 1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZINB") {
		neg.grad <- function(p, y,ee,xx, zz, ...)
		{
			.C("ZinbGrad",
				as.double(p),
				as.double(y),
				as.double(ee),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(ncol(xx) + ncol(zz) + 1),
				PACKAGE="zicounts")$ans
		}
	}
	if(!is.null(sub))
		datset <- data[with(data, eval(sub)),,drop=FALSE]
	else datset <- data
	xx <- model.matrix(x, datset)
	zz <- model.matrix(z, datset)
	y <- as.matrix(model.response(model.frame(resp, datset)))


	if(!is.null(offset))
	  ee <- exp(model.matrix( ~ eval(offset) - 1, datset)) 
	else   ee <- matrix(1,nrow=dim(xx)[1])
         
        
	size <- nrow(xx)
	# if parm not given estimate the initial value
	if(is.null(parm)) {
		reg0 <- glm.fit(xx, y, family = poisson())
		zi0 <- glm.fit(zz, 1 - pmin(y, 1), family = binomial())
	}
	if(distrname == "Poisson") {
		if(is.null(parm))
			parm <- c(reg0$coef)
		pnames <- paste(dimnames(xx)[[2]], "x", sep = "")
		distrname.name <- "Poisson"
		llac <- paste("Formula = ", as.character(resp)[2], " ~ ", as.character(
			x)[2])
		if(length(parm) < (ncol(xx)))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx)))
			stop("Error in 'parm': more starting values specified")
	}
	else if(distrname == "ZIP") {
		if(is.null(parm))
			parm <- c(reg0$coef, zi0$coef)
		pnames <- c(paste(dimnames(xx)[[2]], "x", sep = ""), paste(dimnames(zz)[[
			2]], "z", sep = ""))
		distrname.name <- "Zero-Inflated Poisson"
		llac <- paste("Formula = ", as.character(resp)[2], " ~ ", as.character(
			x)[2], "|", as.character(z)[2])
		if(length(parm) < (ncol(xx) + ncol(zz)))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx) + ncol(zz)))
			stop("Error in 'parm': more starting values specified")
	}
	else if(distrname == "NB") {
		if(is.null(parm))
			parm <- c(reg0$coef, var(y)/mean(y))
		pnames <- c(paste(dimnames(xx)[[2]], "x", sep = ""), "log(tau)")
		distrname.name <- "Negative Binomial"
		llac <- paste("Formula = ", as.character(resp)[2], " ~ ", as.character(
			x)[2])
		if(length(parm) < (ncol(xx) + 1))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx) + 1))
			stop("Error in 'parm': more starting values specified")
	}
	else if(distrname == "ZINB") {
		if(is.null(parm))
			parm <- c(reg0$coef, zi0$coef, var(y)/mean(y))
		pnames <- c(paste(dimnames(xx)[[2]], "x", sep = ""), paste(dimnames(zz)[[
			2]], "z", sep = ""), "log(tau)")
		distrname.name <- "Zero-Inflated Negative Binomial"
		llac <- paste("Formula = ", as.character(resp)[2], " ~ ", as.character(
			x)[2], "|", as.character(z)[2])
		if(length(parm) < (ncol(xx) + ncol(zz) + 1))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx) + ncol(zz) + 1))
			stop("Error in 'parm': more starting values specified")
	}
	z0 <- optim(par = parm, fn = neg.like, gr = neg.grad, hessian = TRUE, method = method,
		control = list(maxit = 250), y = y, xx = xx, zz = zz, ee=ee)
	t <- 1
	while(z0$conv != 0) {
		z0 <- optim(z0$par, fn = neg.like, gr = neg.grad, hessian = TRUE, method = 
			method, control = list(maxit = 250), y = y, xx = xx, zz = zz,ee=ee)
		t <- t + 1
		if(t == ntries)
			break
	}
	if(!is.null(offset))
	   llac <- paste(llac,"(with offest)")
	else   llac <- llac
	
	cov <- solve(z0$hessian)
	se <- sqrt(diag(cov))
	names(se) <- pnames
	estimate <- array(z0$par)
	names(estimate) <- pnames
	corr <- cov/(se %o% se)
	colnames(corr) <- pnames
	rownames(corr) <- pnames
	z1 <- list(distrname = distrname.name, call = llac, coefficients = estimate, se = 
		se, maxlike = z0$value, corr = corr, aic = 2 * z0$value + 2 * length(
		estimate), bic = 2 * z0$value + log(size) * length(estimate), df = size -
		length(estimate), counts = z0$counts, convergence = z0$convergence, ntries
		 = t, message = z0$message, hessian = z0$hessian, method = method, data = 
		y)
	class(z1) <- c("zinb")
	return(z1)
}

#---------------------------------------------------------------------------------------
#           C-e-n-s-o-r-e-d       Count Data  Regression model             
#---------------------------------------------------------------------------------------
zicensor <- function(parm = NULL, resp = y ~ ., upper = r ~ ., x =  ~ 1, z =  ~ 1, data,
	distrname = "ZINB", sub = NULL, ntries = 5, method = "BFGS", ...)
{
	# -log(likelihood)               
	if(distrname == "Poisson") {
		neg.like <- function(p, y, r, xx, ...)
		{
			.C("CensPois",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZIP") {
		neg.like <- function(p, y, r, xx, zz, ...)
		{
			.C("CensZip",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "NB") {
		neg.like <- function(p, y, r, xx, ...)
		{
			.C("CensNegbin",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZINB") {
		neg.like <- function(p, y, r, xx, zz, ...)
		{
			.C("CensZinb",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(1),
				PACKAGE="zicounts")$ans
		}
	}
	# gradient-log(likelihood)               
	if(distrname == "Poisson") {
		grad.like <- function(p, y, r, xx, ...)
		{
			.C("CensPoisGrad",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(ncol(xx)),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZIP") {
		grad.like <- function(p, y, r, xx, zz, ...)
		{
			.C("CensZipGrad",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(ncol(xx) + ncol(zz)),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "NB") {
		grad.like <- function(p, y, r, xx, ...)
		{
			.C("CensNegbinGrad",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				ans = double(ncol(xx) + 1),
				PACKAGE="zicounts")$ans
		}
	}
	else if(distrname == "ZINB") {
		grad.like <- function(p, y, r, xx, zz, ...)
		{
			.C("CensZinbGrad",
				as.double(p),
				as.double(y),
				as.double(r),
				as.double(xx),
				as.integer(nrow(xx)),
				as.integer(ncol(xx)),
				as.double(zz),
				as.integer(nrow(zz)),
				as.integer(ncol(zz)),
				ans = double(ncol(xx) + ncol(zz) + 1),
				PACKAGE="zicounts")$ans
		}
	}
	# create a subset of the data
	if(!is.null(sub))
		datset <- data[with(data, eval(sub)),,drop=FALSE]
	else datset <- data
	xx <- model.matrix(x, datset)
	zz <- model.matrix(z, datset)
	y <- as.matrix(model.response(model.frame(resp, datset)))
	r <- as.matrix(model.response(model.frame(upper, datset)))
	size <- nrow(xx)
	# if parm not given estimate the initial value
	if(is.null(parm)) {
		reg0 <- glm.fit(xx, y, family = poisson())
		zi0 <- glm.fit(zz, 1 - pmin(y, 1), family = binomial())
	}
	if(distrname == "Poisson") {
		if(is.null(parm))
			parm <- c(reg0$coef)
		pnames <- paste(dimnames(xx)[[2]], "x", sep = "")
		distrname.name <- "Censored Poisson"
		llac <- paste("Formula = (", as.character(resp)[2], ",", as.character(
			upper)[2], ") ~ ", as.character(x)[2])
		if(length(parm) < (ncol(xx)))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx)))
			stop("Error in 'parm': more starting values specified")
	}
	else if(distrname == "ZIP") {
		if(is.null(parm))
			parm <- c(reg0$coef, zi0$coef)
		pnames <- c(paste(dimnames(xx)[[2]], "x", sep = ""), paste(dimnames(
			zz)[[2]], "z", sep = ""))
		distrname.name <- "Censored Zero-Inflated Poisson"
		llac <- paste("Formula = (", as.character(resp)[2], ",", as.character(
			upper)[2], ") ~ ", as.character(x)[2], "|", as.character(z)[
			2])
		if(length(parm) < (ncol(xx) + ncol(zz)))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx) + ncol(zz)))
			stop("Error in 'parm': more starting values specified")
	}
	else if(distrname == "NB") {
		if(is.null(parm))
			parm <- c(reg0$coef, var(y)/mean(y))
		pnames <- c(paste(dimnames(xx)[[2]], "x", sep = ""), "log(tau)")
		distrname.name <- "Censored Negative Binomial"
		llac <- paste("Formula = (", as.character(resp)[2], ",", as.character(
			upper)[2], ") ~ ", as.character(x)[2])
		if(length(parm) < (ncol(xx) + 1))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx) + 1))
			stop("Error in 'parm': more starting values specified")
	}
	else if(distrname == "ZINB") {
		if(is.null(parm))
			parm <- c(reg0$coef, zi0$coef, var(y)/mean(y))
		pnames <- c(paste(dimnames(xx)[[2]], "x", sep = ""), paste(dimnames(
			zz)[[2]], "z", sep = ""), "log(tau)")
		distrname.name <- "Censored Zero-Inflated Negative Binomial"
		llac <- paste("Formula = (", as.character(resp)[2], ",", as.character(
			upper)[2], ") ~ ", as.character(x)[2], "|", as.character(z)[
			2])
		if(length(parm) < (ncol(xx) + ncol(zz) + 1))
			stop("Error in 'parm': less starting values specified")
		if(length(parm) > (ncol(xx) + ncol(zz) + 1))
			stop("Error in 'parm': more starting values specified")
	}
	z0 <- optim(par = parm, fn = neg.like, gr = grad.like, hessian = TRUE, method = 
		method, control = list(maxit = 250), y = y, r = r, xx = xx, zz = zz)
	t <- 1
	while(z0$conv != 0) {
		z0 <- optim(z0$par, fn = neg.like, gr = grad.like, hessian = TRUE, method
			 = method, control = list(maxit = 250), y = y, r = r, xx = xx,
			zz = zz)
		t <- t + 1
		if(t == ntries)
			break
	}
	cov <- solve(z0$hessian)
	se <- sqrt(diag(cov))
	names(se) <- pnames
	estimate <- array(z0$par)
	names(estimate) <- pnames
	corr <- cov/(se %o% se)
	colnames(corr) <- pnames
	rownames(corr) <- pnames
	z1 <- list(distrname = distrname.name, call = llac, coefficients = estimate,
		se = se, maxlike = z0$value, corr = corr, aic = 2 * z0$value + 2 * 
		length(estimate), bic = 2 * z0$value + log(size) * length(estimate),
		df = size - length(estimate), counts = z0$counts, convergence = z0$
		convergence, ntries = t, message = z0$message, hessian = z0$hessian,
		method = method, data = y)
	class(z1) <- c("czinb")
	return(z1)
}

# expected value of zinb
#-----------------------
ezinb <- function(coeff = NULL, resp = y ~ ., upper = r ~ ., x =  ~ 1, z =  ~ 1, data=NULL,
	sub = NULL)
{
	# create a subset of the data
	if(!is.null(sub))
		datset <- data[with(data, eval(sub)),,drop=FALSE]
	else datset <- data
	xx <- model.matrix(x, datset)
	zz <- model.matrix(z, datset)
	y <- as.matrix(model.response(model.frame(resp, datset)))
	r <- as.matrix(model.response(model.frame(upper, datset)))
	if(length(coeff) < (ncol(xx) + ncol(zz) + 1))
		stop("Error in 'parm': less starting values specified")
	if(length(coeff) > (ncol(xx) + ncol(zz) + 1))
		stop("Error in 'parm': more starting values specified")
	out <- .C("EZinb",
		as.double(coeff),
		as.double(y),
		as.double(r),
		as.double(xx),
		as.integer(nrow(xx)),
		as.integer(ncol(xx)),
		as.double(zz),
		as.integer(nrow(zz)),
		as.integer(ncol(zz)),
		ans = double(nrow(xx)),
		PACKAGE="zicounts")$ans
	return(out)
}

#------------------------------------------------------------------
#           utility functions:  print and summary of "zinb" class
#------------------------------------------------------------------
print.zinb <- function(x, digits = max(3, getOption("digits") - 3),...)
{
	cat("\n")
	cat(x$distrname, "\n")
	cat("\nCall: \n", x$call, "\n\n")
	if(x$distrname == "Zero-Inflated Poisson") {
		cat("Subscripts:\n")
		cat("\t-- x: coefficients of Poisson part\n")
		cat("\t-- z: coefficients of zero-inflated part\n\n")
	}
	if(x$distrname == "Zero-Inflated Negative Binomial") {
		cat("\nSubscripts:\n")
		cat("\t-- x: coefficients of of negative binomial part\n")
		cat("\t-- z: coefficients of zero-inflated part\n\n")
	}
	if(length(coef(x))) {
		cat("Coefficients")
		cat(":\n")
		print.default(format(x$coefficients, digits = digits), print.gap = 2,
			quote = FALSE)
	}
	else cat("No coefficients\n\n")
	cat("\nDegrees of Freedom:", x$df, "\n")
	cat("Deviance:    ", format(signif(2 * x$max, digits)), "\t AIC:", format(signif(
		x$aic, digits)), "\t BIC:", format(signif(x$bic, digits)), "\n")
	cat("\nNumber of iterations", x$method, "is called before convergence:", x$
		ntries, "\n")
	invisible(x)
}

 summary.zinb <- function(object, digits = max(3, getOption("digits") - 3),...)
{
	tvalue <- object$coeff/object$se
	dn <- c("Estimate", "Std. Error")
	pvalue <- 2 * pnorm( - abs(tvalue))
	cl <- t(apply(cbind(object$coeff, object$se), 1, function(x)
	c(x[1] - 1.96 * x[2], x[1] + 1.96 * x[2])))
	coef.table <- cbind(object$coeff, object$se, cl, tvalue, pvalue)
	dimnames(coef.table) <- list(names(object$coeff), c(dn, "lower", "upper", "z value",
		"Pr(>|z|)"))
	ans <- list()
	ans$coefficients <- coef.table
	ans$distrname <- object$distrname
	ans$call <- object$call
	ans$aic <- object$aic
	ans$bic <- object$bic
	ans$iter <- object$iter
	ans$max <- object$max
	ans$df <- object$df
	ans$ntries <- object$ntries
	ans$cov <- object$cov
	ans$method <- object$method
	ans$corr <- object$corr
	class(ans) <- c("zinb")
	ans
}

 print.czinb <- function(x, digits = max(3, getOption("digits") - 3),...)
{
	cat("\n")
	cat(x$distrname, "\n")
	cat("\nCall: \n", x$call, "\n\n")
	if(x$distrname == "Censored Zero-Inflated Poisson") {
		cat("Subscripts:\n")
		cat("\t-- x: coefficients of censored Poisson part \n")
		cat("\t-- z: coefficients of censored zero-inflated part\n\n")
	}
	if(x$distrname == "Censored Zero-Inflated Negative Binomial") {
		cat("\nSubscripts:\n")
		cat("\t-- x: coefficients of censored negative binomial part\n")
		cat("\t-- z: coefficients of censored zero-inflated part\n\n")
	}
	if(length(coef(x))) {
		cat("Coefficients")
		cat(":\n")
		print.default(format(x$coefficients, digits = digits), print.gap = 2,
			quote = FALSE)
	}
	else cat("No coefficients\n\n")
	cat("\nDegrees of Freedom:", x$df, "\n")
	cat("Deviance:    ", format(signif(2 * x$max, digits)), "\t AIC:", format(signif(
		x$aic, digits)), "\t BIC:", format(signif(x$bic, digits)), "\n")
	cat("\nNumber of iterations", x$method, "is called before convergence:", x$
		ntries, "\n")
	invisible(x)
}

 summary.czinb <- function(object, digits = max(3, getOption("digits") - 3),...)
{
	tvalue <- object$coeff/object$se
	dn <- c("Estimate", "Std. Error")
	pvalue <- 2 * pnorm( - abs(tvalue))
	cl <- t(apply(cbind(object$coeff, object$se), 1, function(x)
	c(x[1] - 1.96 * x[2], x[1] + 1.96 * x[2])))
	coef.table <- cbind(object$coeff, object$se, cl, tvalue, pvalue)
	dimnames(coef.table) <- list(names(object$coeff), c(dn, "lower", "upper", "z value",
		"Pr(>|z|)"))
	ans <- list()
	ans$coefficients <- coef.table
	ans$distrname <- object$distrname
	ans$call <- object$call
	ans$aic <- object$aic
	ans$bic <- object$bic
	ans$iter <- object$iter
	ans$max <- object$max
	ans$df <- object$df
	ans$ntries <- object$ntries
	ans$cov <- object$cov
	ans$method <- object$method
	ans$corr <- object$corr
	class(ans) <- c("czinb")
	ans
}


#----------------------------------------------------------------
# prints the fitted values for the model without covariates      
#----------------------------------------------------------------
prnt.zinb.int <- function(object, digits = 1)
{
	# making table and a summary for a model without covariates
	y <- object$data
	xx <- table(y)
	xy <- as.numeric(names(xx))
	nx <- length(xx)
	fit <- 1:nx
	cat("\n\t", object$distrname, "\n")
	p <- 1/(1 + exp( - object$coeff[2]))
	L <- exp(object$coeff[1])
	t <- exp(object$coeff[3])
	mu <- L * (1 - p)
	v <- L * (1 - p) * (1 + L * p + L/t)
	fit[1] <- (p + (1 - p) * exp( - t * log(1 + L/t))) * sum(xx)
	for(i in 2:nx)
		fit[i] <- (1 - p) * exp(lgamma(xy[i] + t) - lgamma(xy[i] + 1) - lgamma(t) - t *
			log(1 + L/t) - xy[i] * log(1 + t/L)) * sum(xx)
	return(list(p = as.vector(p), lambda = as.vector(L), tau = as.vector(t), mean = 
		as.vector(mu), variance = as.vector(v), caries.free = fit[1]/length(y), fit = 
		data.frame(y = xy, n.obs = array(xx), n.fit = round(fit, digits)), aic = object$aic))
}


 .First.lib <- function(lib, pkg)
{
	library.dynam("zicounts", pkg, lib)
}