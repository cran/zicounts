#include <math.h>
double digamma(double*);

extern "C"{ 
void PoisNLL(double* p, double* y, double* ee,double* x,int* nrowx, int* ncolx, double* nll)
{
// negative log-likelihood of poisson model
	int i, j,  ncol = *ncolx, n=*nrowx;
	double loglambda,neglik,lambda;
	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*n]*p[j];
	 	 lambda = exp(loglambda)*ee[i];
	 	 loglambda = log(lambda);

		neglik += lambda - y[i]*loglambda + lgamma(y[i] +1);
	}
	*nll = neglik;
}


void ZipNLL(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* nll)
{
// negative log-likelihood of zero-inflated Poisson
	int i, j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double loglambda,neglik,logitp,lambda,pp;

	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
	 pp = 1/(1+exp(-logitp));
	 if(y[i]==0)
	   neglik += -log(pp + (1-pp)*exp(-lambda));
	    else
	     neglik += -log(1-pp) + exp(loglambda) - y[i]*loglambda + lgamma(y[i] +1);
	}
	*nll = neglik;
}

void NegbinNLL(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx, double* nll)
{
// negative log-likelihood of negative binomial
	int i, j,  ncol = *ncolx, n=*nrowx;
	double loglambda,lambda,neglik,logtau, tau;
	logtau =p[ncol];
	tau = exp(logtau);
	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*n]*p[j];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
		neglik += lgamma(y[i]+1) + lgamma(tau) - lgamma(tau + y[i]) -
		tau * log(tau/(tau + lambda)) - y[i]* log(lambda/(tau + lambda));
	}
	*nll = neglik;
}


void ZinbNLL(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* nll)
{
// negative log-likelihood of zero-inflated negative binomial
	int i, j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double loglambda,neglik,logitp,lambda,pp,logtau, tau;
	logtau =p[colx+colz];
	tau = exp(logtau);

	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
	 pp = 1/(1+exp(-logitp));
	 if(y[i]==0)
	   neglik += -log(pp + (1-pp)*exp(tau * log(tau/(tau + lambda))));
	    else
	     neglik += -log(1-pp) + lgamma(y[i]+1) + lgamma(tau) - lgamma(tau + y[i]) -
		      tau * log(tau/(tau + lambda)) - y[i]* log(lambda/(tau + lambda));
	}
	*nll = neglik;
}

void CensPois(double* p, double* y, double* r, double* x, int* nrowx, int* ncolx, double* nll)
{
// negative log-likelihood  of censored poisson model
	int i, j,  ncol = *ncolx, n=*nrowx;
	double lik,yt,loglambda,neglik, lambda;
	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	 lik = 0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*n]*p[j];
	 lambda = exp(loglambda);
	for(yt=y[i]; yt<=r[i]; yt++){
		   lik += exp(-lambda + yt*loglambda -lgamma(yt +1));
	   }
     neglik += -log(lik);
	}
	*nll = neglik;
}

void CensZinb(double* p, double* y, double* r,double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* nll)
{
// negative log-likelihood of censored zero-inflated negative binomial
	int i,j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double yt,loglambda,neglik,lik,logitp,lambda,pp,logtau, tau;
	logtau =p[colx+colz];
	tau = exp(logtau);

	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;
	lik = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda);
	 pp = 1/(1+exp(-logitp));
	 for(yt=y[i]; yt<=r[i]; yt++){
	 if(yt==0)
	   lik += pp + (1-pp)*exp(tau * log(tau/(tau + lambda)));
	    else
	     lik += (1-pp)*exp(lgamma(tau + yt) + tau * log(tau/(tau + lambda)) +
                           yt*log(lambda/(tau + lambda))  - lgamma(yt+1) - lgamma(tau));
	//     lik += exp(log(1-pp) - lgamma(yt+1) - lgamma(tau) + lgamma(tau + yt) +
	//	      tau * log(tau/(tau + lambda)) + yt* log(lambda/(tau + lambda)));
		      }
		neglik += -log(lik);
	}
	*nll = neglik;
}



void CensNegbin(double* p, double* y, double* r, double* x, int* nrowx, int* ncolx, double* nll)
{
// negative log-likelihood of censored negative binomial
	int i, j,  ncol = *ncolx, n=*nrowx;
	double yt, lik, loglambda,lambda,neglik,logtau, tau;
	logtau =p[ncol];
	tau = exp(logtau);
	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	lik = 0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*n]*p[j];
	        lambda = exp(loglambda);
	  for(yt=y[i]; yt<=r[i]; yt++){
	        lik += exp(-lgamma(yt+1) - lgamma(tau) + lgamma(tau + yt) +
	 	tau * log(tau/(tau + lambda)) + yt* log(lambda/(tau + lambda)));
	   }
	   neglik += -log(lik);
	}
	*nll = neglik;
}


void CensZip(double* p, double* y, double* r, double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* nll)
{
// negative log-likelihood of zero-inflated Poisson
	int i, j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double yt, lik,loglambda,neglik,logitp,lambda,pp;

	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;
	lik = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda);
	 pp = 1/(1+exp(-logitp));
	 for(yt=y[i]; yt<=r[i]; yt++)
	{
	 if(yt==0)
	    lik  +=  pp + (1-pp)*exp(-lambda);
	    else
	     lik += exp(log(1-pp) - lambda + yt*loglambda - lgamma(yt +1));
	  }
	  neglik += -log(lik);
	}
	*nll = neglik;
  }

void EZinb(double* p, double* y, double* r,double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* xpcti)
{
// expected value of the  zero-inflated negative binomial
	int i,j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double xpct,yt,loglambda,prob,sumprob,logitp,lambda,pp,logtau, tau;
	logtau =p[colx+colz];
	tau = exp(logtau);

	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda);
	 pp = 1/(1+exp(-logitp));
  
	xpct = sumprob = 0;
	 for(yt=y[i]; yt<=r[i]; yt++){
	 if(yt==0) {
	    prob = pp + (1-pp)*exp(tau * log(tau/(tau + lambda)));
	    sumprob += prob;
	    xpct += 0;
	   }
	    else {
	      prob = (1-pp)*exp(lgamma(tau + yt) + tau * log(tau/(tau + lambda)) +
                           yt*log(lambda/(tau + lambda))  - lgamma(yt+1) - lgamma(tau));
              sumprob += prob;
              xpct += prob*yt;
              }
	}
	xpcti[i] = xpct/sumprob;
     }
  }

void PoisGrad(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx,double* grad)
{
// first derivative: negative log-likelihood of censored zero-inflated negative binomial
	int i,j, n=*nrowx, colx = *ncolx;
	double loglambda,lambda;
	for(j=0; j<colx;j++)  grad[j] = 0; // initialize grad

	for(i=0; i<n; i++)
	{
	loglambda =0;
	loglambda =0;
	
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
	    	 for(j=0; j<colx;j++)
	    	 grad[j]         += x[i + j*n]*(lambda-y[i]);
         }
   
}
 
void ZipGrad(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* grad)
{
// first derivative: negative log-likelihood of censored zero-inflated negative binomial
	int i,j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double loglambda,logitp,lambda,pp,eG;
	int ncoefpl = *ncolx + *ncolz;
	for(j=0; j<ncoefpl;j++)  grad[j] = 0; // initialize grad

	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;

	loglambda =0;
	logitp = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);

	 eG = exp(logitp);
         pp = 1/(1+ exp(-logitp));


	 if(y[i]==0){
	 	 for(j=0; j<colx;j++)
	 	        grad[j] += x[i + j*n]*lambda/(1 + exp(lambda)*eG) ;
	 	 for(k=0; k<colz;k++)
	 	     grad[colx + k] += z[i + k*n]*(-pp*(exp(lambda)-1)/(1 + exp(lambda)*eG));

	 }   else {
	    	 for(j=0; j<colx;j++)
	    	 grad[j]         += x[i + j*n]*(lambda-y[i]);
	    	 for(k=0; k<colz;k++)
	    	 grad[colx + k]  += z[i + k*n]*pp; 

         }
   
   }
}
 

void NegbinGrad(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx, double* grad)
{
// first derivative: negative log-likelihood of censored zero-inflated negative binomial
	int i,j, n=*nrowx, colx = *ncolx;
	double loglambda,lambda,logtau, tau_yt,tau;
	for(j=0; j<colx;j++)  grad[j] = 0; // initialize grad
	logtau =p[colx];
	tau = exp(logtau);

	for(i=0; i<n; i++)
	{
	loglambda =0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
	    	 for(j=0; j<colx;j++)
	    	 //grad[j]         += x[i + j*n]*(1-(tau + y[i])/(lambda + tau))*tau;
                   grad[j]         += x[i + j*n]*tau*(lambda-y[i])/(lambda + tau);
          tau_yt = tau + y[i];
          grad[colx] +=   -1  + (tau + y[i])/(lambda + tau) + log((lambda + tau)/tau)
                               + digamma(&tau) - digamma(&tau_yt);
   }
  grad[colx]  = tau*grad[colx];
}




void ZinbGrad(double* p, double* y, double* ee, double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* grad)
{
// first derivative: negative log-likelihood of censored zero-inflated negative binomial
	int i,j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
	double loglambda,logitp,lambda,pp,eG,logtau, tau_yt,tau;
	int ncoefpl = *ncolx + *ncolz + 1;
	for(j=0; j<ncoefpl;j++)  grad[j] = 0; // initialize grad
	logtau =p[colx+colz];
	tau = exp(logtau);

	for(i=0; i<n; i++)
	{
	loglambda =0;
	logitp = 0;

	loglambda =0;
	logitp = 0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
	 eG = exp(logitp);
         pp = 1/(1+ exp(-logitp));


	 if(y[i]==0){
	 	 for(j=0; j<colx;j++)
	 	        grad[j] += x[i + j*n]*(lambda*tau)/((lambda + tau)*(1 +eG*(exp(tau*log((lambda + tau)/tau))))) ;
	 	 for(k=0; k<colz;k++)
	 	     grad[colx + k] += z[i + k*n]*(-1/(1+eG) + 1/(1 +eG*exp(tau*log((lambda + tau)/tau))));
	    	 grad[colx+colz] += (-lambda + (lambda + tau)* log ((lambda + tau)/tau)) /
                              ((lambda + tau)*(1 + eG*(exp(tau*log((lambda + tau)/tau)))));

	 }   else {
	    	 for(j=0; j<colx;j++)
	    	 grad[j]         += x[i + j*n]*(1-(tau + y[i])/(lambda + tau))*tau;
	    	 for(k=0; k<colz;k++)
	    	 grad[colx + k]  += z[i + k*n]*pp; 

          tau_yt = tau + y[i];
          grad[colx+colz] +=   -1  + (tau + y[i])/(lambda + tau) + log((lambda + tau)/tau)
                               + digamma(&tau) - digamma(&tau_yt);
   }
   
   }
 grad[colx+colz]  = tau*grad[colx+colz];
}

void CensPoisGrad(double* p, double* y, double* r,double* x,
        int* nrowx, int* ncolx,double* grad)
{
// first derivative: negative log-likelihood of censored negative binomial
    int i,j, n=*nrowx, colx = *ncolx;
    double yt,loglambda,lambda,lik,sumlik;
    
    double grad2[colx];
    for(j=0; j<colx;j++)  grad[j] = 0; // initialize grad
         
    for(i=0; i<n; i++)
    {
    loglambda =0;

     for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
     lambda = exp(loglambda);
 
      sumlik = 0;
     for(j=0; j<colx;j++)  grad2[j] = 0; // initialize grad2
     for(yt=y[i]; yt<=r[i]; yt++){

             lik = exp(-lambda + yt*loglambda -lgamma(yt +1));
               
             for(j=0; j<colx;j++)
                  grad2[j]   += lik*x[i + j*n]*(lambda-yt);
 
     sumlik += lik;  
        } // yt
      for(j=0; j<colx;j++)  grad[j] += grad2[j]/sumlik;
     
    } //i
    
  }  
  
void CensZipGrad(double* p, double* y, double* r,double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* grad)
{
// first derivative: negative log-likelihood of censored zip
    int i,j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
    double yt,loglambda,logitp,lambda,pp,eG,lik,sumlik;
    
    int ncoefpl = *ncolx + *ncolz;
    double grad2[ncoefpl];
    for(j=0; j<ncoefpl;j++)  grad[j] = 0; // initialize grad
         
    for(i=0; i<n; i++)
    {
    loglambda =0;
    logitp = 0;

     for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
     for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
     lambda = exp(loglambda);
     eG = exp(logitp);
         pp = 1/(1+ exp(-logitp));

      sumlik = 0;
     for(j=0; j<ncoefpl;j++)  grad2[j] = 0; // initialize grad2
     for(yt=y[i]; yt<=r[i]; yt++){

     if(yt==0){
         lik =  pp + (1-pp)*exp(-lambda);;
         //lik =  pp + (1-pp)*pow((lambda + tau)/tau,-tau);
         for(j=0; j<colx;j++)
                    grad2[j] += lik*x[i + j*n]*lambda/(1 + exp(lambda)*eG);
         for(k=0; k<colz;k++)
             grad2[colx + k] += lik*z[i + k*n]*(-pp*(exp(lambda)-1)/(1 + exp(lambda)*eG));
                                    

        }
        else {
             lik = (1-pp)*exp(- lambda + yt*loglambda - lgamma(yt +1));
               
             for(j=0; j<colx;j++)
             grad2[j]         += lik*x[i + j*n]*(lambda-y[i]);
             for(k=0; k<colz;k++)
             grad2[colx + k]  += lik*z[i + k*n]*pp;

        }

       sumlik += lik;  
        } // yt
      for(j=0; j<ncoefpl;j++)  grad[j] += grad2[j]/sumlik;
     
    } //i
    
  } // fun
 
void CensNegbinGrad(double* p, double* y, double* r,double* x,
        int* nrowx, int* ncolx,double* grad)
{
// first derivative: negative log-likelihood of censored negative binomial
    int i,j, n=*nrowx, colx = *ncolx;
    double yt,loglambda,lambda,logtau, tau_yt,tau, lik,sumlik;
    
    int ncoefpl = *ncolx + 1;
    double grad2[ncoefpl];
    for(j=0; j<ncoefpl;j++)  grad[j] = 0; // initialize grad
    logtau =p[colx];
    tau = exp(logtau);
         
    for(i=0; i<n; i++)
    {
    loglambda =0;

     for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
     lambda = exp(loglambda);
 
      sumlik = 0;
     for(j=0; j<ncoefpl;j++)  grad2[j] = 0; // initialize grad2
     for(yt=y[i]; yt<=r[i]; yt++){

             lik = exp(lgamma(tau + yt) + tau * log(tau/(tau + lambda)) +
                                   yt*log(lambda/(tau + lambda))  - lgamma(yt+1) - lgamma(tau));
               
             for(j=0; j<colx;j++)
                  grad2[j]   += lik*x[i + j*n]*tau*(lambda-yt)/(lambda + tau);
 
          tau_yt = tau + yt;
          grad2[colx] +=  lik*(-1  + (tau + yt)/(lambda + tau) + log((lambda + tau)/tau)
                               + digamma(&tau) - digamma(&tau_yt));
       sumlik += lik;  
        } // yt
      for(j=0; j<ncoefpl;j++)  grad[j] += grad2[j]/sumlik;
     
    } //i
    
  grad[colx]  = tau*grad[colx];
  }  
  
  
void CensZinbGrad(double* p, double* y, double* r,double* x, int* nrowx, int* ncolx,
            double* z, int* nrowz, int* ncolz,double* grad)
{
// first derivative: negative log-likelihood of censored zero-inflated negative binomial
    int i,j,k, n=*nrowx, colx = *ncolx, colz=*ncolz;
    double yt,loglambda,logitp,lambda,pp,eG,logtau, tau_yt,tau, lik,sumlik;
    
    int ncoefpl = *ncolx + *ncolz + 1;
    double grad2[ncoefpl];
    for(j=0; j<ncoefpl;j++)  grad[j] = 0; // initialize grad
    logtau =p[colx+colz];
    tau = exp(logtau);
         
    for(i=0; i<n; i++)
    {
    loglambda =0;
    logitp = 0;

    loglambda =0;
    logitp = 0;
     for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
     for(k=0; k<colz;k++) logitp    += z[i + k*n]*p[colx + k];
     lambda = exp(loglambda);
     eG = exp(logitp);
         pp = 1/(1+ exp(-logitp));

      sumlik = 0;
     for(j=0; j<ncoefpl;j++)  grad2[j] = 0; // initialize grad2
     for(yt=y[i]; yt<=r[i]; yt++){

     if(yt==0){
         lik =  pp + (1-pp)*exp(tau * log(tau/(tau + lambda)));
         //lik =  pp + (1-pp)*pow((lambda + tau)/tau,-tau);
         for(j=0; j<colx;j++)
                    grad2[j] += lik*x[i + j*n]*(lambda*tau)/((lambda + tau)*(1 +eG*(exp(tau*log((lambda + tau)/tau))))) ;
         for(k=0; k<colz;k++)
             grad2[colx + k] += lik*z[i + k*n]*(-1/(1+eG) + 1/(1 +eG*exp(tau*log((lambda + tau)/tau))));
             grad2[colx+colz] += lik*(-lambda + (lambda + tau)* log ((lambda + tau)/tau)) /
                              ((lambda + tau)*(1 + eG*(exp(tau*log((lambda + tau)/tau)))));
                             

        }

        else {
             lik = (1-pp)*exp(lgamma(tau + yt) + tau * log(tau/(tau + lambda)) +
                                   yt*log(lambda/(tau + lambda))  - lgamma(yt+1) - lgamma(tau));
               
             for(j=0; j<colx;j++)
             grad2[j]         += lik*x[i + j*n]*(1-(tau + yt)/(lambda + tau))*tau;
             for(k=0; k<colz;k++)
             grad2[colx + k]  += lik*z[i + k*n]*pp;

          tau_yt = tau + yt;
          grad2[colx+colz] +=  lik*( -1  + (tau + yt)/(lambda + tau) + log((lambda + tau)/tau)
                              + digamma(&tau) - digamma(&tau_yt));
        }

       sumlik += lik;  
        } // yt
      for(j=0; j<ncoefpl;j++)  grad[j] += grad2[j]/sumlik;
     
    } //i
    
  grad[colx+colz]  = tau*grad[colx+colz];
  } // fun
 
} // extrn C

double digamma(double *x)
{
  const double h=1.0e-6;
  return ((lgamma(*x + h) - lgamma(*x - h))/(2 * h));
}

