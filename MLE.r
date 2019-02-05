#MLE Function

mle <- function (data, distr)
{
if (!is.character(distr)) 
	stop("distr must be a character string naming a distribution")
	else 
			distname <- distr
			
		#n is the length of data
		#m holds the mean of data
		#variance contains the variance of the data			
	
		{
			n <- length(data)
			m <- mean(data)
			v <- (n - 1)/n*var(data)
			sd=sqrt(v)
		}

	if (distname == "norm") {

			normalF <- function(parvec) {
			  # Log of likelihood of a normal distribution
			  
			  sum ( -((0.5* log(sd)) - 0.5*(((data - m)^2)/(sd))) )
			}
			MLE = 	optim(c(1,1), # initial values for mu and sigma
					fn = normalF, # function to maximize
					method = "L-BFGS-B", # this method lets set lower bounds
					lower = 0.00001, # lower limit for parameters
					control = list(fnscale = -1), # maximize the function
					hessian = T # calculate Hessian matricce because we will need for confidence intervals
					)
			
			return(MLE$par)
		}
		
	if (distname == "pois") {
			lambda=m
			
			# Log of likelihood of a poisson distribution
			poisson.lk <- function(lambda, y){
			log.lk <- sum(dpois(y, lambda, log=TRUE))
			return(-log.lk)
		}
			
			#MLE
			#MLE <- optim(c=m, poisson.lk, y=data, method="BFGS")
			return(lambda)
		}
	if (distname == "beta" ) {
			aux<-m*(1-m)/v - 1
			shape1 <- m*aux
			shape2 <- (1-m)*aux
			pars=c(shape1,shape2)
			# Log of likelihood of a beta distribution
			beta.lk <- function(pars, y){
			a <- pars[1]; b <- pars[2]
			log.lk <- sum(dbeta(y, a, b, log=TRUE))
			return(-log.lk)
		}
			MLE <- optim(c(1,1), beta.lk, y=data)
			
			return(MLE$par)
	   }
	   
	   
	if (distname == "bern" ) {
			i=0
			for (z in 1:n){
				if (data[z]==1)
					{i=i+1}
			}
			prob= i/n  
			pars=c(prob)
			# Log of likelihood of a bernoulli distribution
			LL <- function(pars, y){
			a <- pars[1]
			log.lk <- sum(dbern(y, a, log=TRUE))
			return(-log.lk)
		}
			
			MLE <- optim(par=prob, LL, y=y, method="BFGS")
			
			return(MLE$par)
	   
	   }
	   
	if (distname == "unif" ) {
			min1 <- m-(3*sd)
			max1 <- m+(3*sd)
			pars=c(min1,max1)
			# Log of likelihood of a uniform distribution
			LL <- function(pars, y){
			a <- pars[1]; b <- pars[2]
			log.lk <- sum(dunif(y, a, b, log=TRUE))
			return(-log.lk)
		}
			
			MLE <- optim(pars, LL, y=y, method="BFGS")
			
			return(MLE$par)
	   
	   }

	if (distname == "chisq" ) {
			dhhh <- mean(data)
			pars=c(dhhh)
			# Log of likelihood of a chi squared distribution
			LL <- function(pars, y){
			a <- pars[1]
			log.lk <- sum(dchisq(y, a, log=TRUE))
			return(-log.lk)
		}
	
			
			MLE <- optim(par=df, LL, y=y, method="BFGS")
			
			return(MLE$par)
	   
		}
	if (distname == "exp") {
			rate=1/m
			pars=c(rate)
			# Log of likelihood of a exponential distribution
			LL <- function(pars, y){
			a <- pars[1]
			log.lk <- sum(dexp(y, a,log=TRUE))
			return(-log.lk)
		}
			
			MLE <- optim(par=rate, LL, y=y, method="BFGS")
			
			return(MLE$par)
	   
			
		}
		
	if (distname == "binom" ) {
			
			#No of trials
			NofT=as.numeric()
			prob=as.numeric()
			
			#mean=NofT*prob & variance = NofT*prob(1-prob)
			
			prob=1-(v/m)
			NofT=m/prob
			pars=c(prob,NofT)
			# Log of likelihood of a binomial distribution
			LL <- function(pars, y){
			a <- pars[1]; b <- pars[2]
			log.lk <- sum(dbinom(y, a, b, log=TRUE))
			return(-log.lk)
		}
			#MLE <- optim(pars, LL, y=data, method="BFGS")
			mle=pars
			return(mle)
	   }
	
		
	if (distname == "gamma" ) {
			data <- data + 1e-6
			s = log(mean(data)) - (sum(log(data)))/length(data)
			estimated_alpha <- ((3 - s) + sqrt( ((s-3)**2) + (24*s) ))/(12*s)
			estimated_beta <- mean(data)/estimated_alpha
			return(c(estimated_alpha, estimated_beta))		
	   }
	   
	if (distname == "geom" ) {
			prob<-if (m>0) 1/(1+m)
					else 
					{NaN
					break}
			pars=c(prob)
			# Log of likelihood of a geometric distribution
			LL <- function(pars, y){
			a <- pars[1]
			log.lk <- sum(dgeom(y, a, log=TRUE))
			return(-log.lk)
		}
			
			
			MLE <- optim(par=prob, LL, y=y, method="BFGS")
			
			return(MLE)
			
			
	   }
	
}	
