#GOODNESS OF FIT PARAMETRIC BOOTSTRAP
gof_test<-function(x,dist,nboot=1000){
	if (dist=="bern"){
		qdist <- get(paste("q", "binom", sep=""), mode = "function")
		rdist <- get(paste("r", "binom", sep=""), mode = "function")
	}
	else{
		qdist <- get(paste("q", dist, sep=""), mode = "function")
		rdist <- get(paste("r", dist, sep=""), mode = "function")
	}
	
	
	n<- length(x)
	input=x
	D.vec <- NULL
	est_par <- NULL
	
	if (dist %in% c("geom", "pois", "exp", "chisq")) {
      par <- mle(x,dist)
	  y <- qdist(c(1:n)/(n+1) ,par)
      D0 <- kst(input, y)
    } 
    else {
      
      par <- mle(x,dist)
      
      
      # Only for Bernoulli distribution
      if (dist %in% c("bern")) {
        par[2] = par[1]
        par[1] = 1
		y <- qdist(c(1:n)/(n+1) ,par[1], par[2])
		D0 <- kst(input, y)
		}
	  
	  if (dist %in% c("binom")) {
        par[1] = floor(par[1])
		y <- qdist(c(1:n)/(n+1) ,par[1], par[2])
		D0 <- kst(input, y)
		}
	  else{
      y <- qdist(c(1:n)/(n+1) ,par[1], par[2])
      D0 <- kst(input, y)
      
	  }
    }
  
	
	
	for(i in 1:nboot) {
    
    if (dist %in% c("geom", "pois", "exp", "chisq")) {
      xstar=rdist(n ,par)
	  par <- mle(xstar,dist)
	  y <- qdist(c(1:n)/(n+1) ,par)
      D.vec <- c(D.vec, kst(xstar, y))
      est_par <- c(est_par, c(par))
      input <- y
    } 
    else {
      
      par <- mle(x,dist)
      
      
      # Only for Bernoulli distribution
      if (dist %in% c("bern")) {
        par[2] = par[1]
        par[1] = 1
		xstar=rdist(n ,par[1], par[2])
		y <- qdist(c(1:n)/(n+1) ,par[1], par[2])
		D.vec <- c(D.vec, kst(xstar, y))
		est_par <- c(est_par, c(par[1], par[2]))
		input <- y
      }
	  
	  if (dist %in% c("binom")) {
        par[1] = floor(par[1])
		xstar=rdist(n ,par[1], par[2])
		y <- qdist(c(1:n)/(n+1) ,par[1], par[2])
		D.vec <- c(D.vec, kst(xstar, y))
		est_par <- c(est_par, c(par[1], par[2]))
		input <- y
	  }
	  else{
	  xstar=rdist(n ,par[1], par[2])
      y <- qdist(c(1:n)/(n+1) ,par[1], par[2])
      D.vec <- c(D.vec, kst(xstar, y))
      est_par <- c(est_par, c(par[1], par[2]))
      input <- y
	  }
    }
  }
  #print(D.vec)
  p_value <- sum(D.vec > D0)/length(D.vec)
  return(p_value)
  
}