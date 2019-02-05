#Function for ttest
	mytfunc<-function(x){

	#Finding mean of the input
	xbar<-mean(x)

	#Finding standard deviation of the input
	sd0<-sd(x)

	#Finding length of the input
	n<-length(x)

	#Finding the t statistic using the formula mean*(sqrt(length of data)/(standard deviation)
	tstat<-xbar/(sd0/sqrt(n))

	#Finding the p-value using pt function of student t distribution
	p0<-pt(tstat,n-1)

	#if ((is.na(p0))|(!is.finite(p0))){ 
	#	p0<-pt(tstat, n-1, lower.tail = FALSE)
	#	if ((is.na(p0))|(!is.finite(p0))){
	#		p0=2 * pt(-abs(tstat), n-1)
	#	}
	#}

	#Returning the p-value & 1-pvalue for overexpressed and underexpressed genes.
	c(p0,1-p0)
	}
	
