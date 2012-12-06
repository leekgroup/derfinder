## getTstats() - this is just a less complex version of eBayes().
## most of this code can be found in limma:::fitFDist
## arguments:
## --fit:  list with elements $sigma, $df.residual, $Amean (if trend=TRUE), $coefficients, and $stdev.unscaled
## --fit (cont): this is usually the output of getLimmaInput (in our pipeline), which is designed to be a more memory-efficient version of lmFit().
## --trend: as in eBayes, if TRUE, allow an intensity trend in the priors of variances (across genes/bps)
## return:
## list with elements $tt (containing the moderated t-statistics for each gene/bp) and $logfchange (containing the log2 fold change in coverage between the groups, as estimated by a linear model)


getTstats = 
function(fit, trend = FALSE){
	require(limma)
	# get d0 and s02 (prior parameters)
	sg2 <- fit$sigma^2
	sg2 <- pmax(sg2,1e-05*median(sg2)) 
	zg <- log(sg2)
	eg <- zg - digamma(fit$df.residual/2)+log(fit$df.residual/2)
	rm(zg);gc();gc()
	if(!trend){
	      ebar <- mean(eg)
	      G <- length(fit$sigma)
      	  d0 <- 2*trigammaInverse(mean((eg-ebar)^2*(G/(G-1)) - trigamma(fit$df.residual/2)))
	      s02 <- exp(ebar+digamma(d0/2)-log(d0/2))
  	}
	if(trend){
	      require(splines)
      	  design <- try(ns(fit$Amean, df = 4, intercept = TRUE), silent = TRUE)
          if (is(design, "try-error")) stop("Problem with Amean covariate; perhaps too few distinct values")
    	  fit1 <- lm.fit(design, eg)
      	  ebar <- fit1$fitted
      	  evar <- mean(fit1$residuals[-(1:fit1$rank)]^2)
      	  rm(fit1);gc();gc()
      	  evar <- evar - mean(trigamma(fit$df.residual/2))
      	  if(evar > 0){
      	    d0 <- 2*trigammaInverse(evar)
      	    s02 <- exp(ebar+digamma(d0/2)-log(d0/2))
      	  }
      	  if(evar <= 0){
      	  	d0 = Inf
      	  	s02 = exp(ebar)
      	  }
      	  rm(ebar,design);gc();gc()
    }
	if(!is.finite(d0)){
		sgtilde <- s02
	}else{sgtilde <- (d0*s02+fit$df.residual*sg2)/(d0+fit$df.residual)}
	tt <- fit$coefficients/(sqrt(sgtilde)*fit$stdev.unscaled)
	logfchange = fit$coefficients
	return(list(tt=tt,logfchange=logfchange))
}