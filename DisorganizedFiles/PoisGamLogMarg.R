library(LearnBayes)
data(soccergoals)
logpoissgamma=function(theta,datapar)
{
  y=datapar$data
  npar=datapar$par
  lambda=exp(theta)
  loglike=log(dgamma(lambda,shape=sum(y)+1,rate=length(y)))
  logprior=log(dgamma(lambda,shape=npar[1],rate=npar[2])*lambda)
  return(loglike+logprior)
}
datapar=list(data=x,par=c(alpha1,beta1))
ggg=seq(-100, 100, length.out = 10000)
for(i in 1:length(ggg)){
  if(logpoissgamma(ggg[i], datapar) != -Inf)
  print(paste(i, ggg[i]))
}


fit1=laplace(logpoissgamma,6.19,datapar)
(logmarg1 = fit1$int)
