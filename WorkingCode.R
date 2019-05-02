####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
BirdDat <- read.csv(text=getURL(gitstring))
####################################################################################################################################

#### MODEL 2 - Beta-Binomial ####
x = BirdDat$RedtailedHawk
ndat = length(x)
sumx = sum(x)
N = 10^4
alpha0 = 1
beta0 = 1
xgrid = seq(0.036,0.039,length = 1000)
plot(xgrid, dbeta(xgrid, sumx + alpha0, ndat*N + beta0 - sumx  ),
     type = "l", lty = 4, lwd=2, col = "black")

alpha = c(60, mean(x),1, 100)
beta = c(60,1,1, 1)

alpha.star =  sumx + alpha
beta.star = ndat*N + beta - sumx 

colors = c("red","blue","green", "yellow")
linetype = c(5,2,3,1)
plot(xgrid,xgrid,type="n",ylab="Posterior Density",
     main="Posterior Beta Distributions", las=1, xlab = expression(theta), ylim = c(0,2100))
for(i in 1:length(alpha)){
  lines(xgrid, dbeta(xgrid, alpha.star[i], beta.star[i]), 
        type = 'l', col=colors[i],lwd=2, lty = linetype[i])
}
legend("topleft", legend=c("alpha = beta = 60", "alpha = 376.1, beta = 1", "alpha = beta = 1", "alpha = 100, beta = 1"),
       lwd=rep(2,5), col=colors, lty = linetype, bty="n", ncol=1)
alpha
beta

#######
N = 10^6
# the log-likelihood is a function of lambda and the data
sterlings_logf = function(n){
  n * log(n) - n
}


logchoose = function(N, x){
  return(sterlings_logf(N) - sterlings_logf(x) - sterlings_logf(N-x))
}



poissonLogLik = function(x,theta){
  # total number of observations
  n = length(x)
  xlarge = x[x>=170]
  xsmall = x[x<=170]
  logfact_small = log(factorial(xsmall))
  logfact_large = xlarge * log(xlarge) - xlarge
  
  # equation
  # logLik = (sum(x)*log(theta)-n*theta - (sum(logfact_small) + sum(logfact_large)))
  logLik = (sum(x)*log(theta)-n*theta - (sum(sterlings_logf(x))))
  return(logLik)
}

pois_ll = poissonLogLik(x, mean(x))

BIC_M1 = -2 * pois_ll - log(ndat)

binomLogLik = function(x, theta){
  n = length(x)
  xlarge = x[x>=150]
  xsmall = x[x<=150]
  lchoose_small = log(choose(N, xsmall))
  lchoose_large = logchoose(N,xlarge)
  
  logLik = sum(x) * log(theta) + (n*N - sum(x)) * log(1 - theta) + (sum(logchoose(N,x)))
  return(logLik)
}

bin.mle = sum(x)/(length(x) * N)
bin_ll = binomLogLik(x, bin.mle)

BIC_M2 = -2 * bin_ll + log(length(x))
print(c(BIC_M1, BIC_M2))

