set.seed(123) #make sure it's reproducable

#PART 1
alb.data = read.csv("STAL data.csv", header=T) #breeding pair data

#Calculates NLL assuming exponential growth
#Input: growth rate (r), initial pop (N0), additional variance (CVadd), real data (alb.data)
#Return: NLL
getNLL = function(r, N0, CVadd, alb.data) {
  
  n.years = nrow(alb.data)
  real.pairs = alb.data$BreedingPairs #actual counts
  est.pairs = numeric(length = n.years) #estimated counts
  years = alb.data$Year
  first.year = years[1]
  
  #fill in estimated data with exponential model
  est.pairs[1] = N0
  for(i in 2:n.years) {
    est.pairs[i] = est.pairs[1]*(1+r)^(years[i]-first.year)
  }
  
  #make shortcuts for constant parts of NLL so equation isn't so long
  first = log(sqrt(0.02^2+CVadd^2))
  denom = 2*(0.02^2+CVadd^2)
  
  #Find NLL
  NLL = 0
  for(i in 1:length(est.pairs)) {
    numer = (log(real.pairs[i]/est.pairs[i]))^2
    NLL = NLL + (first+(numer/denom))
  }
  
  return(NLL)
}

getNLL(0.07, 10, 0.05, alb.data) #check if it's right

#PART 2

#Calculates negative log of prior on r
#Input: lower bound of r prior (r.prior.lower), upper bound of r prior (r.prior.upper)
#Return: NL of r prior
rNLprior = function(r.prior.lower, r.prior.upper) {
  rNL = -log(1/(r.prior.upper-r.prior.lower))
  return(rNL)
}

#Calculates negative log of prior on N0
#Input: lower bound of N0 prior (N0.prior.lower), upper bound of N0 prior (N0.prior.upper)
#Return: NL of N0 prior
N0NLprior = function(N0.prior.lower, N0.prior.upper) {
  N0NL = -log(1/(N0.prior.upper-N0.prior.lower))
  return(N0NL)
}

#Calculates negative log of prior on CVadd
#Input: lower bound of CVadd prior (CVadd.prior.lower), upper bound of CVadd prior (CVadd.prior.upper)
#Return: NL of CVadd prior
CVaddNLprior = function(CVadd.prior.lower, CVadd.prior.upper) {
  CVaddNL = -log(1/(CVadd.prior.upper-CVadd.prior.lower))
  return(CVaddNL)
}

NLpriors = rNLprior(-0.1,0.2) + N0NLprior(0,50) + CVaddNLprior(0,0.3) #store sum of NL priors

#PART 3

#Runs the Markov Chain Monte Carlo method. Takes initial parameters and bounces around
#until it settles in the region of high posterior probability.
#Input: number of draws (ndraws), file of counts (filename), initial growth rate (rinit),
#       initial count (N0), initial additional variance (CVaddinit)
#Return: matrix of accepted draws (r*, N0*, CVadd*, X*) where X* is total NLL
runMCMC = function(ndraws, filename, rinit, N0init, CVaddinit) {
  
  posterior = matrix(nrow = ndraws, ncol = 4) #matrix of accepted draws
 
  Xstarinit = getNLL(rinit, N0init, CVaddinit, alb.data) + NLpriors #initial Xstar
  posterior[1,] = c(rinit, N0init, CVaddinit, Xstarinit) #initial draw
  #loop through all the draws
  for(i in 2:ndraws) {
    rstar = posterior[(i-1), 1] + runif(1, -0.01, 0.01)
    N0star = posterior[(i-1), 2] + runif(1, -2, 2)
    CVaddstar = posterior[(i-1), 3] + runif(1, -0.05, 0.05)
    
    Xstar = getNLL(rstar, N0star, CVaddstar, alb.data) + NLpriors
    
    ratio = exp(posterior[(i-1), 4]-Xstar)
    randNum = runif(1, 0, 1)
    #if the random number is less that the ratio, accept draw
    if(randNum < ratio) {
      posterior[i,] = c(rstar, N0star, CVaddstar, Xstar)
    } else { #reject draw and start from same spot as before
      posterior[i,] = posterior[i-1,]
    }
  }
  return(posterior)
}

#PART 4
posterior = runMCMC(100000, "STAL data.csv", 0.03, 10, 0.05)

#thin the posterior to reduce autocorrelation and make it a manageable size
temp = posterior[20000:100000,]
thin.Post = temp[seq(1,80000,40),]

#PART 5

#find distributions for pairs in 2014, individuals, and delta bw 2014 & 2015
dist.pairs = thin.Post[,2]*(1+thin.Post[,1])^(2014-1954)
dist.pop = 7.1*dist.pairs
dist.diff = dist.pop*thin.Post[,1]

#STUFF FOR RESULTS SECTION

#quantiles for parameters
r.quant = quantile(thin.Post[,1], probs = c(0.025, 0.5, 0.975))
N0.quant = quantile(thin.Post[,2], probs = c(0.025, 0.5, 0.975))
CVadd.quant = quantile(thin.Post[,3], probs = c(0.025, 0.5, 0.975))

#histograms for parameters
hist(thin.Post[,1], main="Growth Rate", xlab="r")
hist(thin.Post[,2], main="Initial Pairs", xlab="N0")
hist(thin.Post[,3], main="Additional Variance", xlab="CVadd")

#quantiles for pairs in 2014, individuals, and delta bw 2014 & 2015
pairs.quant = quantile(dist.pairs, probs=c(0.025, 0.5, 0.975))
pop.quant = quantile(dist.pop, probs=c(0.025, 0.5, 0.975))
diff.quant = quantile(dist.diff, probs=c(0.025, 0.5, 0.975))

