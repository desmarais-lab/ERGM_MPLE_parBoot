library(statnet)
#memory.limit(size=10000)

model1<- ergm(scnetwork~edges)
summary(model1)


model2<- ergm(scnetwork~edges+ nodecov('salience'),  estimate= "MPLE")
summary(model2)


model3<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year'),  estimate= "MPLE")
summary(model3)

model4<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2),  estimate= "MPLE")
summary(model4)

# by adding triangle, the simulation takes a lot of more time (about 15mins)
model5<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+triangle,  estimate= "MPLE"))# 1026.4
summary(model5) #degenerate

model51<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+kstar(3)+triangle,  estimate= "MPLE"))# 105.94
summary(model51) #degenerate

model52<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+triangle,  estimate= "MPLE"))# 896 
summary(model52) #degenerate

model53<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2),  estimate= "MPLE"))# 181
summary(model53)

model54<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+kstar(3),  estimate= "MPLE"))# 75.31
summary(model54) # degenerate

#lets see how long it takes to estimate the MCMCMLE 

model6<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2))) # 1019.72 did not converge after 20 iterations

# increase max number of iterations
model611<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2), control=control.ergm(MCMLE.maxit = 30))) #1284.66

model612<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2), control=control.ergm(MCMLE.maxit = 50)))
# still not converging



model61<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+triangle)) # 2260.84
# MCMCMLE is not converging

model62<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+kstar(3)+triangle))# 987.28

# try gwesp

model6<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+kstar(3)+ gwesp(1, fixed=TRUE), estimate = "MPLE")
# degenerate

model6<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+kstar(3)+ gwesp(0, fixed=TRUE), estimate = "MPLE")
# degenerate

model6<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+kstar(2)+ gwesp(0, fixed=TRUE), estimate = "MPLE")
# degenerate


# include absdiff(year)

model7<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year') +kstar(2), estimate = "MPLE")
# degenerate
#Warning messages:
  #1: In ergm.pl(Clist = Clist, Clist.miss = Clist.miss, m = m, theta.offset = init,  :
  #   Too many unique dyads. MPLE is approximate, and MPLE standard errors are suspect.
  #2: In ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init,  :
  #   non-integer #successes in a binomial glm!
  #3: In ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init,  :
  #   non-integer #successes in a binomial glm!
  #4: In eval(expr, envir, enclos) :
  #   non-integer #successes in a binomial glm!

# include gwdegree instead of kstar(2)

model8<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year') +gwdegree(1,fixed=TRUE), estimate = "MPLE")
# degenerate
 #Warning messages:
#  1: In ergm.pl(Clist = Clist, Clist.miss = Clist.miss, m = m, theta.offset = init,  
#     Too many unique dyads. MPLE is approximate, and MPLE standard errors are suspect.
#  2: In ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init,  :
#     non-integer #successes in a binomial glm!
#  3: In ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init,  :
#     non-integer #successes in a binomial glm!
#  4: In eval(expr, envir, enclos) :
#     non-integer #successes in a binomial glm!


model9<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ gwdegree(1,fixed=TRUE), estimate = "MPLE")
#Warning messages:
  #1: In ergm.pl(Clist = Clist, Clist.miss = Clist.miss, m = m, theta.offset = init,  :
  #   Too many unique dyads. MPLE is approximate, and MPLE standard errors are suspect.
  #2: In ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init,  :
  #    non-integer #successes in a binomial glm!
  #3: In ergm.mple(Clist, Clist.miss, m, MPLEtype = MPLEtype, init = init,  :
  #    non-integer #successes in a binomial glm!
  #4: In eval(expr, envir, enclos) :
  #   non-integer #successes in a binomial glm!
  #5: In dev.resids(y, mu, weights) :
  #   Reached total allocation of 8107Mb: see help(memory.size)
  #6: In dev.resids(y, mu, weights) :
  #   Reached total allocation of 8107Mb: see help(memory.size)
  #7: In mu.eta(eta) :
  #   Reached total allocation of 8107Mb: see help(memory.size)
  #8: In mu.eta(eta) :
  #   Reached total allocation of 8107Mb: see help(memory.size)

model10<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ gwdegree(0,fixed=TRUE), estimate = "MPLE")
# works

# 
set.seed(1)
model11<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), estimate = "MPLE")
summary(model11) # works
sys.time.model11<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), estimate = "MPLE"))
sys.time.model11


# try again with more iterations
set.seed(1)
model11.mle<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), 
                   control=control.ergm(MCMLE.maxit = 250,MCMC.samplesize=32784, init=c(-24.55, 1.438, 0.004732, -0.03458, -3.457 ),
                                        MCMLE.MCMC.precision=0.01))
summary(model11.mle) # didn't converge after 1000 iterations
sys.time.model11.mle<- system.time(ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), control=control.ergm(MCMLE.maxit = 100)))
sys.time.model11.mle



model12<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(1,fixed=TRUE), estimate = "MPLE")
# degenerates

model13<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE)+ gwesp(0,fixed=TRUE), estimate = "MPLE")
summary(model13) # degenerates

model14<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(1,fixed=TRUE)+ gwesp(1,fixed=TRUE), estimate = "MPLE")
summary(model14) # degenerate

# try isolates
set.seed(1)
model15<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ 
                 gwdegree(0,fixed=TRUE)+isolates, estimate = "MPLE")
summary(model15) # isolates NA

set.seed(1)
model16<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ 
                 gwdegree(1,fixed=TRUE)+isolates, estimate = "MPLE")
# degenerates

set.seed(1)
model17<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ 
                 gwesp(0,fixed=TRUE)+isolates, estimate = "MPLE")
# degenerates

set.seed(1)
model18<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ 
                 gwesp(1,fixed=TRUE)+isolates, estimate = "MPLE")
# degenerate

set.seed(1)
model19<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE)+ triadcensus, estimate = "MPLE")
summary(model11) # works

set.seed(1)
model20.mle<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE),verbose=2, 
                   control=control.ergm(MCMLE.maxit = 250,MCMC.samplesize=8192,  MCMLE.termination = "Hotelling"))

# increase interval
set.seed(1)
model20.mle<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), 
                   control=control.ergm(MCMLE.maxit = 150,MCMC.interval=8192, MCMLE.termination = "Hotelling"))

###########################################
### standardize year
###########################################


v<-get.vertex.attribute(scnetwork, "year")
m<- mean(v)
s<- sd(v)
v.standard<- (v-m)/s

scnetwork<-set.vertex.attribute(scnetwork,"year.center", v.standard)

set.seed(1)
model21.mle<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year.center')+ absdiff('year.center')+ gwdegree(0,fixed=FALSE), 
                   control=control.ergm(MCMC.interval=8192, MCMLE.termination = "Hotelling"))

# year.center.square

v.standard<- get.vertex.attribute(scnetwork, "year.center")
v.standard.square<- v.standard^2
scnetwork<-set.vertex.attribute(scnetwork,"year.center.square", v.standard.square)

model20.mle<- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year.center') + absdiff('year.center') + absdiff('year.center.square') + dsp(0) + esp(0)+isolates,
                   estimate = "MPLE") 

## smaller network

model54.1<- ergm(scnetwork54~edges+ nodecov('salience')+ nodecov('year.center') + absdiff('year.center') + 
                   absdiff('year.center.square') + dsp(0) + esp(0)+isolates, estimate = "MPLE") 






##### parallel

library(statnet)
library(doParallel)
library(foreach)


registerDoParallel(cores=2)


# get data
load(file="supreme.RData")

mple <- ergm(scnetwork~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), estimate = "MPLE")
summary(mple)

#main function
boot.single.iteration.5.core <- function( erg.obj, simulation){ #simulation must be dividable by 5
  
  
  #create empty matrices for bootstrap estimates
  coef <- matrix(0, 5, simulation)
  rownames(coef) <- c("edges", "salience", "year", "absdiff", "gwdegree")
  Std.err <- matrix(0,5,simulation)
  rownames(Std.err) <- c("edges", "salience", "year", "absdiff", "gwdegree")
  
  #simulate network from true values
  net.mple <- simulate.ergm(mple, nsim=1)#, control=control.simulate.ergm(MCMC.burnin=5000000, MCMC.interval=100) ) # simulate network
  
  #estimate coefficients
  sim.mple <- ergm(net.mple~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), estimate = "MPLE")
  
  # loop for numbers of simulations in each iteration
  
  
  coef1<-foreach(i=1:simulation, .combine=cbind)%dopar%{
    
    net<- simulate.ergm(sim.mple, nsim=1) # changed it to 1
    library(ergm)
    #estimate coefficients
    bootest <- ergm(net~edges+ nodecov('salience')+ nodecov('year')+ absdiff('year')+ gwdegree(0,fixed=TRUE), estimate = "MPLE")
    coef(bootest)
    
  }
  
  
  # bootstrap.CI: if coef(mple) is in BootCI, save 1 otherwise 0
  # create empty matrix
  #bootstrap.CI <- matrix(0,4,2)
  #colnames(bootstrap.CI)<- c("upper", "lower")
  #rownames(bootstrap.CI)<- c("edges", "salience", "year", "kstar2")
  
  # calculate the bootstrap intervals
  #bootstrap.CI[1,]<-quantile(coef[1,], probs = c(0.975, 0.025)) # edges
  #bootstrap.CI[2,]<-quantile(coef[2,], probs = c(0.975, 0.025)) # salience
  #bootstrap.CI[3,]<-quantile(coef[3,], probs = c(0.975, 0.025)) # year
  #bootstrap.CI[4,]<-quantile(coef[4,], probs = c(0.975, 0.025)) # kstar2
  
  
  return(coef1)
}


set.seed(455)
boot <- boot.single.iteration.5.core(mple,20)

ret.list<- list(boot, mple)

save(ret.list,file="coverage_parallel_small1.RData")


library(plyr)
count(degree(scnetwork, gmode = "graph"))
