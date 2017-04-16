library(statnet)

# Code for the Coverage Probability Barplots. The Bootstrap CI's always seem to cover the 'true' coefficients


##################################
# Faux_mesa_high, 205 nodes
##################################
#qw<-ergm(faux.magnolia.high ~ edges +nodematch("Sex")+esp(0))
#summary(qw)

# get data
data(faux.magnolia.high)

mple <- try(ergm(faux.magnolia.high ~ edges +nodematch("Sex")+gwesp(0.5,fixed=TRUE), estimate = "MPLE" )) #
summary(mple)

#main function
bootSE.faux <- function( erg.obj, simulation, iters){
  
  # create empty matrices for logistic coverage
  logistic.coverage<- matrix(0,3, iters)
  rownames(logistic.coverage)<- c("edges", "sex", "gwesp")
  
  # create empty matrices for bootstrap coverage
  bootstrap.coverage<- matrix(0,3,iters)
  rownames(bootstrap.coverage)<- c("edges", "sex", "gwesp")
  
  # create empty matrices for estimated MPLE
  mple.sim<- matrix(0,3,simulation*iters)
  rownames(mple.sim)<- c("edges", "sex", "gwesp")
  
  # create empty matrices for estimated MCMCMLE
  mcmcmle.sim<- matrix(0,3,iters)
  rownames(mcmcmle.sim)<- c("edges", "sex", "gwesp")
  
  # create empty matrices for logistic coverage
  mcmcmle.coverage<- matrix(0,3, iters)
  rownames(mcmcmle.coverage)<- c("edges", "sex", "gwesp")
  
  # create vector for density
  density <- c()
  
  net.mple <- simulate.ergm(mple, nsim=iters, 
             control=control.simulate.ergm(MCMC.burnin=5000000, MCMC.interval=3000000)) 
  
  for(j in 1:iters){ # loop for number of iterations
    
    print(j)
    
    #create empty matrices for bootstrap estimates
    coef <- matrix(0, 3, simulation)
    rownames(coef) <- c("edges", "sex", "gwesp")
    Std.err <- matrix(0,3,simulation)
    rownames(Std.err) <- c("edges", "sex", "gwesp")
    
    #simulate network from true values
    # simulate network
    
    #estimate coefficients
    sim.mple <- ergm(net.mple[[j]]~edges+ nodematch("Sex")+gwesp(0.5,fixed=TRUE), estimate= "MPLE")
    
    #############
    # MCMCMLE
    ############
    
    bootest.mcmcmle <- try(ergm(net.mple[[j]]~edges+ nodematch("Sex")+gwesp(0.5,fixed=TRUE))) #estimate coefficients
    
    # save values
    if(summary(bootest.mcmcmle)[2]=="try-error"){ # in case the ERGM degenerates => safe NA
      mcmcmle.sim[,j] <-NA
    }
    else{ # otherwise safe estimates
      mcmcmle.sim[,j] <-coef(bootest.mcmcmle)  
    } 
    
    # calculate MCMCMLE CI for each simulation
    
    # create empty matrix to save values
    mcmcmle.CI <- matrix(0,3,2)
    colnames(mcmcmle.CI)<- c("2.5%", "97.5%")
    rownames(mcmcmle.CI)<- c("edges", "sex", "gwesp")
    
    mcmcmle.CI[1,]<- confint(bootest.mcmcmle, level=0.95)[1,] # 95% CI for edges
    mcmcmle.CI[2,]<- confint(bootest.mcmcmle, level=0.95)[2,] # 95% CI for Sex
    mcmcmle.CI[3,]<- confint(bootest.mcmcmle, level=0.95)[3,] # 95% CI for gwesp
    
    # determine whether CI includes the actual value => we get a TRUE/FALSE vector
    true.false.mcmcmle<-cbind( mcmcmle.CI[,2] >= coef(mple) & mcmcmle.CI[,1] <= coef(mple) ) # TRUE/FALSE
    
    # turn TRUE/FALSE vector into a vector of 1 and 0
    for(m in 1:3){ # save 1 if coef(mple) lies in logistic CI
      if(true.false.mcmcmle[m,1]==TRUE){
        mcmcmle.coverage[m,j]=1
        
      } # end if
    } # end for loop m
    
    #########################################
    # calculate logistic CI for each simulation
    ##########################################
    
    
    # create empty matrix to save values
    log.CI <- matrix(0,3,2)
    colnames(log.CI)<- c("2.5%", "97.5%")
    rownames(log.CI)<- c("edges", "sex", "gwesp")
    
    log.CI[1,]<- confint(sim.mple, level=0.95)[1,] # 95% CI for edges
    log.CI[2,]<- confint(sim.mple, level=0.95)[2,] # 95% CI for Sex
    log.CI[3,]<- confint(sim.mple, level=0.95)[3,] # 95% CI for gwesp
    
    # determine whether CI includes the actual value => we get a TRUE/FALSE vector
    true.false.log<-cbind( log.CI[,2] >= coef(mple) & log.CI[,1] <= coef(mple) ) # TRUE/FALSE
    
    # turn TRUE/FALSE vector into a vector of 1 and 0
    for(m in 1:3){ # save 1 if coef(mple) lies in logistic CI
      if(true.false.log[m,1]==TRUE){
        logistic.coverage[m,j]=1
      }}
    
    net<- simulate.ergm(sim.mple, nsim=simulation,
          control=control.simulate.ergm(MCMC.burnin=5000000, MCMC.interval=3000000))
    
    # parametric bootstrap
    for ( i in 1:simulation){  # loop for numbers of simulations in each iteration
      
      
      #save(net, file=paste((j-1)*simulation+i, '_sim_net.Rdata', sep='') )
      
      # calculate density of network
      
      density[(j-1)*simulation+i]<- network.density(net[[i]])
      
      #estimate coefficients
      bootest <- try(ergm(net[[i]]~edges+ nodematch("Sex")+gwesp(0.5,fixed=TRUE), estimate= "MPLE"))
      
      # save values
      if(summary(bootest)[2]=="try-error"){ # in case the ERGM degenerates => safe NA
        mple.sim[,(j-1)*simulation+i] <-NA
        coef[,i]<- NA
        #Std.err[,i] <- NA
      }
      else{ # otherwise safe estimates
        mple.sim[,(j-1)*simulation+i] <-summary(bootest)$coef[,1]   
        coef[,i]<-coef(bootest)
        #Std.err[,i] <- summary(bootest)$coef[,2]
      }
      
      
      
      
    } # end for loop i
    
    
    
    #use the coefficient estimates from the inner loop (loop i) to calculate the bootstrap intervals
    
    # save coef and Std.err in list
    theta.hat <- list()
    theta.hat[[1]] <- coef
    theta.hat[[2]] <- Std.err
    
    names(theta.hat)=c("Bootstrap estimates","Bootstrap Standard Error") # name list
    
    # bootstrap.CI: if coef(mple) is in BootCI, save 1 otherwise 0
    # create empty matrix
    bootstrap.CI <- matrix(0,3,2)
    colnames(bootstrap.CI)<- c("upper", "lower")
    rownames(bootstrap.CI)<- c("edges", "sex", "gwesp")
    
    # calculate the bootstrap intervals
    bootstrap.CI[1,]<-quantile(theta.hat[[1]][1,], probs = c(0.975, 0.025)) # edges
    bootstrap.CI[2,]<-quantile(theta.hat[[1]][2,], probs = c(0.975, 0.025)) # Sex
    bootstrap.CI[3,]<-quantile(theta.hat[[1]][3,], probs = c(0.975, 0.025)) # gwesp
    
    # determine whether CI includes the actual value => we get a TRUE/FALSE vector
    true.false.boot<-cbind( bootstrap.CI[,1] >= coef(mple) & bootstrap.CI[,2] <= coef(mple) ) # TRUE/FALSE
    
    # turn TRUE/FALSE vector into a vector of 1 and 0
    for(k in 1:3){ # save 1 if coef(mple) lies in boot CI
      if(true.false.boot[k,1]==TRUE){
        bootstrap.coverage[k,j]=1
      } # end if
    } # end for loop k
    
    
    
  } #end for loop j
  
  # save results in a list and return list
  result<- list()
  result[[1]]<-mple.sim
  result[[2]]<-bootstrap.coverage
  result[[3]]<-logistic.coverage
  result[[4]]<-mcmcmle.sim
  result[[5]]<-mcmcmle.coverage
  result[[6]]<-density
  
  names(result)=c("Simulated MPLE Estimates","Bootstrap coverage","Logistic Coverage", "Simulated MLE Estimates", "MCMCMLE Coverage", "Density")
  
  return(result)
}



# simulation MPLE
set.seed(5555)
sim=100
iter=20
faux.results.mple2 <- bootSE.faux(mple, sim,iter)

# coef=Inf appears
faux.results.mple[[1]][!is.finite(faux.results.mple[[1]])] <- NA
faux.results.mple[[2]][!is.finite(faux.results.mple[[2]])] <- NA
faux.results.mple[[3]][!is.finite(faux.results.mple[[3]])] <- NA
faux.results.mple[[4]][!is.finite(faux.results.mple[[4]])] <- NA
faux.results.mple[[5]][!is.finite(faux.results.mple[[5]])] <- NA


## Coverage probability plot

# create empty matrix for barplot values
par(mfrow=c(1,1))
counts<- matrix(0,3,3)
colnames(counts)<- c("Edges", "Sex", "Gwesp")
rownames(counts)<- c("Bootstrap", "MCMCMLE", "Logistic")

# fill matrix with corresponding values
counts[1,1] <- mean(faux.results.mple[[2]][1,]) # edges bootstrap
counts[1,2] <- mean(faux.results.mple[[2]][2,]) # sex bootstrap
counts[1,3] <- mean(faux.results.mple[[2]][3,]) # gwesp bootstrap

counts[3,1] <- mean(faux.results.mple[[3]][1,]) # edges logistic
counts[3,2] <- mean(faux.results.mple[[3]][2,]) # sex logistic
counts[3,3] <- mean(faux.results.mple[[3]][3,]) # gwesp logistic

counts[2,1] <- mean(faux.results.mple[[5]][1,]) # edges mcmcmle
counts[2,2] <- mean(faux.results.mple[[5]][2,]) # sex mcmcmle
counts[2,3] <- mean(faux.results.mple[[5]][3,]) # gwesp mcmcmle

# Plot
barplot(counts, beside=TRUE, main="Coverage Probability Faux Magnolia High", ylim=c(0,1),legend = rownames(counts),args.legend = list(x="bottomleft", cex=1.5),
         cex=1.5, cex.axis=1.5,cex.names = 2, cex.main=2)
abline(h=0.95,col="grey", lty=2, lwd=2)

# Histogram of the density of simulated networks
hist(faux.results.mple[[6]], main="Density of Simulated Networks", xlab="Density")
dens.faux.magnolia=network.density(faux.magnolia.high)
abline(v=dens.faux.magnolia,col="red")

## Bias plot

boxplot(faux.results.mple[[1]][1,]-mple$coef[1], faux.results.mple[[1]][2,]-mple$coef[2], 
        faux.results.mple[[1]][3,]-mple$coef[3],
        names=c("edges", "Sex" ,"GWESP"), main="MPLE", ylim=c(-0.4,0.4))
abline(a=0, b=0, col="lightgrey", lty=3)

## Density time line

plot(faux.results.mple[[6]], main="Density of Simulated Networks", ylab="Density")


## Histogram of Bootstrap iterations

for(i in 1:10){
  par(mfrow=c(3,1),oma=c(0,0,2,0))
  
  #edges
  hist(faux.results.mple[[1]][1,(100*(i-1)+1):(i*100)], main="Edges", xlab="")
  abline(v=mple$coef[1],col="red")
  
  #Sex
  hist(faux.results.mple[[1]][2,(100*(i-1)+1):(i*100)], main="Sex", xlab="")
  abline(v=mple$coef[2],col="red")
  
  #gwesp
  hist(faux.results.mple[[1]][3,(100*(i-1)+1):(i*100)], main="Gwesp", xlab="")
  abline(v=mple$coef[3],col="red")
  
  title(paste("Bootstrap Iteration",i), outer=TRUE)
}


#######################################
#######################################
## system time MCMCMLE vs MPLE

bootMPLE <- function(net, simulation){
  for ( j in 1:simulation){  # loop for numbers of simulations in each iteration
    
    #estimate coefficients
    bootest <- try(ergm(net[[j]]~edges+ nodematch("Sex")+gwesp(0.5,fixed=TRUE), estimate= "MPLE"))
    
  }
}


mple <- try(ergm(faux.magnolia.high ~ edges +nodematch("Sex")+gwesp(0.5,fixed=TRUE), estimate = "MPLE" )) #
summary(mple)

# empty vectors to store system times

mcmle<- c()
boot.mple<- c()
simulation=500
iter=100

set.seed(5555)
net.mple <- simulate.ergm(mple, nsim=iter )

for ( i in 1:iter){
  print(i)
  
  # mcmcmle system time
  a<- system.time(try(ergm(net.mple[[i]]~edges+ nodematch("Sex")+gwesp(0.5,fixed=TRUE))))
  mcmle[i]<- a[3]
  
  #### boot MPLE
  
  #estimate coefficients
  sim.mple <- ergm(net.mple[[i]]~edges+ nodematch("Sex")+gwesp(0.5,fixed=TRUE), estimate= "MPLE")
  
  # simulate.networks
  net<- simulate.ergm(sim.mple, nsim=simulation)
  
  # save computing time
  b<- system.time(bootMPLE(net, simulation))
  boot.mple[i]<- b[3]
  
}



# median run time for MCMCMLE
time.mle<- c()
time.simulation<- c()


for (i in 41:50){
  print(i)
  set.seed(i*1000)
  b<- system.time(ergm(faux.magnolia.high ~ edges +nodematch("Sex")+gwesp(0.25,fixed=TRUE), control=control.ergm(MCMLE.maxit=100,
                  MCMC.samplesize = 8000, MCMC.interval=5000), eval.loglik=FALSE))
  time.mle[i]<- b[3]
  
  #get numbers of iterations
  set.seed(i*1000)
  model<- ergm(faux.magnolia.high ~ edges +nodematch("Sex")+gwesp(0.25,fixed=TRUE), control=control.ergm(MCMLE.maxit=100,
               MCMC.samplesize = 8000, MCMC.interval=5000 ), eval.loglik=FALSE)
  time.simulation[i]<-b[3]/(model$iterations)
  
}
median(time.mle) # 119.07
median(time.simulation) #39.83 => need 500 networks for bootMPLE => 39.83/16=2.49

# median run time for MPLE
time.mple<- c()

set.seed(123)
for (i in 1:100){
  print(i)
  b<- system.time(ergm(faux.magnolia.high ~ edges +nodematch("Sex")+gwesp(0.5,fixed=TRUE), 
                       control=control.ergm(MCMLE.maxit=100), estimate="MPLE"))
  time.mple[i]<-b[3]
}
median(time.mple) # 0.51

# time for simulating network 

# simulating 1024 networks
sim.time<- system.time(simulate(mple, nsim=1024)) # 42.08

