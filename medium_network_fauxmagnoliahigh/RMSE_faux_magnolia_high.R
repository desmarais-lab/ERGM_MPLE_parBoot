library(statnet)

###################################
###################################
# using faux.magnolia.high data
###################################
###################################

data("faux.magnolia.high")

plot(faux.magnolia.high)

set.seed(123)
model=ergm(faux.magnolia.high~edges+nodefactor("Sex")+gwesp(0,fixed=TRUE))
summary(model)
#get true coefficients
coefficient <- summary(model)$coef[,1]

# function to calculate RMSE for MPLE, MCMCMLE with samplesize= 10000, 3000, 1500, 800, 400, 200, 100, 50
main.function.faux.magnolia.RMSE <- function(iters,  coefs=coefficient){
  # size = number of nodes in network
  
  set.seed(999)
  
  l=length(coefs)
  
  #creating empty matrices to store results
  coef.matrix.mple <- matrix(0,l,iters)
  se.matrix.mple <- matrix(0,l,iters)
  
  coef.matrix.10000 <- matrix(0,l,iters)
  se.matrix.10000 <- matrix(0,l,iters)
  
  coef.matrix.3000 <- matrix(0,l,iters)
  se.matrix.3000 <- matrix(0,l,iters)
  
  coef.matrix.1500 <- matrix(0,l,iters)
  se.matrix.1500 <- matrix(0,l,iters)

  coef.matrix.800 <- matrix(0,l,iters)
  se.matrix.800 <- matrix(0,l,iters)
  
  coef.matrix.400 <- matrix(0,l,iters)
  se.matrix.400 <- matrix(0,l,iters)
  
  coef.matrix.200 <- matrix(0,l,iters)
  se.matrix.200 <- matrix(0,l,iters)
  
  coef.matrix.100 <- matrix(0,l,iters)
  se.matrix.100 <- matrix(0,l,iters)
  
  coef.matrix.50 <- matrix(0,l,iters)
  se.matrix.50 <- matrix(0,l,iters)
  

  # density vector
  dens <- c()
  
  
  # start loop
  for(i in 1:iters){
    print(i)
    
    # simulate new network 
    g.sim <- simulate(faux.magnolia.high ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE), nsim=1, coef=coefs #true coefficients
                      , control=control.simulate(MCMC.burnin=1000000,MCMC.interval=1000))
    
    #density
    dens[i]<- (sum(g.sim[,])/2)/(1461^2-1461)
    
    # MPLE
    mple <- ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE), estimate = "MPLE")
    coef.matrix.mple[,i] <- mple$coef
    se.matrix.mple[,i] <- summary(mple)$coef[,2]
    
    # MCMCMLE
    
    ## sample size=10000
    
    samplesize=10000 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){
      coef.matrix.10000[,i] <-NA
      se.matrix.10000[,i] <- NA
    }
    else{
      coef.matrix.10000[,i] <-summary(mcmcmle)$coef[,1]
      se.matrix.10000[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    ## sample size=3000
    
    samplesize=3000 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){
      coef.matrix.3000[,i] <-NA
      se.matrix.3000[,i] <- NA
    }
    else{
      coef.matrix.3000[,i] <-summary(mcmcmle)$coef[,1]
      se.matrix.3000[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    ## sample size=1500
    
    samplesize=1500 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){
      coef.matrix.1500[,i] <-NA
      se.matrix.1500[,i] <- NA
    }
    else{
      coef.matrix.1500[,i] <-summary(mcmcmle)$coef[,1]
      se.matrix.1500[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    
    ## sample size=800
    
    samplesize=800 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){ #if model degenerates, then summary(mcmcmle)[2]=="try-error"
      coef.matrix.800[,i] <-NA
      se.matrix.800[,i] <- NA
    }
    else{
      coef.matrix.800[,i] <-summary(mcmcmle)$coef[,1] # otherwise save results
      se.matrix.800[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    
    ## samplesize=400
    
    samplesize=400 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){ #if model degenerates, then summary(mcmcmle)[2]=="try-error"
      coef.matrix.400[,i] <-NA
      se.matrix.400[,i] <- NA
    }
    else{
      coef.matrix.400[,i] <-summary(mcmcmle)$coef[,1] # otherwise save results
      se.matrix.400[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    ## samplesize=200
    
    samplesize=200 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){ #if model degenerates, then summary(mcmcmle)[2]=="try-error"
      coef.matrix.200[,i] <-NA
      se.matrix.200[,i] <- NA
    }
    else{
      coef.matrix.200[,i] <-summary(mcmcmle)$coef[,1] #otherwise save results
      se.matrix.200[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    ## samplesize=100
    
    samplesize=100 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){#if model degenerates, then summary(mcmcmle)[2]=="try-error"
      coef.matrix.100[,i] <-NA
      se.matrix.100[,i] <- NA
    }
    else{
      coef.matrix.100[,i] <-summary(mcmcmle)$coef[,1] # otherwise save results
      se.matrix.100[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    ## sample size=50
    
    samplesize=50 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){#if model degenerates, then summary(mcmcmle)[2]=="try-error"
      coef.matrix.50[,i] <-NA
      se.matrix.50[,i] <- NA
    }
    else{
      coef.matrix.50[,i] <-summary(mcmcmle)$coef[,1] # otherwise save results
      se.matrix.50[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    
    
    
  }
  
  results <- list() #create list to return results
  results[[1]] <- coef.matrix.mple
  results[[2]] <- se.matrix.mple
  
  
  results[[3]] <- coef.matrix.10000 
  results[[4]] <- se.matrix.10000 
  
  results[[5]] <- coef.matrix.3000 
  results[[6]] <- se.matrix.3000 
  
  results[[7]] <- coef.matrix.1500 
  results[[8]] <- se.matrix.1500 
  
  results[[9]] <- coef.matrix.800
  results[[10]] <- se.matrix.800
  
  results[[11]] <- coef.matrix.400
  results[[12]] <- se.matrix.400 
  
  results[[13]] <- coef.matrix.200 
  results[[14]] <- se.matrix.200 
  
  results[[15]] <- coef.matrix.100 
  results[[16]] <- se.matrix.100 
  
  results[[17]] <- coef.matrix.50 
  results[[18]] <- se.matrix.50 

  
  results[[19]] <- dens 
  
  names(results)=c("MPLE Coefficient estimates","MPLE Standard Error",
                   "10000 samplesize coefficient",
                   "10000 samplesize variance","3000 samplesize coefficient",
                   "3000 samplesize variance","1500 samplesize coefficient",
                   "1500 samplesize variance",
                   "800 samplesize coefficient",
                   "800 samplesize variance","400 samplesize coefficient",
                   "400 samplesize variance","200 samplesize coefficient",
                   "200 samplesize variance","100 samplesize coefficient",
                   "100 samplesize variance","50 samplesize coefficient",
                   "50 samplesize variance","Density Vector") 
  
  return(results)
  
}

# run simulation with 500 iterations
faux.magnolia.results <- main.function.faux.magnolia.RMSE(100, coefs= coefficient) 



# RMSE

# coef=Inf appears, thus replace Inf with "NA"
faux.magnolia.results[[1]][!is.finite(faux.magnolia.results[[1]])] <- NA
faux.magnolia.results[[2]][!is.finite(faux.magnolia.results[[2]])] <- NA
faux.magnolia.results[[3]][!is.finite(faux.magnolia.results[[3]])] <- NA
faux.magnolia.results[[4]][!is.finite(faux.magnolia.results[[4]])] <- NA
faux.magnolia.results[[5]][!is.finite(faux.magnolia.results[[5]])] <- NA
faux.magnolia.results[[6]][!is.finite(faux.magnolia.results[[6]])] <- NA
faux.magnolia.results[[7]][!is.finite(faux.magnolia.results[[7]])] <- NA
faux.magnolia.results[[8]][!is.finite(faux.magnolia.results[[8]])] <- NA
faux.magnolia.results[[9]][!is.finite(faux.magnolia.results[[9]])] <- NA
faux.magnolia.results[[10]][!is.finite(faux.magnolia.results[[10]])] <- NA
faux.magnolia.results[[11]][!is.finite(faux.magnolia.results[[11]])] <- NA
faux.magnolia.results[[12]][!is.finite(faux.magnolia.results[[12]])] <- NA
faux.magnolia.results[[13]][!is.finite(faux.magnolia.results[[13]])] <- NA
faux.magnolia.results[[14]][!is.finite(faux.magnolia.results[[14]])] <- NA
faux.magnolia.results[[15]][!is.finite(faux.magnolia.results[[15]])] <- NA
faux.magnolia.results[[16]][!is.finite(faux.magnolia.results[[16]])] <- NA
faux.magnolia.results[[17]][!is.finite(faux.magnolia.results[[17]])] <- NA
faux.magnolia.results[[18]][!is.finite(faux.magnolia.results[[18]])] <- NA

## calculate RMSE

#edges
rmse.mple.edges <-  sqrt(sum((faux.magnolia.results[[1]][1,]-coefficient[1])^2,na.rm = TRUE)) 
rmse.mple.edges

rmse.edges.10000 <- sqrt(sum((faux.magnolia.results[[3]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.10000

rmse.edges.3000 <- sqrt(sum((faux.magnolia.results[[5]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.3000

rmse.edges.1500 <- sqrt(sum((faux.magnolia.results[[7]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.1500

rmse.edges.800 <- sqrt(sum((faux.magnolia.results[[9]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.800

rmse.edges.400 <- sqrt(sum((faux.magnolia.results[[11]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.400

rmse.edges.200 <- sqrt(sum((faux.magnolia.results[[13]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.200

rmse.edges.100 <- sqrt(sum((faux.magnolia.results[[15]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.100

rmse.edges.50 <- sqrt(sum((faux.magnolia.results[[17]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.edges.50




# nodefactor("Sex")

rmse.mple.sex <-  sqrt(sum((faux.magnolia.results[[1]][2,]-coefficient[2])^2,na.rm = TRUE)) 
rmse.mple.sex

rmse.sex.10000 <- sqrt(sum((faux.magnolia.results[[3]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.10000

rmse.sex.3000 <- sqrt(sum((faux.magnolia.results[[5]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.3000

rmse.sex.1500 <- sqrt(sum((faux.magnolia.results[[7]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.1500

rmse.sex.800 <- sqrt(sum((faux.magnolia.results[[9]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.800

rmse.sex.400 <- sqrt(sum((faux.magnolia.results[[11]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.400

rmse.sex.200 <- sqrt(sum((faux.magnolia.results[[13]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.200

rmse.sex.100 <- sqrt(sum((faux.magnolia.results[[15]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.100

rmse.sex.50 <- sqrt(sum((faux.magnolia.results[[17]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.sex.50


# gwesp(0,fixed=T)

rmse.mple.gwesp <-  sqrt(sum((faux.magnolia.results[[1]][3,]-coefficient[3])^2,na.rm = TRUE)) 
rmse.mple.gwesp

rmse.gwesp.10000 <- sqrt(sum((faux.magnolia.results[[3]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.10000

rmse.gwesp.3000 <- sqrt(sum((faux.magnolia.results[[5]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.3000

rmse.gwesp.1500 <- sqrt(sum((faux.magnolia.results[[7]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.1500

rmse.gwesp.800 <- sqrt(sum((faux.magnolia.results[[9]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.800

rmse.gwesp.400 <- sqrt(sum((faux.magnolia.results[[11]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.400

rmse.gwesp.200 <- sqrt(sum((faux.magnolia.results[[13]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.200

rmse.gwesp.100 <- sqrt(sum((faux.magnolia.results[[15]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.100

rmse.gwesp.50 <- sqrt(sum((faux.magnolia.results[[17]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.gwesp.50



# plotting
par(mfrow=c(2,2))
rel.rmse.edges <- c(rmse.edges.50/rmse.mple.edges,rmse.edges.100/rmse.mple.edges, # relative RMSE
                    rmse.edges.200/rmse.mple.edges,rmse.edges.400/rmse.mple.edges,rmse.edges.800/rmse.mple.edges,
                    rmse.edges.1500/rmse.mple.edges,rmse.edges.3000/rmse.mple.edges,rmse.edges.10000/rmse.mple.edges)
log.rel.rmse.edges <- log(rel.rmse.edges) #log relative RMSE
plot(log.rel.rmse.edges, type="l", main="Faux.Magnolia.High, Edges", xaxt = "n",  xlab="MCMC Sample Size", ylab="Log Relative RMSE")
points(log.rel.rmse.edges, pch=20)
axis(1, at=1:8, labels=c(50,100,200,400, 800, 1500, 3000, 10000))
abline(h=0,lty=2, col="grey")

#Sex
rel.rmse.sex <- c(rmse.sex.50/rmse.mple.sex,rmse.sex.100/rmse.mple.sex, #relative RMSE
                  rmse.sex.200/rmse.mple.sex,rmse.sex.400/rmse.mple.sex,rmse.sex.800/rmse.mple.sex,
                  rmse.sex.1500/rmse.mple.sex,rmse.sex.3000/rmse.mple.sex,rmse.sex.10000/rmse.mple.sex)
log.rel.rmse.sex <- log(rel.rmse.sex) #log relative RMSE
plot(log.rel.rmse.sex, type="l",  main="Faux.Magnolia.High, Nodefactor(Sex)", xaxt = "n", xlab="MCMC Sample Size", ylab="Log Relative RMSE")
points(log.rel.rmse.sex, pch=20)
axis(1, at=1:8, labels=c(50,100,200,400, 800, 1500, 3000, 10000))
abline(h=0,lty=2, col="grey")

#gwesp
rel.rmse.gwesp <- c(rmse.gwesp.50/rmse.mple.gwesp,rmse.gwesp.100/rmse.mple.gwesp, #relative RMSE
                    rmse.gwesp.200/rmse.mple.gwesp, rmse.gwesp.400/rmse.mple.gwesp,rmse.gwesp.800/rmse.mple.gwesp,
                    rmse.gwesp.1500/rmse.mple.gwesp,rmse.gwesp.3000/rmse.mple.gwesp,rmse.gwesp.10000/rmse.mple.gwesp)
log.rel.rmse.gwesp <- log(rel.rmse.gwesp) #log relative RMSE
plot(log.rel.rmse.gwesp, type="l", main="Faux.Magnolia.High, Gwesp", xaxt = "n", xlab="MCMC Sample Size", ylab="Log Relative RMSE")
points(log.rel.rmse.gwesp, pch=20)
axis(1, at=1:8, labels=c(50,100,200, 400, 800, 1500, 3000, 10000))
abline(h=0,lty=2, col="grey")


# Histogram of the density of simulated networks
hist(faux.magnolia.results[[19]], main="Density of Simulated Networks", xlab="Density")
dens.faux.magnolia=974/(1461^2-1461)
abline(v=dens.faux.mesa,col="red")
