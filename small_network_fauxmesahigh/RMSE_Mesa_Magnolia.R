# calculating RMSE for faux mesa high and faux magnolia high

# first Faux mesa high

library(statnet)

data("faux.mesa.high")

plot(faux.mesa.high)

set.seed(123)
model=ergm(faux.mesa.high~edges+nodefactor("Sex")+gwesp(0,fixed=TRUE))
summary(model)
#get true coefficients
coefficient <- summary(model)$coef[,1]

# function to calculate RMSE for MPLE, MCMCMLE with samplesize= 1500, 800, 400, 200, 100, 50, 20
main.function.faux.mesa.RMSE <- function(iters,  coefs=coefficient){
  # size = number of nodes in network
  
  set.seed(999)
  
  l=length(coefs) # number of coefficients
  
  #creating empty matrices to store results
  coef.matrix.mple <- matrix(0,l,iters)
  se.matrix.mple <- matrix(0,l,iters)
  
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
  
  coef.matrix.25 <- matrix(0,l,iters)
  se.matrix.25 <- matrix(0,l,iters)
  
  
  # start loop
  for(i in 1:iters){
    print(i)
    
    # simulate new network 
    g.sim <- simulate(faux.mesa.high ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE), nsim=1, coef=coefs #true coefficients
                      , control=control.simulate(MCMC.burnin=1000000,MCMC.interval=1000))
    
    # MPLE
    mple <- ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE), estimate = "MPLE")
    coef.matrix.mple[,i] <- mple$coef
    se.matrix.mple[,i] <- summary(mple)$coef[,2]
    
    # MCMCMLE
    
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
    
    
    ## sample size=25
    
    samplesize=25 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){
      coef.matrix.25[,i] <-NA
      se.matrix.25[,i] <- NA
    }
    else{
      coef.matrix.25[,i] <-summary(mcmcmle)$coef[,1]
      se.matrix.25[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    
  }
  
  results <- list() #create list to return results
  results[[1]] <- coef.matrix.mple
  results[[2]] <- se.matrix.mple
  
  results[[3]] <- coef.matrix.1500 
  results[[4]] <- se.matrix.1500 
  
  results[[5]] <- coef.matrix.800
  results[[6]] <- se.matrix.800
  
  results[[7]] <- coef.matrix.400
  results[[8]] <- se.matrix.400 
  
  results[[9]] <- coef.matrix.200 
  results[[10]] <- se.matrix.200 
  
  results[[11]] <- coef.matrix.100 
  results[[12]] <- se.matrix.100 
  
  results[[13]] <- coef.matrix.50 
  results[[14]] <- se.matrix.50 
  
  results[[15]] <- coef.matrix.25 
  results[[16]] <- se.matrix.25 
  
  
  names(results)=c("MPLE Coefficient estimates","MPLE Standard Error",
                   "1500 samplesize coefficient", "1500 samplesize variance",
                   "800 samplesize coefficient",
                   "800 samplesize variance","400 samplesize coefficient",
                   "400 samplesize variance","200 samplesize coefficient",
                   "200 samplesize variance","100 samplesize coefficient",
                   "100 samplesize variance","50 samplesize coefficient",
                   "50 samplesize variance","25 samplesize coefficient",
                   "25 samplesize variance") 
  
  return(results)
  
}

# run simulation with 500 iterations
faux.mesa.results <- main.function.faux.mesa.RMSE(500, coefs= coefficient) 


# RMSE

# coef=Inf appears, thus replace Inf with "NA"
faux.mesa.results[[1]][!is.finite(faux.mesa.results[[1]])] <- NA
faux.mesa.results[[2]][!is.finite(faux.mesa.results[[2]])] <- NA
faux.mesa.results[[3]][!is.finite(faux.mesa.results[[3]])] <- NA
faux.mesa.results[[4]][!is.finite(faux.mesa.results[[4]])] <- NA
faux.mesa.results[[5]][!is.finite(faux.mesa.results[[5]])] <- NA
faux.mesa.results[[6]][!is.finite(faux.mesa.results[[6]])] <- NA
faux.mesa.results[[7]][!is.finite(faux.mesa.results[[7]])] <- NA
faux.mesa.results[[8]][!is.finite(faux.mesa.results[[8]])] <- NA
faux.mesa.results[[9]][!is.finite(faux.mesa.results[[9]])] <- NA
faux.mesa.results[[10]][!is.finite(faux.mesa.results[[10]])] <- NA
faux.mesa.results[[11]][!is.finite(faux.mesa.results[[11]])] <- NA
faux.mesa.results[[12]][!is.finite(faux.mesa.results[[12]])] <- NA
faux.mesa.results[[13]][!is.finite(faux.mesa.results[[13]])] <- NA
faux.mesa.results[[14]][!is.finite(faux.mesa.results[[14]])] <- NA
faux.mesa.results[[15]][!is.finite(faux.mesa.results[[15]])] <- NA
faux.mesa.results[[16]][!is.finite(faux.mesa.results[[16]])] <- NA


## calculate RMSE

#edges
rmse.mesa.mple.edges <-  sqrt(sum((faux.mesa.results[[1]][1,]-coefficient[1])^2,na.rm = TRUE)) 
rmse.mesa.mple.edges

rmse.mesa.edges.1500 <- sqrt(sum((faux.mesa.results[[3]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.1500

rmse.mesa.edges.800 <- sqrt(sum((faux.mesa.results[[5]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.800

rmse.mesa.edges.400 <- sqrt(sum((faux.mesa.results[[7]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.400

rmse.mesa.edges.200 <- sqrt(sum((faux.mesa.results[[9]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.200

rmse.mesa.edges.100 <- sqrt(sum((faux.mesa.results[[11]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.100

rmse.mesa.edges.50 <- sqrt(sum((faux.mesa.results[[13]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.50

rmse.mesa.edges.25 <- sqrt(sum((faux.mesa.results[[15]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.mesa.edges.25



# nodefactor("Sex")

rmse.mesa.mple.sex <-  sqrt(sum((faux.mesa.results[[1]][2,]-coefficient[2])^2,na.rm = TRUE)) 
rmse.mesa.mple.sex

rmse.mesa.sex.1500 <- sqrt(sum((faux.mesa.results[[3]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.1500

rmse.mesa.sex.800 <- sqrt(sum((faux.mesa.results[[5]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.800

rmse.mesa.sex.400 <- sqrt(sum((faux.mesa.results[[7]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.400

rmse.mesa.sex.200 <- sqrt(sum((faux.mesa.results[[9]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.200

rmse.mesa.sex.100 <- sqrt(sum((faux.mesa.results[[11]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.100

rmse.mesa.sex.50 <- sqrt(sum((faux.mesa.results[[13]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.50

rmse.mesa.sex.25 <- sqrt(sum((faux.mesa.results[[15]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.mesa.sex.25


# gwesp(0,fixed=T)

rmse.mesa.mple.gwesp <-  sqrt(sum((faux.mesa.results[[1]][3,]-coefficient[3])^2,na.rm = TRUE)) 
rmse.mesa.mple.gwesp

rmse.mesa.gwesp.1500 <- sqrt(sum((faux.mesa.results[[3]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.1500

rmse.mesa.gwesp.800 <- sqrt(sum((faux.mesa.results[[5]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.800

rmse.mesa.gwesp.400 <- sqrt(sum((faux.mesa.results[[7]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.400

rmse.mesa.gwesp.200 <- sqrt(sum((faux.mesa.results[[9]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.200

rmse.mesa.gwesp.100 <- sqrt(sum((faux.mesa.results[[11]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.100

rmse.mesa.gwesp.50 <- sqrt(sum((faux.mesa.results[[13]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.50

rmse.mesa.gwesp.25 <- sqrt(sum((faux.mesa.results[[15]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.mesa.gwesp.25



######################
## faux magnolia high
######################




data("faux.magnolia.high")

plot(faux.magnolia.high)

set.seed(123)
model=ergm(faux.magnolia.high~edges+nodefactor("Sex")+gwesp(0,fixed=TRUE))
summary(model)
#get true coefficients
coefficient <- summary(model)$coef[,1]

# function to calculate RMSE for MPLE, MCMCMLE with samplesize= 1500, 800, 400, 200, 100, 50, 25
main.function.faux.magnolia.RMSE <- function(iters,  coefs=coefficient){
  # size = number of nodes in network
  
  set.seed(999)
  
  l=length(coefs)
  
  #creating empty matrices to store results
  coef.matrix.mple <- matrix(0,l,iters)
  se.matrix.mple <- matrix(0,l,iters)
  
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
  
  coef.matrix.25 <- matrix(0,l,iters)
  se.matrix.25 <- matrix(0,l,iters)
  
  
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
    
    ## sample size=25
    
    samplesize=25 # MCMC samplesize
    
    mcmcmle <- try(ergm(g.sim ~ edges+nodefactor("Sex")+gwesp(0,fixed=TRUE),control=control.ergm(MCMC.samplesize=samplesize)))
    if(summary(mcmcmle)[2]=="try-error"){#if model degenerates, then summary(mcmcmle)[2]=="try-error"
      coef.matrix.25[,i] <-NA
      se.matrix.25[,i] <- NA
    }
    else{
      coef.matrix.25[,i] <-summary(mcmcmle)$coef[,1] # otherwise save results
      se.matrix.25[,i] <- summary(mcmcmle)$coef[,2]
    }
    
    
  }
  
  results <- list() #create list to return results
  results[[1]] <- coef.matrix.mple
  results[[2]] <- se.matrix.mple

  results[[3]] <- coef.matrix.1500 
  results[[4]] <- se.matrix.1500 
  
  results[[5]] <- coef.matrix.800
  results[[6]] <- se.matrix.800
  
  results[[7]] <- coef.matrix.400
  results[[8]] <- se.matrix.400 
  
  results[[9]] <- coef.matrix.200 
  results[[10]] <- se.matrix.200 
  
  results[[11]] <- coef.matrix.100 
  results[[12]] <- se.matrix.100 
  
  results[[13]] <- coef.matrix.50 
  results[[14]] <- se.matrix.50 
  
  results[[15]] <- coef.matrix.25 
  results[[16]] <- se.matrix.25 
  
  
  results[[17]] <- dens 
  
  names(results)=c("MPLE Coefficient estimates","MPLE Standard Error",
                   "1500 samplesize coefficient",
                   "1500 samplesize variance",
                   "800 samplesize coefficient",
                   "800 samplesize variance","400 samplesize coefficient",
                   "400 samplesize variance","200 samplesize coefficient",
                   "200 samplesize variance","100 samplesize coefficient",
                   "100 samplesize variance","50 samplesize coefficient",
                   "50 samplesize variance","25 samplesize coefficient",
                   "25 samplesize variance","Density Vector") 
  
  return(results)
  
}

# run simulation with 500 iterations
faux.magnolia.results <- main.function.faux.magnolia.RMSE(500, coefs= coefficient) 


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

## calculate RMSE

#edges
rmse.magno.mple.edges <-  sqrt(sum((faux.magnolia.results[[1]][1,]-coefficient[1])^2,na.rm = TRUE)) 
rmse.magno.mple.edges

rmse.magno.edges.1500 <- sqrt(sum((faux.magnolia.results[[3]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.1500

rmse.magno.edges.800 <- sqrt(sum((faux.magnolia.results[[5]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.800

rmse.magno.edges.400 <- sqrt(sum((faux.magnolia.results[[7]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.400

rmse.magno.edges.200 <- sqrt(sum((faux.magnolia.results[[9]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.200

rmse.magno.edges.100 <- sqrt(sum((faux.magnolia.results[[11]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.100

rmse.magno.edges.50 <- sqrt(sum((faux.magnolia.results[[13]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.50

rmse.magno.edges.25 <- sqrt(sum((faux.magnolia.results[[15]][1,]-coefficient[1])^2,na.rm = TRUE))
rmse.magno.edges.25


# nodefactor("Sex")

rmse.magno.mple.sex <-  sqrt(sum((faux.magnolia.results[[1]][2,]-coefficient[2])^2,na.rm = TRUE)) 
rmse.magno.mple.sex

rmse.magno.sex.1500 <- sqrt(sum((faux.magnolia.results[[3]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.1500

rmse.magno.sex.800 <- sqrt(sum((faux.magnolia.results[[5]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.800

rmse.magno.sex.400 <- sqrt(sum((faux.magnolia.results[[7]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.400

rmse.magno.sex.200 <- sqrt(sum((faux.magnolia.results[[9]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.200

rmse.magno.sex.100 <- sqrt(sum((faux.magnolia.results[[11]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.100

rmse.magno.sex.50 <- sqrt(sum((faux.magnolia.results[[13]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.50

rmse.magno.sex.25 <- sqrt(sum((faux.magnolia.results[[15]][2,]-coefficient[2])^2,na.rm = TRUE))
rmse.magno.sex.25


# gwesp(0,fixed=T)

rmse.magno.mple.gwesp <-  sqrt(sum((faux.magnolia.results[[1]][3,]-coefficient[3])^2,na.rm = TRUE)) 
rmse.magno.mple.gwesp

rmse.magno.gwesp.1500 <- sqrt(sum((faux.magnolia.results[[3]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.1500

rmse.magno.gwesp.800 <- sqrt(sum((faux.magnolia.results[[5]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.800

rmse.magno.gwesp.400 <- sqrt(sum((faux.magnolia.results[[7]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.400

rmse.magno.gwesp.200 <- sqrt(sum((faux.magnolia.results[[9]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.200

rmse.magno.gwesp.100 <- sqrt(sum((faux.magnolia.results[[11]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.100

rmse.magno.gwesp.50 <- sqrt(sum((faux.magnolia.results[[13]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.50

rmse.magno.gwesp.25 <- sqrt(sum((faux.magnolia.results[[15]][3,]-coefficient[3])^2,na.rm = TRUE))
rmse.magno.gwesp.25



#### Plotting

par(xpd=FALSE)
par(mfrow=c(2,2))
par(mar=c(5.1, 4.1, 4.1, 2.1))
rel.rmse.mesa.edges <- c(rmse.mesa.edges.25/rmse.mesa.mple.edges,rmse.mesa.edges.50/rmse.mesa.mple.edges,rmse.mesa.edges.100/rmse.mesa.mple.edges, # relative RMSE
                    rmse.mesa.edges.200/rmse.mesa.mple.edges,rmse.mesa.edges.400/rmse.mesa.mple.edges,rmse.mesa.edges.800/rmse.mesa.mple.edges,
                    rmse.mesa.edges.1500/rmse.mesa.mple.edges)
rel.rmse.magnolia.edges <- c(rmse.magno.edges.25/rmse.magno.mple.edges,rmse.magno.edges.50/rmse.magno.mple.edges,rmse.magno.edges.100/rmse.magno.mple.edges, # relative RMSE
                    rmse.magno.edges.200/rmse.magno.mple.edges,rmse.magno.edges.400/rmse.magno.mple.edges,rmse.magno.edges.800/rmse.magno.mple.edges,
                    rmse.magno.edges.1500/rmse.magno.mple.edges)
log.rel.mesa.rmse.edges <- log(rel.rmse.mesa.edges) #log relative RMSE
log.rel.magnolia.rmse.edges <- log(rel.rmse.magnolia.edges) #log relative RMSE
plot(log.rel.mesa.rmse.edges, type="l", main="Edges", xaxt = "n", xlab="MCMC Sample Size", cex.main=2,cex.axis=1.3, cex.lab=1.5,ylab="Log Relative RMSE", ylim=c(-0.15, 0.65))
lines(log.rel.magnolia.rmse.edges, type="l", lty=2)
points(log.rel.mesa.rmse.edges, pch=20)
points(log.rel.magnolia.rmse.edges, pch=20)
axis(1, at=1:7, labels=c(25,50,100,200,400, 800, 1500), cex.axis=1.3)
abline(h=0,lty=2, col="grey")

#Sex
rel.rmse.mesa.sex <- c(rmse.mesa.sex.25/rmse.mesa.mple.sex,rmse.mesa.sex.50/rmse.mesa.mple.sex,rmse.mesa.sex.100/rmse.mesa.mple.sex, # relative RMSE
                         rmse.mesa.sex.200/rmse.mesa.mple.sex,rmse.mesa.sex.400/rmse.mesa.mple.sex,rmse.mesa.sex.800/rmse.mesa.mple.sex,
                         rmse.mesa.sex.1500/rmse.mesa.mple.sex)
rel.rmse.magnolia.sex <- c(rmse.magno.sex.25/rmse.magno.mple.sex,rmse.magno.sex.50/rmse.magno.mple.sex,rmse.magno.sex.100/rmse.magno.mple.sex, # relative RMSE
                             rmse.magno.sex.200/rmse.magno.mple.sex,rmse.magno.sex.400/rmse.magno.mple.sex,rmse.magno.sex.800/rmse.magno.mple.sex,
                             rmse.magno.sex.1500/rmse.magno.mple.sex)
log.rel.mesa.rmse.sex <- log(rel.rmse.mesa.sex) #log relative RMSE
log.rel.magnolia.rmse.sex <- log(rel.rmse.magnolia.sex) #log relative RMSE
plot(log.rel.mesa.rmse.sex, type="l", main="Sex", xaxt = "n", xlab="MCMC Sample Size",cex.main=2,cex.axis=1.3, cex.lab=1.5, ylab="Log Relative RMSE", ylim=c(-0.25, 0.75))
lines(log.rel.magnolia.rmse.sex, type="l", lty=2)
points(log.rel.mesa.rmse.sex, pch=20)
points(log.rel.magnolia.rmse.sex, pch=20)
axis(1, at=1:7, labels=c(25,50,100,200,400, 800, 1500), cex.axis=1.3)
abline(h=0,lty=2, col="grey")

#gwesp
rel.rmse.mesa.gwesp <- c(rmse.mesa.gwesp.25/rmse.mesa.mple.gwesp,rmse.mesa.gwesp.50/rmse.mesa.mple.gwesp,rmse.mesa.gwesp.100/rmse.mesa.mple.gwesp, # relative RMSE
                       rmse.mesa.gwesp.200/rmse.mesa.mple.gwesp,rmse.mesa.gwesp.400/rmse.mesa.mple.gwesp,rmse.mesa.gwesp.800/rmse.mesa.mple.gwesp,
                       rmse.mesa.gwesp.1500/rmse.mesa.mple.gwesp)
rel.rmse.magnolia.gwesp <- c(rmse.magno.gwesp.25/rmse.magno.mple.gwesp,rmse.magno.gwesp.50/rmse.magno.mple.gwesp,rmse.magno.gwesp.100/rmse.magno.mple.gwesp, # relative RMSE
                           rmse.magno.gwesp.200/rmse.magno.mple.gwesp,rmse.magno.gwesp.400/rmse.magno.mple.gwesp,rmse.magno.gwesp.800/rmse.magno.mple.gwesp,
                           rmse.magno.gwesp.1500/rmse.magno.mple.gwesp)
log.rel.mesa.rmse.gwesp <- log(rel.rmse.mesa.gwesp) #log relative RMSE
log.rel.magnolia.rmse.gwesp <- log(rel.rmse.magnolia.gwesp) #log relative RMSE
plot(log.rel.mesa.rmse.gwesp, type="l", main="Gwesp", xaxt = "n", xlab="MCMC Sample Size",cex.main=2,cex.axis=1.3, cex.lab=1.5, ylab="Log Relative RMSE", ylim=c(-0.1, 2.7))
lines(log.rel.magnolia.rmse.gwesp, type="l", lty=2)
points(log.rel.mesa.rmse.gwesp, pch=20)
points(log.rel.magnolia.rmse.gwesp, pch=20)
axis(1, at=1:7, labels=c(25,50,100,200,400, 800, 1500), cex.axis=1.3)
abline(h=0,lty=2, col="grey")

par(mar=c(0, 0, 0, 0))
# c(bottom, left, top, right)
plot.new()
legend("center", legend=c("Mesa High","Magnolia High"), lty = c(1,2),lwd=2,
       bty ="n", cex=2)

