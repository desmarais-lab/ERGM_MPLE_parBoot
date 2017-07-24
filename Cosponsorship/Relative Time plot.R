# comparison plot MCMLE vs MPLE

# time in seconds

library(statnet)

#############################
# mesa high
data("faux.mesa.high")

time.mple.mesa<- c()
set.seed(1)
for(i in 1:100){
  b<- system.time(ergm(faux.mesa.high~edges+nodefactor("Sex")+gwesp(0.25,fixed=TRUE), estimate="MPLE"))
  time.mple.mesa[i]<- b[3]                  
}

median(time.mple.mesa)


# MCMLE runtime

time.mcmle.mesa<- c() # for storing the time

set.seed(1)
for(i in 1:20){
  print(i)
  b<- system.time(ergm(faux.mesa.high~edges+nodefactor("Sex")+gwesp(0.25,fixed=TRUE), 
                       eval.loglik = FALSE, control = control.ergm(MCMC.interval = 2000,MCMC.samplesize=2000)))
  time.mcmle.mesa[i]<- b[3]                  
}

median(time.mcmle.mesa)

time.mple.mesa.med<- median(time.mple.mesa)
time.mcmle.mesa.med<- median(time.mcmle.mesa)
time.sim.mesa.med<- median(time.mcmle.mesa)/5 #takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)

boot.mple.mesa<-c()
boot.mple.mesa[1]<- 500*time.mple.mesa.med/1 + time.sim.mesa.med
boot.mple.mesa[2]<- 500*time.mple.mesa.med/2 + time.sim.mesa.med
boot.mple.mesa[3]<- 500*time.mple.mesa.med/3 + time.sim.mesa.med
boot.mple.mesa[4]<- 500*time.mple.mesa.med/4 + time.sim.mesa.med
boot.mple.mesa[5]<- 500*time.mple.mesa.med/5 + time.sim.mesa.med
boot.mple.mesa[6]<- 500*time.mple.mesa.med/6 + time.sim.mesa.med
boot.mple.mesa[7]<- 500*time.mple.mesa.med/7 + time.sim.mesa.med
boot.mple.mesa[8]<- 500*time.mple.mesa.med/8 + time.sim.mesa.med

#relative time

rel.time.mesa<- boot.mple.mesa/time.mcmle.mesa.med

# where does the ratio converge to?
more.cores <- 1:100
boot.mple.mesa.long<- 500*time.mple.mesa.med/more.cores + time.sim.mesa.med
rel.mesa.long<- boot.mple.mesa.long/time.mcmle.mesa.med
rel.mesa.long # 0.218

########################################
# magnolia high
data("faux.magnolia.high")

time.mple.magno<- c()
set.seed(1)
for(i in 1:100){
  b<- system.time(ergm(faux.magnolia.high~edges+nodefactor("Sex")+gwesp(0.25,fixed=TRUE), estimate="MPLE"))
  time.mple.magno[i]<- b[3]                  
}

median(time.mple.magno)

# MCMLE runtime

time.mcmle.magno<- c() # for storing the time

set.seed(1)
for(i in 1:5){
  print(i)
  b<- system.time(ergm(faux.magnolia.high~edges+nodefactor("Sex")+gwesp(0.25,fixed=TRUE), 
                       eval.loglik = FALSE,control=control.ergm(MCMC.interval = 5000,MCMC.samplesize=8000)))
  time.mcmle.magno[i]<- b[3]                  
}

median(time.mcmle.magno)



time.mple.magno.med<- median(time.mple.magno)
time.mcmle.magno.med <- median(time.mcmle.magno)
time.sim.magno.med<- median(time.mcmle.magno)/6 #takes 3 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)

boot.mple.magnolia<-c()
boot.mple.magnolia[1]<- 500*time.mple.magno.med/1 + time.sim.magno.med
boot.mple.magnolia[2]<- 500*time.mple.magno.med/2 + time.sim.magno.med
boot.mple.magnolia[3]<- 500*time.mple.magno.med/3 + time.sim.magno.med
boot.mple.magnolia[4]<- 500*time.mple.magno.med/4 + time.sim.magno.med
boot.mple.magnolia[5]<- 500*time.mple.magno.med/5 + time.sim.magno.med
boot.mple.magnolia[6]<- 500*time.mple.magno.med/6 + time.sim.magno.med
boot.mple.magnolia[7]<- 500*time.mple.magno.med/7 + time.sim.magno.med
boot.mple.magnolia[8]<- 500*time.mple.magno.med/8 + time.sim.magno.med

#relative time

rel.time.magnolia<- boot.mple.magnolia/time.mcmle.magno.med

# where does the ratio converge to?
more.cores <- 1:100
boot.mple.magnolia.long<- 500*time.mple.magno.med/more.cores + time.sim.magno.med
rel.magnolia.long<- boot.mple.magnolia.long/time.mcmle.magno.med
rel.magnolia.long # 0.188



###################################################
# relative time cosponsorship

load(file="cosponsorship.RData")

# MPLE runtime

time.mple.co<- c()
set.seed(1)
for(i in 1:100){
  b<- system.time(ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), estimate="MPLE"))
  time.mple.co[i]<- b[3]                  
}

median(time.mple.co)


# MCMLE runtime

time.mcmle.co<- c() # for storing the time

set.seed(1)
for(i in 1:5){
  print(i)
  b<- system.time(ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), 
                       eval.loglik = FALSE, control = control.ergm(MCMC.interval = 30000, MCMC.samplesize=10000)))
  time.mcmle.co[i]<- b[3]                  
}

median(time.mcmle.co)

time.mple.co.med<- median(time.mple.co)
time.mcmle.co.med <- median(time.mcmle.co)
time.sim.co.med<- median(time.mcmle.co)/5 #takes 2 iterations, MCMC samplesize is quadrupled in last iteration

boot.mple.co<-c()
boot.mple.co[1]<- 500*time.mple.co.med/1 + time.sim.co.med
boot.mple.co[2]<- 500*time.mple.co.med/2 + time.sim.co.med
boot.mple.co[3]<- 500*time.mple.co.med/3 + time.sim.co.med
boot.mple.co[4]<- 500*time.mple.co.med/4 + time.sim.co.med
boot.mple.co[5]<- 500*time.mple.co.med/5 + time.sim.co.med
boot.mple.co[6]<- 500*time.mple.co.med/6 + time.sim.co.med
boot.mple.co[7]<- 500*time.mple.co.med/7 + time.sim.co.med
boot.mple.co[8]<- 500*time.mple.co.med/8 + time.sim.co.med

#relative time

rel.time.co<- boot.mple.co/time.mcmle.co.med

# where does the ratio converge to?
more.cores <- 1:100
boot.mple.co.long<- 500*time.mple.co.med/more.cores + time.sim.co.med
rel.co.long<- boot.mple.co.long/time.mcmle.co.med
rel.co.long #0.215



################################
## Plots
#################################



# plot
plot(1:8, rel.mesa.long[1:8], main="Relative Time", xlab="Number of Cores", ylab="Relative Time", type="l", xaxt="n", 
     col="blue", lwd=2, cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=c(0.2,2.3))
axis(1, 1:8, label=1:8, cex.axis=1.5)
lines(1:8, rel.magnolia.long[1:8], col="red", lwd=2)
lines(1:8, rel.co.long[1:8], col="darkgreen", lwd=2)
abline(h=1,col="grey", lty=2, lwd=2)

legend("topright", inset=c(.05, .01), legend=c("Faux Mesa High","Faux Magnolia High", "Cosponsorship"),lty = 1,lwd=2, col = c("blue", "red", "darkgreen"), cex=1.5, bty="n")
# 7x10 inches


##################################
# plot for presentation

# plot
plot(1:8, rel.time.mesa, main="Relative Time", xlab="Number of Cores", ylab="Relative Time", type="l", xaxt="n", 
     col="blue", lwd=2, cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=c(0.2,2.3))
axis(1, 1:8, label=1:8, cex.axis=1.5)
lines(1:8, rel.time.magnolia, col="red", lwd=2)
abline(h=1,col="grey", lty=2, lwd=2)

legend("topright", inset=c(.05, .01), legend=c("Mesa High","Magnolia High"),lty = 1,lwd=2, col = c("blue", "red"), cex=1.5, bty="n")
# 7x10 inches