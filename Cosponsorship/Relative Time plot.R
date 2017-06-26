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
                       eval.loglik = FALSE, control = control.ergm(MCMC.interval = 2000,MCMC.samplesize=1000)))
  time.mcmle.mesa[i]<- b[3]                  
}

median(time.mcmle.mesa)


time.mple<- median(time.mple.mesa)
time.mcmle<- median(time.mcmle.mesa)
time.sim<- median(time.mcmle.mesa)/5 #takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)

boot.mple.mesa<-c()
boot.mple.mesa[1]<- 500*time.mple/1 + time.sim
boot.mple.mesa[2]<- 500*time.mple/5 + time.sim
boot.mple.mesa[3]<- 500*time.mple/10 + time.sim
boot.mple.mesa[4]<- 500*time.mple/25 + time.sim
boot.mple.mesa[5]<- 500*time.mple/50 + time.sim
boot.mple.mesa[6]<- 500*time.mple/100 + time.sim
boot.mple.mesa[7]<- 500*time.mple/250 + time.sim
boot.mple.mesa[8]<- 500*time.mple/500 + time.sim

#relative time

rel.time.mesa<- boot.mple.mesa/time.mcmle

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
for(i in 1:4){
  print(i)
  b<- system.time(ergm(faux.magnolia.high~edges+nodefactor("Sex")+gwesp(0.25,fixed=TRUE), 
                       eval.loglik = FALSE,control=control.ergm(MCMC.interval = 5000,MCMC.samplesize=8000)))
  time.mcmle.magno[i]<- b[3]                  
}

median(time.mcmle.magno)



time.mple<- median(time.mple.magno)
time.mcmle <- median(time.mcmle.magno)
time.sim<- median(time.mcmle.magno)/5 #takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)
boot.mple.magnolia<-c()
boot.mple.magnolia[1]<- 500*time.mple/1+ time.sim
boot.mple.magnolia[2]<- 500*time.mple/5+ time.sim
boot.mple.magnolia[3]<- 500*time.mple/10+ time.sim
boot.mple.magnolia[4]<- 500*time.mple/25+ time.sim
boot.mple.magnolia[5]<- 500*time.mple/50+ time.sim
boot.mple.magnolia[6]<- 500*time.mple/100+ time.sim
boot.mple.magnolia[7]<- 500*time.mple/250+ time.sim
boot.mple.magnolia[8]<- 500*time.mple/500+ time.sim

#relative time

rel.time.magnolia<- boot.mple.magnolia/time.mcmle




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
for(i in 1:10){
  print(i)
  b<- system.time(ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), 
                       eval.loglik = FALSE, control = control.ergm(MCMC.interval = 30000, MCMC.burnin = 300000,MCMC.samplesize=20000)))
  time.mcmle.co[i]<- b[3]                  
}

median(time.mcmle.co)

time.mple<- median(time.mple.co)
time.mcmle<-median(time.mcmle.co)
time.sim <- median(time.mcmle.co)/5 # takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)
boot.mple.co<-c()
boot.mple.co[1]<- 500*time.mple/1 + time.sim
boot.mple.co[2]<- 500*time.mple/5+ time.sim
boot.mple.co[3]<- 500*time.mple/10+ time.sim
boot.mple.co[4]<- 500*time.mple/25+ time.sim
boot.mple.co[5]<- 500*time.mple/50+ time.sim
boot.mple.co[6]<- 500*time.mple/100+ time.sim
boot.mple.co[7]<- 500*time.mple/250+ time.sim
boot.mple.co[8]<- 500*time.mple/500+ time.sim


#relative time

rel.time.co<- boot.mple.co/time.mcmle


# plot
par(mfrow=c(1,1))

plot(1:8, rel.time.mesa, main="Relative Time", xlab="Number of Cores", ylab="Relative Time", type="l", xaxt="n", 
     col="blue", lwd=2, cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=c(0,2.3))
axis(1, 1:8, label=c(1,5,10,25, 50,100,250,500))
lines(1:8, rel.time.magnolia, col="red", lwd=2)
lines(1:8, rel.time.co, col="darkgreen", lwd=2)
abline(h=1,col="grey", lty=2, lwd=2)

legend("topright", inset=c(.05, .01), legend=c("Faux Mesa High","Faux Magnolia High", "Cosponsorship"),
       lty = 1,lwd=2, col = c("blue", "red", "darkgreen"), cex=1.5, bty="n")




################################
## relative time less cores
#################################

# comparison plot MCMLE vs MPLE

# mesa high

time.mple<- median(time.mple.mesa)
time.mcmle<- median(time.mcmle.mesa)
time.sim<- median(time.mcmle.mesa)/5 #takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)

boot.mple.mesa<-c()
boot.mple.mesa[1]<- 500*time.mple/1 + time.sim
boot.mple.mesa[2]<- 500*time.mple/2 + time.sim
boot.mple.mesa[3]<- 500*time.mple/3 + time.sim
boot.mple.mesa[4]<- 500*time.mple/4 + time.sim
boot.mple.mesa[5]<- 500*time.mple/5 + time.sim
boot.mple.mesa[6]<- 500*time.mple/6 + time.sim
boot.mple.mesa[7]<- 500*time.mple/7 + time.sim
boot.mple.mesa[8]<- 500*time.mple/8 + time.sim

#relative time

rel.time.mesa<- boot.mple.mesa/time.mcmle

# magnolia high

time.mple<- median(time.mple.magno)
time.mcmle <- median(time.mcmle.magno)
time.sim<- median(time.mcmle.magno)/5 #takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)
boot.mple.magnolia<-c()
boot.mple.magnolia[1]<- 500*time.mple/1+ time.sim
boot.mple.magnolia[2]<- 500*time.mple/2+ time.sim
boot.mple.magnolia[3]<- 500*time.mple/3+ time.sim
boot.mple.magnolia[4]<- 500*time.mple/4+ time.sim
boot.mple.magnolia[5]<- 500*time.mple/5+ time.sim
boot.mple.magnolia[6]<- 500*time.mple/6+ time.sim
boot.mple.magnolia[7]<- 500*time.mple/7+ time.sim
boot.mple.magnolia[8]<- 500*time.mple/8+ time.sim

#relative time

rel.time.magnolia<- boot.mple.magnolia/time.mcmle


# relative time cosponsorship

time.mple<- median(time.mple.co)
time.mcmle<-median(time.mcmle.co)
time.sim <- median(time.mcmle.co)/5 # takes 2 iterations, MCMC samplesize is quadrupled in last iteration

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)
boot.mple.co<-c()
boot.mple.co[1]<- 500*time.mple/1 + time.sim
boot.mple.co[2]<- 500*time.mple/2+ time.sim
boot.mple.co[3]<- 500*time.mple/3+ time.sim
boot.mple.co[4]<- 500*time.mple/4+ time.sim
boot.mple.co[5]<- 500*time.mple/5+ time.sim
boot.mple.co[6]<- 500*time.mple/6+ time.sim
boot.mple.co[7]<- 500*time.mple/7+ time.sim
boot.mple.co[8]<- 500*time.mple/8+ time.sim

#relative time

rel.time.sc<- boot.mple.co/time.mcmle


# plot
plot(1:8, rel.time.mesa, main="Relative Time", xlab="Number of Cores", ylab="Relative Time", type="l", xaxt="n", 
     col="blue", lwd=2, cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=c(0.2,2.3))
axis(1, 1:8, label=1:8, cex.axis=1.5)
lines(1:8, rel.time.magnolia, col="red", lwd=2)
lines(1:8, rel.time.co, col="darkgreen", lwd=2)
abline(h=1,col="grey", lty=2, lwd=2)

legend("topright", inset=c(.05, .01), legend=c("Faux Mesa High","Faux Magnolia High", "Cosponsorship"),lty = 1,lwd=2, col = c("blue", "red", "darkgreen"), cex=1.5, bty="n")
# 7x10 inches