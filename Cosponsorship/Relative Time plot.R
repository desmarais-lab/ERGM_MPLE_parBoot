# comparison plot MCMLE vs MPLE

# time in seconds

# mesa high

time.mple<- 0.05
time.mcmle<- 12.59
time.sim<- 1.573 # simulation time

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

# magnolia high

time.mple<- 0.51
time.mcmle <- 119.07
time.sim<- 2.49 # simulation time

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


# relative time scnetwork

time.mple<- 3218
time.mcmle<- 1123200
time.sim <- 10800

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)
boot.mple.sc<-c()
boot.mple.sc[1]<- 500*time.mple/1 + time.sim
boot.mple.sc[2]<- 500*time.mple/5+ time.sim
boot.mple.sc[3]<- 500*time.mple/10+ time.sim
boot.mple.sc[4]<- 500*time.mple/25+ time.sim
boot.mple.sc[5]<- 500*time.mple/50+ time.sim
boot.mple.sc[6]<- 500*time.mple/100+ time.sim
boot.mple.sc[7]<- 500*time.mple/250+ time.sim
boot.mple.sc[8]<- 500*time.mple/500+ time.sim

#relative time

rel.time.sc<- boot.mple.sc/time.mcmle



# relative time cosponsorship


time.mple<- 3.245
time.mcmle<- 100.12
time.sim <- 20.024 # takes 2 iterations, MCMC samplesize is quadrupled in last iteration

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

legend("topright", inset=c(.05, .01), legend=c("Faux Mesa High","Faux Magnolia High", "Cosponsorship"),lty = 1,lwd=2, col = c("blue", "red", "darkgreen"), cex=1.5, bty="n")


################################
## relative time less cores
#################################

# comparison plot MCMLE vs MPLE

# mesa high

time.mple<- 0.05
time.mcmle<- 12.59
time.sim<- 1.573 # simulation time

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

time.mple<- 0.51
time.mcmle <- 119.07
time.sim<- 2.49 # simulation time

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


# relative time scnetwork

time.mple<- 3218
time.mcmle<- 1123200
time.sim <- 199680

# calculate bootMPLE time for 500 iterations and different number of cores (1, 5, 10, 25, 50, 100, 250, 500)
boot.mple.sc<-c()
boot.mple.sc[1]<- 500*time.mple/1 + time.sim
boot.mple.sc[2]<- 500*time.mple/2+ time.sim
boot.mple.sc[3]<- 500*time.mple/3+ time.sim
boot.mple.sc[4]<- 500*time.mple/4+ time.sim
boot.mple.sc[5]<- 500*time.mple/5+ time.sim
boot.mple.sc[6]<- 500*time.mple/6+ time.sim
boot.mple.sc[7]<- 500*time.mple/7+ time.sim
boot.mple.sc[8]<- 500*time.mple/8+ time.sim

#relative time

rel.time.sc<- boot.mple.sc/time.mcmle


# plot
plot(1:8, rel.time.mesa, main="Relative Time", xlab="Number of Cores", ylab="Relative Time", type="l", xaxt="n", 
     col="blue", lwd=2, cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=c(0.2,2.3))
axis(1, 1:8, label=1:8, cex.axis=1.5)
lines(1:8, rel.time.magnolia, col="red", lwd=2)
lines(1:8, rel.time.sc, col="darkgreen", lwd=2)
abline(h=1,col="grey", lty=2, lwd=2)

legend("topright", inset=c(.05, .01), legend=c("Faux Mesa High","Faux Magnolia High", "Supreme Court"),lty = 1,lwd=2, col = c("blue", "red", "darkgreen"), cex=1.5, bty="n")
# 7x10 inches