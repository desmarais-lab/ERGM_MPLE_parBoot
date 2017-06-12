#load(file="supreme_additional_data.RData")

library(ergm)

sponsorshipMatrix <- (as.matrix(read.csv("./cosponsorship2010/house_matrices/108_housematrix.txt",stringsAsFactors=F,header=F))==1)

bills <- read.csv("./cosponsorship2010/house_bills.txt",stringsAsFactors=F,header=F)

bills <- bills[which(bills[,1]==108),]


sponsorParty <- rep(NA,ncol(sponsorshipMatrix))

members <- read.csv("./cosponsorship2010/house_members/108_house.txt",stringsAsFactors=F,header=F)

memberData <- read.csv("SH.csv",stringsAsFactors=F)
memberData <- memberData[memberData$congress==108,]

for(i in 1:length(sponsorParty)){
    idSponsori <- as.numeric(gsub(" ","",members[which(sponsorshipMatrix[,i]),3]))
    if(length(idSponsori)>0){
        sponsorParty[i] <- memberData$party[match(idSponsori,memberData$ids)]
    }
    
}

cosponsorshipMatrix <- (as.matrix(read.csv("./cosponsorship2010/house_matrices/108_housematrix.txt",stringsAsFactors=F,header=F))==2)


allSponsorMat <- sponsorshipMatrix + cosponsorshipMatrix

cosponsors <- apply(allSponsorMat,2,sum)

corMat <- cor(allSponsorMat)

remove <- which(is.na(sponsorParty))
remove <- c(remove,which(cosponsors<6),which(bills[,2]!="HR"))

sponsorParty <- sponsorParty[-remove]
corMat <- corMat[-remove,-remove]

set.seed(123)
threshMat <- matrix(0.225+0.775*runif(nrow(corMat)^2),nrow(corMat),nrow(corMat))

cosponsorNetwork <- network(1*(corMat>threshMat),dir=F)

set.vertex.attribute(cosponsorNetwork,"sponsorParty",sponsorParty)


# get MPLE
set.seed(5555)
mple<- ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), estimate="MPLE")
summary(mple)

# simulate 500 networks
set.seed(1234)
sim.mple<- simulate.ergm(mple, nsim=500,  statsonly = TRUE, control=control.simulate.ergm(MCMC.interval=200000))
# get observed statistics by simulating 25 times with MCMC.interval=MCMC.brunin=1. Then, take the mode

d<- dim(sim.mple)[2] # save number of columns

# calculate observed statistics by simulating 25 networks, each with burnin=1. In the next step we take the mode and treat is as the observed statistic
set.seed(1)
net.stats <- numeric(d)
A <- matrix(0,25,d)
for(i in 1:25){
    set.seed(i*1000)
    simulated.network = simulate.ergm(mple, nsim=1, control=control.simulate.ergm(MCMC.burnin=1, MCMC.interval=1), statsonly = TRUE)
    A[i,]<- simulated.network
}
for(i in 1:d){
    net.stats[i] <- names(sort(-table(A[,i])))[1]
}

par(mfrow=c(2,2))
for(i in 1:3){
    plot(sim.mple[,i], ylab="Value")
    abline(h=net.stats[i])
}

### Code from here on needs to be revised ###

# get mode and save observed stats
edges.obs<- net.stats[1]
edges.obs
sponsorParty.obs<- net.stats[2]
sponsorParty.obs
altkstar.obs<- net.stats[3]
altkstar.obs



# plotting
pdf(file="MPLEDegeneracyCheck.pdf",width=4,height=8)
par(mfrow=c(3,2))

plot(sim.mple[,1], ylab="Value", type="l", main="Edges")
abline(h=net.stats[1])

hist(sim.mple[,1], main="Edges", xlab="Value") 
abline(v=edges.obs, col='red')

plot(sim.mple[,2], ylab="Value", type="l", main="Sponsor Party")
abline(h=net.stats[2])

hist(sim.mple[,2], main="Sponsor Party", xlab="Value") 
abline(v=sponsorParty.obs, col='red')

plot(sim.mple[,3], ylab="Value", type="l", main="Alternating K-Star")
abline(h=net.stats[3])

hist(sim.mple[,3], main="Alternating K-Star", xlab="Value") 
abline(v=altkstar.obs, col='red')


dev.off()

save.image(file="MPLE_add_data_hist_half.RData")


# now MCMLE

set.seed(5555)
mcmle<- ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), eval.loglik = FALSE, control = control.ergm(MCMC.interval = 30000, MCMC.burnin = 300000))
summary(mcmle)

mcmc.diagnostics(mcmle)



# MPLE runtime

time.mple<- c()
set.seed(1)
for(i in 1:100){
  b<- system.time(ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), estimate="MPLE"))
  time.mple[i]<- b[3]                  
}

median(time.mple)


# MCMLE runtime

time.mcmle<- c()
set.seed(1)
for(i in 1:5){
  b<- system.time(ergm(cosponsorNetwork~edges+nodematch("sponsorParty")+altkstar(0.4975,fixed=T), 
                       eval.loglik = FALSE, control = control.ergm(MCMC.interval = 30000, MCMC.burnin = 300000)))
  time.mcmle[i]<- b[3]                  
}

median(time.mcmle)

