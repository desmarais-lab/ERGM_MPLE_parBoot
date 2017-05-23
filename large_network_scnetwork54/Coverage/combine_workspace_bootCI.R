# loading different workspaces

load("supreme54.RData")

load("coverage_parallel_scnetwork1.RData")
boot1<- boot

load("coverage_parallel_scnetwork2.RData")
boot2<- boot

load("coverage_parallel_scnetwork3.RData")
boot3<- boot

load("coverage_parallel_scnetwork4.RData")
boot4<- boot

rm(boot)

# combine results
boot<-cbind(boot1, boot2, boot3, boot4)
rm(boot1,boot2, boot3, boot4)

# bootstrap.CI: if coef(mple) is in BootCI, save 1 otherwise 0
# create empty matrix
bootstrap.CI <- matrix(0,7,2)
colnames(bootstrap.CI)<- c("upper", "lower")
rownames(bootstrap.CI)<- c("edges", "nodecov.salience", "nodecov.year.center","absdiff.year.center", "absdiff.year.center.square", "dsp0", "esp0")

# calculate the bootstrap intervals
bootstrap.CI[1,]<-quantile(boot[1,], probs = c(0.975, 0.025)) # edges
bootstrap.CI[2,]<-quantile(boot[2,], probs = c(0.975, 0.025)) # salience
bootstrap.CI[3,]<-quantile(boot[3,], probs = c(0.975, 0.025)) # nodecov.year
bootstrap.CI[4,]<-quantile(boot[4,], probs = c(0.975, 0.025)) # absdiff.year
bootstrap.CI[5,]<-quantile(boot[5,], probs = c(0.975, 0.025)) # absdiff.year.square
bootstrap.CI[6,]<-quantile(boot[6,], probs = c(0.975, 0.025)) # dsp0
bootstrap.CI[7,]<-quantile(boot[7,], probs = c(0.975, 0.025)) # esp0

mple
bootstrap.CI
