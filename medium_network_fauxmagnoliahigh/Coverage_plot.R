library(statnet)
# plot coverage probability plots next to each other

par(mfrow=c(1,2))

# Mesa High
load("Coverage_mesahigh2.RData")

# coef=Inf appears
faux.results.mple[[1]][!is.finite(faux.results.mple[[1]])] <- NA
faux.results.mple[[2]][!is.finite(faux.results.mple[[2]])] <- NA
faux.results.mple[[3]][!is.finite(faux.results.mple[[3]])] <- NA
faux.results.mple[[4]][!is.finite(faux.results.mple[[4]])] <- NA
faux.results.mple[[5]][!is.finite(faux.results.mple[[5]])] <- NA


## Coverage probability plot

# create empty matrix for barplot values
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
barplot(counts, beside=TRUE, main="Coverage Probability Faux Mesa High", ylim=c(0,1),
        legend = rownames(counts),args.legend = list(x="bottomleft", cex=1.3), cex.axis=1.5,
        cex.names = 2, cex.main=1.5)
abline(h=0.95,col="grey", lty=2, lwd=2)


# Magnolia High
load("Coverage_magnoliahigh2.RData")


# coef=Inf appears
faux.results.mple[[1]][!is.finite(faux.results.mple[[1]])] <- NA
faux.results.mple[[2]][!is.finite(faux.results.mple[[2]])] <- NA
faux.results.mple[[3]][!is.finite(faux.results.mple[[3]])] <- NA
faux.results.mple[[4]][!is.finite(faux.results.mple[[4]])] <- NA
faux.results.mple[[5]][!is.finite(faux.results.mple[[5]])] <- NA


## Coverage probability plot

# create empty matrix for barplot values
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
barplot(counts, beside=TRUE, main="Coverage Probability Faux Magnolia High", ylim=c(0,1),
        #legend = rownames(counts),args.legend = list(x="bottomleft", cex=1.3), 
        cex.axis=1.5, cex.names = 2, cex.main=1.5)
abline(h=0.95,col="grey", lty=2, lwd=2)



####### Bias plots

par(mfrow=c(1,2))
# Faux Mesa High

# Mesa High
load("Coverage_mesahigh2.RData")

# coef=Inf appears
faux.results.mple[[1]][!is.finite(faux.results.mple[[1]])] <- NA
faux.results.mple[[2]][!is.finite(faux.results.mple[[2]])] <- NA
faux.results.mple[[3]][!is.finite(faux.results.mple[[3]])] <- NA
faux.results.mple[[4]][!is.finite(faux.results.mple[[4]])] <- NA
faux.results.mple[[5]][!is.finite(faux.results.mple[[5]])] <- NA

# get data
data(faux.mesa.high)

# calculate the true values
mple <- try(ergm(faux.mesa.high ~ edges +nodematch("Sex")+gwesp(0,fixed=TRUE), estimate = "MPLE" )) #

boxplot(faux.results.mple[[1]][1,]-mple$coef[1], faux.results.mple[[1]][2,]-mple$coef[2], 
        faux.results.mple[[1]][3,]-mple$coef[3],outline=FALSE,
        names=c("Edges", "Sex" ,"Gwesp"), main="Bias Faux Mesa High",
        cex.axis=1.8, cex.names = 2.5, cex.main=1.8)
abline(a=0, b=0, col="grey50", lty=2, lwd=3)


# Magnolia High
load("Coverage_magnoliahigh2.RData")

# coef=Inf appears
faux.results.mple[[1]][!is.finite(faux.results.mple[[1]])] <- NA
faux.results.mple[[2]][!is.finite(faux.results.mple[[2]])] <- NA
faux.results.mple[[3]][!is.finite(faux.results.mple[[3]])] <- NA
faux.results.mple[[4]][!is.finite(faux.results.mple[[4]])] <- NA
faux.results.mple[[5]][!is.finite(faux.results.mple[[5]])] <- NA

# get data
data(faux.magnolia.high)

# calculate the true values
mple <- try(ergm(faux.magnolia.high ~ edges +nodematch("Sex")+gwesp(0,fixed=TRUE), estimate = "MPLE" )) #

boxplot(faux.results.mple[[1]][1,]-mple$coef[1], faux.results.mple[[1]][2,]-mple$coef[2], 
        faux.results.mple[[1]][3,]-mple$coef[3], outline=FALSE,
        names=c("Edges", "Sex" ,"Gwesp"), main="Bias Faux Magnolia High",
        cex.axis=1.8, cex.names = 2.5, cex.main=1.8)
abline(a=0, b=0, col="grey50", lty=2, lwd=3)
