##Load the packages
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,coda,readxl,patchwork,lmerTest)

#Read in the data
sealisos <- read.csv(here('Processed Data/Seal_Data.csv'))%>%
  mutate(otherID=case_when(
    substr(sealisos$name, 1, 1) == 'D' ~ substr(sealisos$name, 1, 4),
    substr(sealisos$name, 1, 1) == '0' ~ paste(substr(sealisos$name, 1, 3),substr(sealisos$name, 5, 5),sep='')))

sealmeta <- read.csv(here('Raw Data/captureSummary.csv'))

sealdata_big <-merge(sealisos,sealmeta)

#Mutate it to add species name and sample type, and pipe to subset out error rows
sidat<-mutate(sealdata_big,Species=case_when(
  substr(sealdata_big$name, 1, 1) == 'D' ~ 'Grey',
  substr(sealdata_big$name, 1, 1) == '0' ~ 'Harbor'),
  Type=case_when(
    substr(sealdata_big$name, 6, 6) == 'R' ~ 'RBC',
    substr(sealdata_big$name, 6, 6) == 'P' ~ 'Plasma',
    substr(sealdata_big$name, 5, 5) == 'P' ~ 'Plasma'
  ))
sidat<-sidat[sidat$C.N.ratio>2,]
sidat<-sidat[!duplicated(sidat$name),]

#Transform the data to SIBERdat, data that matches SIBER's expected format:
siberdat<-subset(sidat, select = c(δ15N..air, δ13C...V.PDB, Species, Type))

#Rename the columns so SIBER doesn't cry
colnames(siberdat)<-c('iso2','iso1','community','group')
siberdat<-siberdat[,c(2,1,4,3)]
write.csv(file='SIBER_Raw_Data.csv',siberdat)
siberdat <- createSiberObject(siberdat)

###DON'T RUN THESE!!! Just want to know how the models were configured
parms <- list()
parms$n.iter <- 1000000   # number of iterations to run the model for
parms$n.burnin <- 500000 # discard the first set of values
parms$n.thin <- 500     # thin the posterior by this many
parms$n.chains <- 5        # run this many chains
parms$save.output = T #Tell it to save the SIBER output
parms$save.dir <- here('SIBER') # Where to save this model

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior.seal <- siberMVN(siberdat, parms, priors)


sea.b.seal <- siberEllipses(ellipses.posterior.seal)

##Could add the density plot here but don't 100% need it

#Checking to see if the burnin is really needed (we have to do it by hand so I'm interested to see what the differences are, if any)
back.ell <- list()
back.ell$Grey.RBC<-ellipses.posterior.seal$Grey.RBC[(nrow(ellipses.posterior.seal$Grey.RBC)-100):nrow(ellipses.posterior.seal$Grey.RBC),]
back.ell$Harbor.RBC<-ellipses.posterior.seal$Harbor.RBC[(nrow(ellipses.posterior.seal$Harbor.RBC)-100):nrow(ellipses.posterior.seal$Harbor.RBC),]
ellipses.posterior$Plasma.Grey<-ellipses.posterior$Plasma.Grey[15001:30000,]
ellipses.posterior$Plasma.Harbor<-ellipses.posterior$Plasma.Harbor[15001:30000,]

back_o <- bayesianOverlap(1,2,back.ell, draws=NULL,p.interval=NULL,n=100,do.plot=F)
front_o <- bayesianOverlap(1,2,ellipses.posterior.seal, draws= 101, p.interval =NULL, n= 100, do.plot=F)
par(mfrow=c(1,2))
hist(back_o$overlap)
hist(front_o$overlap)

t.test(back_o$overlap,front_o$overlap)
##Looks like burnin can kind of fuck off, but I don't think that's best practices for MCMC
#We'll do it live
ellipses.posterior.seal$Grey.RBC<-ellipses.posterior.seal$Grey.RBC[5001:10000,]
ellipses.posterior.seal$Harbor.RBC<-ellipses.posterior.seal$Harbor.RBC[5001:10000,]
ellipses.posterior.seal$Grey.Plasma<-ellipses.posterior.seal$Grey.Plasma[5001:10000,]
ellipses.posterior.seal$Harbor.Plasma<-ellipses.posterior.seal$Harbor.Plasma[5001:10000,]

#Now here the bayesianOverlap code is lowkey broken, but there are a few fixes
#Set the p.interval to NULL; this computes the standard isotopic ellipse (which ends up being close to 40%)
#And then the groups are like just 1:4 as per a call of ellipses.posterior.seal$ in the console
