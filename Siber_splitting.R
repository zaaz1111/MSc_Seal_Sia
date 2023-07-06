##Use librarian or pacman on startup to manage all packages 
##Figure out how to use a package on startup
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,coda)

#Read in the data
sidat<-read.csv(here('Processed Data/Seal_Data.csv'),head=T)

#Mutate it to add species name and sample type, and pipe to subset out error rows
sidat<-mutate(sidat,Species=case_when(
  substr(sidat$name, 1, 1) == 'D' ~ 'Grey',
  substr(sidat$name, 1, 1) == '0' ~ 'Harbor'),
  Type=case_when(
    substr(sidat$name, 6, 6) == 'R' ~ 'RBC',
    substr(sidat$name, 6, 6) == 'P' ~ 'Plasma',
    substr(sidat$name, 5, 5) == 'P' ~ 'Plasma'
  ))%>%
  subset(sidat$C.N.ratio>2)

#Make subsets for Plasma and RBC
sidatP<-subset(sidat[sidat$Type=='Plasma',])
sidatRBC<-subset(sidat[sidat$Type=='RBC',])

#Subset the species data to see if the sample types are significantly different
sidatharb<-subset(sidat[sidat$Species=='Harbor',])
sidatgrey<-subset(sidat[sidat$Species=='Grey',])

#Transform the data to SIBERdat, data that matches SIBER's expected format:
siberdat<-subset(sidat, select = c(δ15N..air, δ13C...V.PDB, Species, Type))

#Siberdat for Harbor Plasma
siberdat_H_P <- sidat[sidat$Species == 'Harbor' & sidat$Type == 'Plasma',]%>%
  subset(select = c(δ15N..air, δ13C...V.PDB, Species, Type))

#Siberdat for Harbor RBC
group.ML_H_P <- groupMetricsML(siberdat_H_P)

parms <- list()
parms$n.iter <- 1000000   # number of iterations to run the model for
parms$n.burnin <- 500000 # discard the first set of values
parms$n.thin <- 500     # thin the posterior by this many
parms$n.chains <- 3        # run this many chains
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
ellipses.posterior_H_P <- siberMVN(siberdat_H_P, parms, priors)

SEA.B_HP <- siberEllipses(ellipses.posterior)


#Look at the diagnostics using Coda
all.files <- dir(parms$save.dir, full.names = TRUE)

# find which ones are jags model files
model.files <- all.files[grep("jags_output", all.files)]

# test convergence for the first one
do.this <- 1

load(model.files[do.this])

gelman.diag(output, multivariate = FALSE)

par(mfrow=c(3,2))
gelman.plot(output, auto.layout = FALSE)
par(mfrow=c(1,1))

