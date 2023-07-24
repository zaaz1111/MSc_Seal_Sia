##Use librarian or pacman on startup to manage all packages 
##Figure out how to use a package on startup
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,coda,readxl,patchwork)

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

#Rename the columns so SIBER doesn't cry
colnames(siberdat)<-c('iso2','iso1','community','group')
siberdat<-siberdat[,c(2,1,4,3)]%>%
  createSiberObject()

colnames(siberdatp)<-c('iso2','iso1','community','group')
siberdatp<-siberdatp[,c(2,1,4,3)]%>%
  createSiberObject()

colnames(siberdatRBC)<-c('iso2','iso1','group','community')
siberdatRBC<-siberdatRBC[,c(2,1,4,4)]%>%
  createSiberObject()

##Redoing the SIBER stuff with the ggplot syntax from https://github.com/AndrewLJackson/SIBER/blob/306aa9dc922b9c73dc8223f423e667f11fad844a/vignettes/Plot-posterior-ellipses.R
#Summary stats (TA, SEA and SEAc) for each group
group.ML <- groupMetricsML(siberdat)#
community.ML <- communityMetricsML(siberdat)
write.csv(group.ML,file='Seal sample Layman metrics.csv')


# options for running jags
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
ellipses.posterior <- siberMVN(siberdat, parms, priors)


# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

#DevinCode go big mode
#Except I'm transposing it here
mu_post <- extractPosteriorMeans(siberdat, ellipses.posterior)
###mu.post[[a]][b,c,d] *************************************
#where a is community (1 = Plasma 2 = RBC)
#b = row, c = isotope (1 = d13C, 2 = d15N),
#and d = individual seal 
#so mu.post[[1]][,2,] will give posteriors for all seal plasma d15N

layman.B <- bayesianLayman(mu_post)
#In this case, layman.B is a list of all posterior draws for layman metrics 
#for both communities. so layman.B[[1]] is community 1, layman.B[[2]] is
#community 2, and so on.

#Make a dataframe of the group metrics plus the mean of SEAb posterior density
SEA.B<-SEA.B%>%
  data.frame()

gml<-group.ML%>%
  data.frame()%>%
  rbind(colMeans(SEA.B))
rownames(gml)[4]<-'SEAb'


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


# I think (?) this is the distribution of ellipse areas from the MCMCs??
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

#Cute little dataframe to ggplot above
sea_dumyframe <- data.frame(c(SEA.B[,1],SEA.B[,2],SEA.B[,3],SEA.B[,4]))%>%
  cbind(m)
colnames(sea_dumyframe)<-c('z','m')
sillygoose <- rep('Grey | Plasma ',6000)
e<-rep('Harbor | Plasma ',6000)
f<-rep('Grey | RBC',6000)
g<-rep('Harbor | RBC',6000)
m<-c(sillygoose,e,f,g)

#Same thing with violins in ggplot
ggplot(sea_dumyframe)+
  geom_violin(aes(x=factor(m),y=z,fill=factor(m)),col='black')+
  geom_boxplot(width=0.1, outlier.shape=NA, aes(x=factor(m),y=z,fill=factor(m)))+
  scale_fill_manual(values=scheme1 ,guide_legend(title='Sample Type | Species'))+
  ylab(expression("Standard Ellipse Area " ('\u2030' ^2) ))+
  xlab(NULL)+
  theme_minimal()

scheme1 <- c('#c2a5cf','#7b3294','#a6dba0','#008837')
scheme2 <- c('#b2abd2','#5e3c99','#fdb863','#e66101')
scheme3 <- c('#dfc27d','#a6611a','#80cdc1','#018571')
scheme4 <- c('#92c5de','#0571b0','#f4a582','#ca0020')
  


#Plot the 95%CI ellipses
ggplot(data = sidat, aes(x = δ13C...V.PDB, y = δ15N..air)) + 
  geom_point(aes(color = Species, shape = Species), size = 2) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  stat_ellipse(aes(Species = interaction(Species, Type), 
                   fill = Species, 
                   color = Species), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=c('#008837','#7b3294'))+
  scale_color_manual(values=c('#008837','#7b3294'))+
  ylim(values=c(11,16))+
  xlim(values=c(-20.5,-16))+
  facet_wrap(~Type, scales = 'free')+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 15),
        axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank())
  

#Calculate overlap:
overlap_p <- maxLikOverlap('1.1', '1.2', siberdat, p = 0.95)

st<-Sys.time()
bayes_overlap_p <- bayesianOverlap(1.1,2.1,ellipses.posterior,draws=10,n=360)
et <- Sys.time()
et-st

bayes_overlap_r <- bayesianOverlap(1.2,2.2,ellipses.posterior,draws=10,n=360)
