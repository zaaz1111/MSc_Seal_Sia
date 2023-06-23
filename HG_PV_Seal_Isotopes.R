##Use librarian or pacman on startup to manage all packages 
##Figure out how to use a package on startup
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags)

#Read in the data
sidat<-read.csv(here('Processed Data/Exploratory_plate_1_Data.csv'),head=T)

#Mutate it to add species name and sample type, and pipe to subset out error rows
sidat<-mutate(sidat,Species=case_when(
         substr(sidat$name, 1, 1) == 'D' ~ 'Grey',
         substr(sidat$name, 1, 1) == '0' ~ 'Harbor'),
         Type=case_when(
           substr(sidat$name, 6, 6) == 'R' ~ 'RBC',
           substr(sidat$name, 6, 6) == 'P' ~ 'Plasma'
         ))%>%
  subset(sidat$C.N.ratio>2)

##Doing SIBER stuff from CRAn
siberdat<-subset(sidat,select=c("δ15N..air","δ13C...V.PDB",'Species','Type'))
colnames(siberdat)<-c("iso2", "iso1", "community", "group")
siberdat<-siberdat[, c(2,1,4,3)]%>%
  createSiberObject()

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")


plotSiberObject(siberdat,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
                legend=legend('topright',inset=0.025,
                       legend=c('Grey Seal','Harbor Seal','RBC','Plasma'),
                       col=c('red','black','grey50','grey50'),pch=c(4,4,17,19),pt.cex=2))
)

group.ML <- groupMetricsML(siberdat)%>%
  data.frame()
print(group.ML)

###Starting a bayesian model?!?!?
parms<-list()
parms$n.iter<-2*10^4
parms$n.burnin<-1*10^3
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2 

priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

#Fitting the ellipses using the specified above parameters (This uses JAGS, something I am still learning)
ellipses.posterior <- siberMVN(siberdat, parms, priors)

SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)



#plot the whole thing
ggplot(sidat)+
  geom_point(aes(δ13C...V.PDB,δ15N..air,color=Species,shape=Type))+
  xlab('δ13C')+
  ylab('δ15N')+
  theme_minimal()

#Exploratory plotting for Plasma
ggplot(sidat[sidat$Type=='Plasma',])+
  geom_point(aes(δ13C...V.PDB,δ15N..air,color=Species))+
  stat_ellipse(aes(δ13C...V.PDB,δ15N..air,group = interaction(Species), 
                   fill = Species, 
                   color = Species), 
               alpha = 0.25, 
               level = 0.95,
               type = "norm",
               geom = "polygon") + 
  xlab('δ13C')+
  ylab('δ15N')+
  theme_minimal()

#Exploratory plotting for RBCs
ggplot(sidat[sidat$Type=='RBC',])+
  geom_point(aes(δ13C...V.PDB,δ15N..air,color=Species))+
  stat_ellipse(aes(δ13C...V.PDB,δ15N..air,group = interaction(Species), 
                   fill = Species, 
                   color = Species), 
               alpha = 0.25, 
               level = 0.95,
               type = "norm",
               geom = "polygon") + 
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030'))+
  theme_minimal()

###Processing fish sampling data
fs<-matrix(nrow=70,ncol=7)
colnames(fs)<-c('Sample Name','Lipid Extracted?','Spp','Tow','Replicate','Plate','Well')
fs<-data.frame(fs)
for(i in 1:nrow(fs)){
  fs$Sample.Name[i]<-paste('FS',i,'LE',sep='')
  fs$Lipid.Extracted.<-TRUE
}

write.csv(fs,file='fish.samples.csv')
fs<-read.csv("C:/Users/zaahi/Documents/fish.samples.csv")#%>%

mutate(fs,Sample.Name=case_when(
  fs$Lipid.Extracted. == T ~ fs$Sample.Name,
  fs$Lipid.Extracted. == F ~ str_sub(fs$Sample.Name, 1, -3)
  ))
