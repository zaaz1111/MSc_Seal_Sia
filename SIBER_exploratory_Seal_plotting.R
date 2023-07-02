##Use librarian or pacman on startup to manage all packages 
##Figure out how to use a package on startup
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse)

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

mixdat_consumer <- write.csv(sidat,file='Consumer_mixing_data.csv')
#Make subsets for Plasma and RBC
sidatP<-subset(sidat[sidat$Type=='Plasma',])
sidatRBC<-subset(sidat[sidat$Type=='RBC',])

#Subset the species data to see if the sample types are significantly different
sidatharb<-subset(sidat[sidat$Species=='Harbor',])
sidatgrey<-subset(sidat[sidat$Species=='Grey',])

#Transform the data to SIBERdat, data that matches SIBER's expected format:
siberdat<-subset(sidat, select = c(δ15N..air, δ13C...V.PDB, Species, Type))
colnames(siberdat)<-c('iso2','iso1','group','community')
siberdat<-siberdat[,c(2,1,3,4)]%>%
  createSiberObject()

##Redoing the SIBER stuff with the ggplot syntax from https://github.com/AndrewLJackson/SIBER/blob/306aa9dc922b9c73dc8223f423e667f11fad844a/vignettes/Plot-posterior-ellipses.R
#Summary stats (TA, SEA and SEAc) for each group
group.ML <- groupMetricsML(siberdat)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

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

# I think (?) this is the distribution of ellipse areas from the MCMCs??
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# how many of the posterior draws do you want?
n.posts <- 10

# decide how big an ellipse you want to draw
p.ell <- 0.95

# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- bind_rows(all_ellipses, .id = "id")

# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

#Begin plotting the SIBER data with each ellipse drawn
ggplot(data = sidat, aes(δ13C...V.PDB, δ15N..air)) +
  geom_point(aes(color = factor(Species):factor(Type)), size = 2)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=15))+
  theme_minimal() + #Optionally, facet wrap them
  #facet_wrap(~factor(Species):factor(Type))+
  geom_polygon(data = ellipse_df,
               mapping = aes(iso1, iso2,
                             color = factor(group):factor(community),
                             fill = NULL),
               fill = NA,
               alpha = 0.2)+
  guides(color=guide_legend(title='Species:Sample Type'))
  

#Now do it but make it look pretty (These ellispes are normally distributed and thus not a 100% 
#reflection of our ellipses, but they'll give us a nice idea and look spicyyyyyy)
#This manipulates ellipse size
p.ell <- 0.95

ggplot(data = sidat, 
       aes(x = δ13C...V.PDB, 
           y = δ15N..air)) + 
  geom_point(aes(color = Species, shape = Type), size = 3,alpha=0.75) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16))+
  stat_ellipse(aes(group = interaction(Species, Type), 
                   fill = Species, 
                   color = Species), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon")+
  theme_minimal()

##ANOVA the d13C and d15N
onewayC_P<-aov(δ13C...V.PDB~Species,data=sidatP)

onewayC_RBC<-aov(δ13C...V.PDB~Species,data=sidatRBC)

onewayN_P<-aov(δ15N..air~Species,data=sidatP)

onewayN_RBC<-aov(δ15N..air~Species,data=sidatRBC)

onewayC_harb<-aov(δ13C...V.PDB~(Type),data=sidatharb)

onewayC_grey<-aov(δ13C...V.PDB~(Type),data=sidatgrey)

onewayN_harb<-aov(δ15N..air~(Type),data=sidatharb)

onewayN_grey<-aov(δ15N..air~(Type),data=sidatgrey)
