##Use librarian or pacman on startup to manage all packages 
##Figure out how to use a package on startup
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,coda,readxl,patchwork,lmerTest)

#Read in the data
sealisos <- read.csv(here('Processed Data/Seal_Data.csv'))%>%
  mutate(otherID=case_when(
    substr(sealisos$name, 1, 1) == 'D' ~ substr(sealisos$name, 1, 4),
    substr(sealisos$name, 1, 1) == '0' ~ paste(substr(sealisos$name, 1, 3),substr(sealisos$name, 5, 5),sep='')),
    Type=case_when(
      substr(sealdata_big$name, 6, 6) == 'R' ~ 'RBC',
      substr(sealdata_big$name, 6, 6) == 'P' ~ 'Plasma',
      substr(sealdata_big$name, 5, 5) == 'P' ~ 'Plasma'
    ))

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

sidatRBC <- sidat[sidat$Type=='RBC',]

sidatP <- sidat[sidat$Type=='Plasma',]

#Transform the data to SIBERdat, data that matches SIBER's expected format:
siberdat<-subset(sidat, select = c(δ15N..air, δ13C...V.PDB, Species, Type))

#Rename the columns so SIBER doesn't cry
colnames(siberdat)<-c('iso2','iso1','community','group')
siberdat<-siberdat[,c(2,1,4,3)]
write.csv(file='SIBER_Raw_Data.csv',siberdat)
siberdat <- createSiberObject(siberdat)

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
ellipses.posterior.seal <- siberMVN(siberdat, parms, priors)


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


gml<-group.ML#%>%
  data.frame()%>%
  rbind(colMeans(SEA.B))
rownames(gml)[4]<-'SEAb'

#DevinCode
#CHECKING MODEL CONVERGENCE
#you must set wd to same spot the saved posterior draw files are
#   for me this is desktop/belugamanuscript (approx) / Siber_convergence
# get a list of all the files in the save directory
all.files <- dir(parms$save.dir, full.names = TRUE)

# find which ones are jags model files
model.files <- all.files[grep("jags_output", all.files)]

# test convergence for the first one
do.this <- 1
#These need to be run separately for 1:nchains. Could loop through all outputs but with the interactive
#plotting functions separate is easier
load(model.files[do.this])

gelman.diag(output, multivariate = FALSE)
gelman.plot(output, auto.layout = TRUE) #auto.layout = TRUE for 6 in one fig

geweke.diag(output, frac1=0.1, frac2=0.5)
geweke.plot(output, frac1=0.1, frac2=0.5)
traceplot(output)

#output is whatever file selected above (do.this <- 3, for example). Need to keep changing
autocorr.diag(output) #using coda to calculate autocorrelation values
autocorr.plot(output) #using coda to plot autocorrelation 
acfplot(output)
effectiveSize(output) #calculate effective sample size (says 40,000, which is total output length...)
# I think (?) this is the distribution of ellipse areas from the MCMCs??
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

#Cute little dataframe to ggplot above
SEA.B <- siberEllipses(ellipses.posterior.seal)

sillygoose <- rep('Grey | Plasma',5000)
e<-rep('Harbor | Plasma',5000)
f<-rep('Grey | RBC',5000)
g<-rep('Harbor | RBC',5000)
m<-c(sillygoose,e,f,g)

sea_dumyframe <- data.frame(c(SEA.B[,1],SEA.B[,2],SEA.B[,3],SEA.B[,4]))%>%
  cbind(m)
colnames(sea_dumyframe)<-c('data','Group')

sea_dumyframe <- mutate(sea_dumyframe, Type = case_when(substr(Group, nchar(Group), nchar(Group)) == 'a' ~ 'Plasma',
                                                        substr(Group, nchar(Group), nchar(Group)) == 'C' ~ 'RBC'),
                        Species = case_when(substr(Group, 1, 1) == 'G'~'Grey',
                                            substr(Group, 1, 1) == 'H'~'Harbour'))

dodge = position_dodge(width = .9)
                        
ggplot(sea_dumyframe)+
  geom_violin(aes(x=Type,y=data,fill=Species), col = 'black', position = dodge)+
  geom_boxplot(width = 0.2, outlier.shape=NA, aes(x=Type,y=data,fill=Species),col='black', position = dodge)+
  scale_fill_manual(values = c('#008837','#7b3294'))+
  ylab(expression("Standard Ellipse Area " ('\u2030' ^2)))+
  xlab('Sample Type')+
  theme_minimal()+
  theme(axis.text = element_text(size=10, color = 'black'),
        axis.text.x = element_text(vjust = -2),
        axis.title = element_text(size=12, color = 'black'),
        legend.text = element_text(size=10, color = 'black'))
ggsave(filename = 'SEA.B Violins.jpeg')
#Same thing with violins in _

ggplot(sea_dumyframe)+
  geom_violin(aes(x=factor(Group),y=data,fill=factor(Group)),col='black')+
  geom_boxplot(width=0.1, outlier.shape=NA, aes(x=factor(Group),y=data,fill=factor(Group)))+
  scale_fill_manual(name = 'Species | Sample Type', values=scheme1 ,guide_legend(title='Sample Type | Species'), labels = c('Grey | Plasma',
                                                                                            'Harbour | Plasma',
                                                                                            'Grey | RBC',
                                                                                            'Harbour | RBC'))+
  scale_x_discrete(labels = c('Grey | Plasma',
                              'Harbour | Plasma',
                              'Grey | RBC',
                              'Harbour | RBC'))+
  ylab(expression("Standard Ellipse Area " ('\u2030' ^2) ))+
  xlab(NULL)+
  facet_wrap(~Type)+
  theme_minimal()

scheme1 <- c('#a6dba0','#c2a5cf','#008837','#7b3294')
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
                   color = Species), 
               fill = NA, 
               level = 0.95,
               type = "norm",
               geom = "polygon",
               lty = 2) + 
  stat_ellipse(aes(Species = interaction(Species, Type), 
                   color = Species, fill = Species), 
               alpha = 0.25,
               level = 0.40,
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
  
siber.object=siberdat
#Calculate overlap:
overlap_p <- maxLikOverlap(ellipse1='1.2',ellipse2='2.2',siber.object, p = NULL)

coords.1 <- addEllipse(siberdat$ML.mu[[1]][, , 1], 
                       siberdat$ML.cov[[1]][, , 1], m = siberdat$sample.sizes[1, 
                                                                                          1], small.sample = TRUE, n = 100, p.interval = 0.95, 
                       ci.mean = FALSE, do.plot = FALSE)
area.1 <- hullArea(coords.1[, 1], coords.1[, 2])

siber.example <- createSiberObject(demo.siber.data) 
ellipse1 <- "1.1" 
ellipse2 <- "1.2"

ellipse95.overlap <- maxLikOverlap('2.2', '1.1', siberdat,p.interval = 0.95)
                                   
ellipses.posterior$RBC.Grey<-ellipses.posterior$RBC.Grey[15001:30000,]
ellipses.posterior$RBC.Harbor<-ellipses.posterior$RBC.Harbor[15001:30000,]
ellipses.posterior$Plasma.Grey<-ellipses.posterior$Plasma.Grey[15001:30000,]
ellipses.posterior$Plasma.Harbor<-ellipses.posterior$Plasma.Harbor[15001:30000,]

SEA.B <- siberEllipses(ellipses.posterior)
ay <- bayesianOverlap(1.1,2.1, ellipses.posterior, draws=NULL, p.interval = NULL, n=100,do.plot=T)
ayprop <- (ay[,3] / (ay[,2] + 
                                                ay[,1] -
                                                ay[,3])
)






maybe<-bayesianOverlap(1.1,1.2, ellipses.posterior, draws=10, p.interval = 0.95, n=100,do.plot=T)
maybeprop <- (maybe[,3]/(maybe[,2]+maybe[,1]-maybe[,3]))

bayesianOverlap(1.1,2.1, ellipses.posterior, draws=10, p.interval = 0.90, n=100,do.plot=T)
bayesianOverlap(1.1,2.1, ellipses.posterior, draws=10, p.interval = 0.50, n=100,do.plot=T)







st<-Sys.time()
bayes_overlap_p <- bayesianOverlap(1.1,2.1, SEA.B, draws=10, p.interval = 0.95,do.plot=T)
et <- Sys.time()
et-st

hist(bayes_overlap_p[,3], 10)

bayes_overlap_r <- bayesianOverlap(ellipse1=1.2,ellipse2=2.2,ellipses.posterior,draws=100,n=360,p.interval=0.95)


bayes.prop.95.over <- (bayes_overlap_p[,3] / (bayes_overlap_p[,2] + 
                                                bayes_overlap_p[,1] -
                                                bayes_overlap_p[,3])
)

ggplot(sidat)+
  geom_point(aes(x=theDateTime,y=δ15N..air,col=location,shape=sex))+
  theme_minimal()

sidatlmmN <- lmerTest::lmer(δ15N..air ~  sex + Type + (1|otherID), data=sidat[sidat$species=='Hg',])
sidatlmmC <- lmerTest::lmer(δ13C...V.PDB ~ Type+(1|otherID), data=sidat[sidat$species=='Pv',])
sidatlmm <- lmerTest::lmer(δ13C...V.PDB ~ δ15N..air+sex+species+Type+region+(1|otherID), data=sidat)

ellipse1=1.2
ellipse2 <- 2.2
coords.1 <- addEllipse(ellipses.posterior[[ellipse1]][i, 
                                                      5:6], matrix(ellipses.posterior[[ellipse1]][i, 1:4], 
                                                                   nrow = 2, ncol = 2), p.interval = 0.95, n = 100, 
                       do.plot = T, small.sample = FALSE)
area.1 <- siberConvexhull(coords.1[, 1], coords.1[, 2])
coords.2 <- addEllipse(ellipses.posterior[[ellipse2]][i, 
                                                      5:6], matrix(ellipses.posterior[[ellipse2]][i, 1:4], 
                                                                   nrow = 2, ncol = 2), p.interval = 0.95, n = 100, 
                       do.plot = F, small.sample = FALSE)
area.2 <- siberConvexhull(coords.2[, 1], coords.2[, 2])

overlap <- abs(spatstat.utils::overlap.xypolygon(list(x = coords.1[, 
                                                                   1], y = coords.1[, 2]), list(x = coords.2[, 1], y = coords.2[, 
                                                                                                                                2])))

#Plot the siber ellipses idfkanymore


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
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = 0.95)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- bind_rows(all_ellipses, .id = "id")
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)
first.plot <- ggplot(data = sidat, aes(δ13C...V.PDB,δ15N..air)) +
  geom_point(aes(color = factor(Type):factor(Species)), size = 2)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=15))
print(first.plot)

second.plot <- first.plot + facet_wrap(~factor(Type):factor(Species))
second.plot

colnames(ellipse_df)<-c('id','iso1','iso2','rep','Species','Type')
third.plot <- second.plot + 
  geom_polygon(data = ellipse_df,
               mapping = aes(iso1, iso2,
                             group = rep,
                             color = factor(Type):factor(Species),
                             fill = NULL),
               fill = NA,
               alpha = 0.2)
print(third.plot)


