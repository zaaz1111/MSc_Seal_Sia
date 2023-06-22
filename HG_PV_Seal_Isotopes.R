##Use librarian or pacman on startup to manage all packages 
##Figure out how to use a package on startup
library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2)

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
  xlab('δ13C')+
  ylab('δ15N')+
  theme_minimal()



