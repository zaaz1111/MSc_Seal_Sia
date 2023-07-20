fish.samples <- read.csv(here('fish.samples.csv'))%>%
  dplyr::select(name:Scientific.Name)
fish.data <- read.csv(here('Processed Data/fish.samples.final.csv'))
mergedfish <- merge(fish.samples,fish.data)
#Make a dataframe with the N from lipid rich samples, and the C from Lipid Extracted samples
fishsourcec <- subset(mergedfish,mergedfish$Lipid.Extracted.,
                      select=c(name,Spp,δ13C))#%>%
fishsourcec <- mutate(fishsourcec,name = substr(fishsourcec$name,1,nchar(fishsourcec$name)-2))

fishsourcen<-subset(mergedfish,!mergedfish$Lipid.Extracted.,
                    select=c(name,δ15N))

fishsourcecomp <- merge(fishsourcec,fishsourcen,by='name')%>%
  subset(fishsourcecomp$δ15N>0)

write.csv(fishsourcecomp,file='Composite_source_data.csv')
###Add in sandeel by hand
####Silly stupid dumb fish re-cleaning hehe xD
Flatfish <- c('Lemon Sole', 'Long Rough Dab','Common Dab')
Scorpionfish <- c('Red Gurnard','Grey Gurnard')
Trisopterus <- c('Blue Whiting','Norway Pout','Poor Cod')
#Reload and pipe data to useable format
fishsource.new <- read.csv(here('Processed data/Composite_source_data.csv'))%>%
  mutate(group = case_when(
    fishsource.new$Spp %in% Flatfish ~ 'Flatfish',
    fishsource.new$Spp %in% Scorpionfish ~ 'Scorpionfish',
    fishsource.new$Spp %in% Trisopterus ~ 'Trisopterus',
    TRUE ~ fishsource.new$Spp
  ))

write.csv(fishsource.new,file='idfkanymore..csv')


unique(fishsource.new$group)
for(i in unique(fishsource.new$group)){
  #hist(fishsource.new$δ13C[fishsource.new$group==i])
  hist(fishsource.new$δ15N[fishsource.new$group==i])
}

###Unioning Fish
#Flatfish
kw_flatc<-kruskal.test(δ13C ~ Spp,
               data=fishsource.new[fishsource.new$group=='Flatfish',])
kw_flatc

kw_flatn<-kruskal.test(δ15N ~ Spp,
               data=fishsource.new[fishsource.new$group=='Flatfish',])
kw_flatn
#Check on FishBase...
#We can union LRD and Lemon sole, and Common Dab


#Scorpionfish
kw_scorpc<-kruskal.test(δ13C ~ Spp,
               data=fishsource.new[fishsource.new$group=='Scorpionfish',])
kw_scorpc

kw_scorpn<-kruskal.test(δ15N ~ Spp,
                data=fishsource.new[fishsource.new$group=='Scorpionfish',])
kw_scorpn
#Can union Red and Grey gurnard!!
#But NOT Bluemouth

#Trisopterus
kw_tric<-kruskal.test(δ13C ~ Spp,
                data=fishsource.new[fishsource.new$group=='Trisopterus',])
kw_tric

kw_trin<-kruskal.test(δ15N ~ Spp,
              data=fishsource.new[fishsource.new$group=='Trisopterus',])
kw_trin
###Unioning Blue Whiting, Poor Cod, and Norway Pout
#It's shaky



  
a <-ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ13C,fill=Spp))+
  theme_minimal()

b <- ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ15N,fill=Spp))+
  theme_minimal()

a+b


c <- ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ13C,fill=group))+
  theme_minimal()

d <- ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ15N,fill=group))+
  theme_minimal()

(a+b)/(c+d)


###Fish Isospace Plotting
sbs <- read.csv(here('Processed data/Composite_source_data.csv'))%>%
  group_by(Spp)%>%
  summarize(count = n(),
            mC = mean(δ13C),
            sdC= sd(δ13C),
            mN = mean(δ15N), 
            sdN = sd(δ15N))


ggplot()+
  geom_errorbar(data = sbs, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN, 
                              ymax = mN + 1.96*sdN), 
                width = 0)+
  geom_errorbarh(data = sbs, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0) + 
  geom_point(data = sbs, aes(x = mC, 
                             y = mN,
                             fill = Spp),
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = T) + 
  scale_fill_viridis_d()+
  theme_minimal()+
  geom_point(data=sidatP, aes(x = δ13C...V.PDB,
                         y = δ15N..air,
                         col=Species))+
  scale_color_manual(values = c('#7b3294','#008837'))