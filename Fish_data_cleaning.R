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
Scorpionfish <- c('Red Gurnard','Grey Gurnard','Squid')
Small_Benthics <- c('Blue Whiting','Norway Pout','Poor Cod')
#Reload and pipe data to useable format
fishsource.new <- read.csv(here('Processed data/Composite_source_data.csv'))%>%
  mutate(group = case_when(
    fishsource.new$Spp %in% Flatfish ~ 'Flatfish',
    fishsource.new$Spp %in% Scorpionfish ~ 'Scorpionfish',
    fishsource.new$Spp %in% Small_Benthics ~ 'Small Benthics',
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


kw_bc <- kruskal.test(δ13C ~ Spp,
               data=fishsource.new[fishsource.new$group=='Small Benthics',])

kw_bc

kw_bn <- kruskal.test(δ15N ~ Spp,
                        data=fishsource.new[fishsource.new$group=='Small Benthics',])

kw_bn
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

ggplot(mergedfish)+
  geom_boxplot(aes(x=Spp,y=δ13C,fill=Lipid.Extracted.))+
  theme_minimal()

  
a <-ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ13C,fill=Spp))+
  theme_minimal()

b <- ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ15N,fill=Spp))+
  theme_minimal()

a+b


c <- ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ13C,fill=group))+
  scale_fill_manual(values=fish_pal)+
  theme_minimal()

d <- ggplot(fishsource.new)+
  geom_boxplot(aes(y=δ15N,fill=group))+
  scale_fill_manual(values=fish_pal)+
  theme_minimal()

(c+d)


###Fish Isospace Plotting: make a grouped dataframe (this will help with creating our errorbars)
sbs <- read.csv(here('idfkanymore..csv'))%>%
  group_by(group)%>%
  summarize(count = n(),      
            mC = mean(δ13C),
            sdC= sd(δ13C),
            mN = mean(δ15N), 
            sdN = sd(δ15N))

#Add in tdf data (adjusting the means and SDs one by one using mutate)
r_tdfs <- read.csv(here('Processed data/Semi_smart_Seal_RBC_TDFS.csv'))
p_tdfs <- read.csv(here('Processed data/Semi_smart_Seal_Plasma_TDFS.csv'))

#Create the sbs corrected with tdfs
sbs<-sbs%>%
  mutate(mCr = mC + r_tdfs$Meand13c,
         sdCr = sdC + r_tdfs$SDd13c,
         mNr = mN + r_tdfs$Meand15n,
         sdNr = sdN + r_tdfs$SDd15n,
         mCp = mC + p_tdfs$Meand13c,
         sdCp = sdC + p_tdfs$SDd13c,
         mNp = mN + p_tdfs$Meand15n,
         sdNp = sdN + p_tdfs$SDd15n)

#Custom plot limits so we don't drop data
xlims <- c(-15,-20.5,-14,-20,-20,-16)
ylims<- c(8.5,16,10,17,7,13)

#Custom color palette from cooler (thanks Emily!!), and based on The Life Aquatic poster
fish_pal <- c('#B3D4E3','#729CAB','#306372','#8FA15A','#EEDE41','#B2A23B','#756534','#8C462F','#A22729')

#Plasma isospace plot
isospaceplot_p <- ggplot()+
  geom_errorbar(data = sbs, 
                mapping = aes(x = mCp, y = mNp,
                              ymin = mNp - 1.96*sdNp, 
                              ymax = mNp + 1.96*sdNp,
                              color = group), 
                width = 0.2)+
  geom_errorbarh(data = sbs, 
                 mapping = aes(x = mCp, y = mNp,
                               xmin = mCp - 1.96*sdCp,
                               xmax = mCp + 1.96*sdCp,
                               color = group),
                 height = 0.2) + 
  geom_point(data = sbs, aes(x = mCp, 
                             y = mNp,
                             col = group), 
             shape = 19, size = 5,
             alpha = 0.7, show.legend = F) + 
  scale_fill_manual(values = fish_pal)+
  geom_point(data=fishsource.new, aes(x = δ13C,
                                      y = δ15N,
                                      color = group), alpha = 0.5)+
  scale_color_manual(values = fish_pal, name = 'Prey Group')+
  guides(color=guide_legend(override.aes = list(alpha = 1)))+
  theme_minimal()+
  ggnewscale::new_scale_color()+
  geom_point(data=sidatP, aes(x = δ13C...V.PDB,
                              y = δ15N..air,
                              col=Species), shape = 17, size = 3, alpha = 0.8, show.legend = F)+
  scale_color_manual(values = c('#7b3294','#008837'))+
  xlim(values=c(min(xlims),max(xlims)))+
  ylim(values=c(min(ylims),max(ylims)))+
  ylab(NULL) +
  xlab(NULL) + 
  theme(axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank())+
  ggtitle('Plasma')

isospaceplot_rbc <- ggplot()+
  geom_errorbar(data = sbs, 
                mapping = aes(x = mCr, y = mNr,
                              ymin = mNr - 1.96*sdNr, 
                              ymax = mNr + 1.96*sdNr,
                              color = group), 
                width = 0.2, show.legend = F)+
  geom_errorbarh(data = sbs, 
                 mapping = aes(x = mCr, y = mNr,
                               xmin = mCr - 1.96*sdCr,
                               xmax = mCr + 1.96*sdCr,
                               color = group),
                 height = 0.2, show.legend = F) + 
  geom_point(data = sbs, aes(x = mCr, 
                             y = mNr,
                             color = group), 
             shape = 19, size = 5,
             alpha = 0.7, show.legend = F) + 
  geom_point(data=fishsource.new, aes(x = δ13C,
                                      y = δ15N,
                                      color = group), alpha = 0.5, show.legend = F)+
  scale_fill_manual(values = fish_pal)+
  scale_color_manual(values = fish_pal)+
  ggnewscale::new_scale_color()+
  theme_minimal()+
  geom_point(data=sidatRBC, aes(x = δ13C...V.PDB,
                         y = δ15N..air,
                         col=Species), shape = 17, size = 3, alpha = 0.8)+
  scale_color_manual(values = c('#7b3294','#008837'), name = 'Seal Species')+
  ylab(NULL) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  xlim(values=c(min(xlims),max(xlims)))+
  ylim(values=c(min(ylims),max(ylims)))+
  theme(axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
    panel.grid = element_blank())+
  ggtitle('RBC')

#Make it look sexy in patchwork
isopatch <- isospaceplot_p/isospaceplot_rbc +
  plot_layout(guides = 'collect')

wrap_elements(isopatch) +
  labs(tag = expression(paste(delta^{15}, "N (\u2030)"))) +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )

