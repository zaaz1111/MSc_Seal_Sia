fish.samples <- read.csv(here('fish.samples.csv'))#%>%
  dplyr::select(name:Scientific.Name)
fish.data <- read.csv(here('Processed Data/fish.samples.final.csv'))
mergedfish <- merge(fish.samples,fish.data)

#Ungrouped fish with Suess-corrected sandeel
fishsource.ungrouped <- read.csv(here('Processed data/Composite_source_data.csv'))
sandy <- fishsource.ungrouped[70:73,]
colnames(sandy) <- c('id','spp','d13c','nitrogen')
sandy <- mutate(sandy,
                year = c(1995,1996,2005,2005),
                region= c("Subpolar North Atlantic","Subpolar North Atlantic",
                          "Subpolar North Atlantic","Subpolar North Atlantic"))%>%
  subset(select=c(id,d13c,year,region))
suesssandy <- SuessR(sandy)
fishsource.ungrouped$δ13C[70:73] <- suesssandy$d13c.cor

write.csv(fishsource.ungrouped,file='ungrouped_source_data.csv')


####Silly stupid dumb fish re-cleaning hehe xD
Flatfish <- c('Lemon Sole', 'Long Rough Dab','Common Dab')
Scorpionfish <- c('Red Gurnard','Grey Gurnard')
Small_Benthics <- c('Blue Whiting','Norway Pout','Poor Cod')
Pelagics <- c('Cod', 'Argentine')
#Reload and pipe data to useable format
fishsource.new <- read.csv(here('Processed data/Composite_source_data.csv'))%>%
  mutate(group = case_when(
    fishsource.new$Spp %in% Flatfish ~ 'Flatfish',
    fishsource.new$Spp %in% Scorpionfish ~ 'Scorpionfish',
    fishsource.new$Spp %in% Small_Benthics ~ 'Small Benthics',
    fishsource.new$Spp %in% Pelagics ~ 'Pelagics',
    TRUE ~ fishsource.new$Spp
  ))

  #Suess correction for sandeel
sandy <- fishsource.new[70:nrow(fishsource.new),]
colnames(sandy) <- c('id','spp','d13c','nitrogen','group')
sandy <- mutate(sandy,
                year = c(1995,1996,2005,2005),
                region= c("Subpolar North Atlantic","Subpolar North Atlantic","Subpolar North Atlantic","Subpolar North Atlantic"))%>%
  subset(select=c(id,d13c,year,region))
suesssandy <- SuessR(sandy)
fishsource.new$δ13C[70:nrow(fishsource.new)] <- suesssandy$d13c.cor

write.csv(fishsource.new,file='grouped_fish_sources.csv')

subtract_grey_rbc <- c('Benthic Sandeel*','Flatfish')
fishsource.new.grey.rbc <- fishsource.new %>%
  mutate(group = case_when(
    fishsource.new$group %in% subtract_grey_rbc ~ NA,
    TRUE ~ fishsource.new$group
  ))%>%
  na.omit()
  
write.csv(fishsource.new.grey.rbc,file='fishsources_subracted_grey_rbc.csv')

subtract_harbor_rbc <- c('Benthic Sandeel*','Flatfish','Small Benthics')
fishsource.new.harbor.rbc <- fishsource.new %>%
  mutate(group = case_when(
    fishsource.new$group %in% subtract_harbor_rbc ~ NA,
    TRUE ~ fishsource.new$group
  ))%>%
  na.omit()

write.csv(fishsource.new.harbor.rbc,file='fishsources_subracted_harbor_rbc.csv')



subtract_grey_plasma <- c('Benthic Sandeel*')
fishsource.new.grey.plasma <- fishsource.new %>%
  mutate(group = case_when(
    fishsource.new$group %in% subtract_grey_plasma ~ NA,
    TRUE ~ fishsource.new$group
  ))%>%
  na.omit()

write.csv(fishsource.new.grey.plasma,file='fishsources_subracted_grey_plasma.csv')


subtract_harbor_rbc <- c('Benthic Sandeel*','Flatfish','Small Benthics')
fishsource.new.harbor.rbc <- fishsource.new %>%
  mutate(group = case_when(
    fishsource.new$group %in% subtract_harbor_rbc ~ NA,
    TRUE ~ fishsource.new$group
  ))%>%
  na.omit()

write.csv(fishsource.new.harbor.rbc,file='fishsources_subracted_harbor_rbc.csv')


#Harbor Plasma
subtract_harbor_plasma <- c('Benthic Sandeel*')
fishsource.new.grey.plasma <- fishsource.new %>%
  mutate(group = case_when(
    fishsource.new$group %in% subtract_harbor_plasma ~ NA,
    TRUE ~ fishsource.new$group
  ))%>%
  na.omit()

write.csv(fishsource.new.grey.plasma,file='fishsources_subracted_harbor_plasma.csv')




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
#Squid also fit in this group, weirdly

kw_pc <- kruskal.test(δ13C ~ Spp,
                      data=fishsource.new[fishsource.new$group=='Pelagics',])
kw_pc

kw_pn <- kruskal.test(δ15N ~ Spp,
                      data=fishsource.new[fishsource.new$group=='Pelagics',])
kw_pn
#Cod and Argentine are isotopically similar



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
sbs <- fishsource.new%>%
  group_by(Spp)%>%
  summarize(count = n(),      
            mC = mean(δ13C),
            sdC= sd(δ13C),
            mN = mean(δ15N), 
            sdN = sd(δ15N))
sbs$sdC[sbs$Spp=='Pelagic Sandeel*']<-0.5
sbs$sdN[sbs$Spp=='Pelagic Sandeel*']<-0.45

#Add in tdf data (adjusting the means and SDs one by one using mutate)
r_tdfs <- read.csv(here('Processed data/Seal_RBC_TDFS_ungrouped.csv'))
p_tdfs <- read.csv(here('Processed data/Seal_Plasma_TDFS_ungrouped.csv'))

p_tdfs_grouped <- read.csv(here('Processed data/Seal_Plasma_TDFs_grouped.csv'))
r_tdfs_grouped <- read.csv(here('Processed data/Seal_RBC_TDFs_grouped.csv'))
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
fish_pal <- c('#007289','#2B8AA1','#55A1B8','#A9D0E7','#BDD8AE','#D0DF74','#F6ED00','#EDB505','#E37C09','#CF0A11','#A12F09','#8A4105','#735300','purple','cyan','orange')

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

