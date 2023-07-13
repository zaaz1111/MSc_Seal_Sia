###Processing fish sampling data
fs<-matrix(nrow=70,ncol=7)
colnames(fs)<-c('Sample Name','Lipid Extracted?','Spp','Tow','Replicate','Plate','Well')
fs<-data.frame(fs)
for(i in 1:nrow(fs)){
  fs$Sample.Name[i]<-paste('FS',i,'LE',sep='')
  fs$Lipid.Extracted.<-TRUE
}

write.csv(fs,file='fish.samples.csv')
fs<-read.csv("C:/Users/zaahi/Documents/fish.samples.csv")%>%
  mutate(fs,Sample.Name=case_when(
    fs$Lipid.Extracted. == T ~ fs$Sample.Name,
    fs$Lipid.Extracted. == F ~ str_sub(fs$Sample.Name, 1, -3)
  ))

fish.samples <- read.csv(here('fish.samples.csv'))
fish.data <- read.csv(here('Processed Data/LE_Fish_CN_Collated.csv'))
fish.merge<-merge(fish.data,fish.samples)%>%
  subset(select=c(1,2,3,9,11))
write.csv(fish.merge, file='Source_mixing_data.csv')

###Attempt to ANOVA the fish groups
fsle<-subset(megafish, subset=megafish$Lipid.Extracted.==T,select=c('name','Spp',
  'δ15N','δ13C'))

daov <- aov(δ15N ~ δ13C, data = fsle)



d13caov <- aov(δ13C~Spp,data=fsle)
hsd13<-data.frame(TukeyHSD(d13caov)$Spp)
hsd13[hsd13$p.adj <= 0.05,]

d15naov <- aov(δ15N~Spp,data=fsle)
hsd15<-data.frame(TukeyHSD(d15naov)$Spp)
hsd15<-hsd15[hsd15$p.adj <= 0.05,]

iso2way <- aov(δ13C ~  Spp + δ15N, data = fsle)
TukeyHSD(iso2way)


gadid <- c()
trisopterus <- c('Norway Pout', 'Poor Cod','Blue Whiting') #Not sure about blue whiting yet
Flatfish <- c('Common Dab','Long Rough Dab', 'Plaice', 'Lemon Sole')
scorpions <- c('Grey Gurnard', 'Red Gurnard')
sandeel <- c()

megafish <- megafish%>%
  mutate(group = case_when(
    megafish$Spp %in% trisopterus ~ 'trisopterus',
    megafish$Spp %in% Flatfish ~ 'Flatfish',
    megafish$Spp %in% scorpions ~'scorpionfish',
    megafish$Spp %in% c(trisopterus,Flatfish,scorpions)==F ~megafish$Spp
  ))
write.csv(megafish[megafish$Lipid.Extracted.==T,],file='Source_mixing_data_grouped.csv')



megafish2<-megafish
megafish2$δ13C <- megafish2$δ13C/mean(megafish2$δ13C)
megafish2$δ15N <- megafish2$δ15N/mean(megafish2$δ15N)


megafish2[megafish2$Lipid.Extracted.==T,]%>%
  pivot_longer(
    cols = c("δ13C", "δ15N"),
    names_to = "Iso",
    values_to = "Val"
  )%>%
#plot it idk
  ggplot()+
    geom_boxplot(aes(x=Spp,y=Val,fill=Iso))+#,color=Lipid.Extracted.))+
    ggtitle('Lipid Extracted Samples')+
    xlab('Species')+
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+
    guides(color=guide_legend(title='Lipid Extracted?'),labels=c('No','Yes'))+
    theme_minimal()





###Compare LE and LR fish tissue
#Load in any dfs with fish data
#separate them by LE and LR
#rename the columns so they have LE/LR naming conventions
#Pull the CN ratios and add to the dfs
#Plot clustered boxplots
fish1 <- read_xlsx(here('Raw Data/Fish and Seal Plate Collated Raw (Plate 2).xlsx'))
fish2 <- read_csv(here('Raw Data/Cleaned_plate_3.csv'))
bigfish <- rbind(fish1,fish2)

megafish<-merge(fish.samples, bigfish)
colnames(megafish)[15]<-'CN_ratio'

#Look at LE vs NLE
ggplot(megafish)+
  geom_boxplot(aes(x=Spp,y=CN_ratio,fill=Lipid.Extracted.))+
  ggtitle('Lipid Extracted Samples')+
  xlab('C/N Ratio')+
  ylab('Species')+
  guides(color=guide_legend(title='Lipid Extracted?'),labels=c('No','Yes'))


ggplot(megafish)+
  geom_boxplot(aes(x=Spp,y=CN_ratio,fill=Lipid.Extracted.))+
  ggtitle('Lipid Extracted Samples')+
  xlab('C/N Ratio')+
  ylab('Species')+
  guides(color=guide_legend(title='Lipid Extracted?'))+
  theme_minimal()

###Fish Isospace Plotting
sbs <- megafish%>%
  group_by(Spp)%>%
  summarize(count = n(),
            mC = mean(δ13C),
            sdC= sd(δ13C),
            mN = mean(δ15N), 
            sdN = sd(δ15N))

#
ggplot(data = megafish, 
                     aes(x = δ13C, 
                         y = δ15N)) + 
  geom_point(aes(color = Spp), size = 5) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_viridis_d()

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
