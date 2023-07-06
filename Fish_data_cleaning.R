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

ggplot(megafish)+
  geom_boxplot(aes(x=Spp,y=CN_ratio,fill=Lipid.Extracted.))+
  ggtitle('Lipid Extracted Samples')+
  xlab('C/N Ratio')+
  ylab('Species')+
  guides(color=guide_legend(title='Lipid Extracted?'),labels=c('No','Yes'))

                            