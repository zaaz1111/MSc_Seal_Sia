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

LR_fish <- c(fish1$name,fish2$name)
LR_fish<-LR_fish[nchar(LR_fish)<=4]

fish<-c(fish1$name[5:20],fish1$name[27:69],fish2$name[1:11])%>%
  data.frame()

colnames(fish)<-'LE_sample_name'
fish<-fish%>%
  mutate(LR_sample_name = case_when(
    substr(fish$LE_sample_name,1,4) %in% LR_fish ~ substr(LE_sample_name,1,4)))
    
le_cn<-{}
for(i in fish$LE_sample_name){
  for(z in bigfish$name){
    if(i==z){
      le_cn<-append(le_cn,bigfish$`C/N`[bigfish$name==z])
    }
  }
}

fish <- cbind(fish,le_cn)

lr_cn<-{}
for(i in na.omit(fish$LR_sample_name)){
  for(z in bigfish$name){
    if(i==z){
      lr_cn<-append(lr_cn,bigfish$`C/N`[bigfish$name==z])
    }
  }
}

index <- 1
lr_cn_ord <- {}
for(i in fish$LR_sample_name){
  if(is.na(i)==F){
    lr_cn_ord <- append(lr_cn_ord, lr_cn[index])
    index <- index + 1
  } else{
    lr_cn_ord <- append(lr_cn_ord, NA)
  }
}

fish <- cbind(fish,lr_cn_ord)

spp <- {}
for(i in fish$LE_sample_name){
  for(z in fish.samples$name){
    if(i==z){
      spp<-append(spp,fish.samples$Spp[fish.samples$name==z])
    }
  }
}
  
ish<-fish[-16,]
fish <- cbind(fish,spp)

