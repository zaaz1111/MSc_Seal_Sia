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
