library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,coda,readxl,patchwork,sf,terra,tidyterra,ggspatial,colorspace)

ukv <-vect('C:/Users/zaahi/Documents/england_is_my_city.shp')

harbor_usage <- rast(here('Harbour_seal_at-sea_usage_geotiffs/PvSeaUsage.tif'))%>%
  project('epsg:32630')%>%
  crop(ukv)
harbor_usage[harbor_usage<1]<-NA
harbor_usage[harbor_usage>100]<-100


grey_usage <- rast(here('Grey_seal_at_sea_usage/HgSeaUsage.tif'))%>%
  project('epsg:32630')%>%
  crop(ukv)
grey_usage[grey_usage<1]<-NA
grey_usage[grey_usage>100]<-100

uk <- st_read('C:/Users/zaahi/Documents/england_is_my_city.shp')%>%
  subset(select = geometry)%>%
  st_transform(crs = 'epsg:32630')

overlap <- raster::intersect(harbor_usage,grey_usage)
overlap[overlap<=0] <- NA

huplot <- ggplot()+
  geom_sf(data=uk, col='black',fill = 'black')+
  geom_spatraster(data=harbor_usage)+
  scale_fill_whitebox_b(pal='bl_yl_rd',guide=guide_colorbar(title='Number of Seals'))+
  ggtitle('Harbour Seal Sea Usage')+
  theme_minimal()

guplot <- ggplot()+
  geom_sf(data=uk, col='black',fill = 'black')+
  geom_spatraster(data=grey_usage)+
  scale_fill_whitebox_b(pal='bl_yl_rd',guide=guide_colorbar(title = 'Number of Seals'))+
  ggtitle('Grey Seal Sea Usage')+
  theme_minimal()

oplot <- ggplot()+
  geom_sf(data=uk, col='black',fill = 'black')+
  geom_spatraster(data=overlap)+
  scale_fill_hypso_d(pal='nordisk-familjebok_hypso',guide=NULL)+
  ggtitle('Harbour/Grey Seal Overlap')+
  theme_minimal()

huplot +oplot+ guplot
