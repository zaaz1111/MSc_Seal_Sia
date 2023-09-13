faffp <- subset(sbs, select = c(Spp, mCp, sdCp, mNp, sdNp))%>%
  mutate(cmax = mCp + 1.96*sdCp,
         cmin = mCp - 1.96*sdCp,
         nmax = mNp + 1.96*sdNp,
         nmin = mNp - 1.96*sdNp)

x <- c(faffp$cmax,faffp$cmin,faffp$mCp,faffp$mCp)
y <- c(faffp$mNp,faffp$mNp,faffp$nmax,faffp$nmin)
hullp <- chull(x,y)

hullpx <- {}
hullpy<- {}
for(i in hullp){
  hullpx <- append(hullpx, x[i])
  hullpy <- append(hullpy, y[i])
}
hullframep <- data.frame(x=hullpx,y=hullpy)
##Need to fix the prey groups in sbs but otherwise living fecking large
ggplot()+
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
             alpha = 0.7, show.legend = T) + 
  scale_fill_manual(values = fish_pal)+
  scale_color_manual(values = fish_pal, name = 'Prey Group')+
  guides(color=guide_legend(override.aes = list(alpha = 1)))+
  theme_minimal()+
  ggnewscale::new_scale_color()+
  geom_point(data=sidatP, aes(x = δ13C...V.PDB,
                              y = δ15N..air,
                              col=Species), shape = 17, size = 3, alpha = 0.5)+
  scale_color_manual(values = c('#7b3294','#008837'), name = 'Seal Species')+
  geom_polygon(data=hullframep,aes(x=x,y=y),fill= 'transparent', color = 'black', lty='dotdash',lwd=1)+
  xlim(values=c(min(hullframe$x),max(hullframe$x)))+
  ylim(values=c(min(hullframe$y),max(hullframe$y)))+
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank(), aspect.ratio = 1)+
  ggtitle('Plasma')+
  facet_wrap(~Species)
##Copy this ^

##To paste here
isoPlasma <- ggplot()+
  geom_errorbar(data = sbs, 
                mapping = aes(x = mCp, y = mNp,
                              ymin = mNp - 1.96*sdNp, 
                              ymax = mNp + 1.96*sdNp,
                              color = Spp), 
                width = 0.2)+
  geom_errorbarh(data = sbs, 
                 mapping = aes(x = mCp, y = mNp,
                               xmin = mCp - 1.96*sdCp,
                               xmax = mCp + 1.96*sdCp,
                               color = Spp),
                 height = 0.2) + 
  geom_point(data = sbs, aes(x = mCp, 
                             y = mNp,
                             col = Spp), 
             size = 5,
             alpha = 0.7, show.legend = T) + 
  scale_fill_manual(values = fish_pal)+
  scale_color_manual(values = fish_pal, name = 'Prey Group')+
  guides(color=guide_legend(override.aes = list(alpha = 1)))+
  theme_minimal()+
  ggnewscale::new_scale_fill()+
  theme_minimal()+
  geom_point(data=sidatRBC,   aes(x = δ13C...V.PDB,
                                  y = δ15N..air, fill=Species),
             size = 3, alpha = 0.5,pch = 21, col = 'black')+
  scale_fill_manual(values = c('#7b3294','#008837'), name = 'Seal Species')+
  guides(color='none')+
  geom_polygon(data=hullframep,aes(x=x,y=y),fill= 'transparent', color = 'black', lty='dotdash',lwd=1)+
  xlim(values=c(min(hullframe$x),max(hullframe$x)))+
  ylim(values=c(min(hullframe$y),max(hullframe$y)))+
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  ggtitle('Plasma')+
  theme(axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank(), aspect.ratio = 1,
        plot.title=element_text(hjust=0.5))+
  facet_wrap(~Species, scales = 'free')

