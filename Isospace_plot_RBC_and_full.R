faffr <- subset(sbs, select = c(Spp, mCr, sdCr, mNr, sdNr))%>%
  mutate(cmax = mCr + 1.96*sdCr,
         cmin = mCr - 1.96*sdCr,
         nmax = mNr + 1.96*sdNr,
         nmin = mNr - 1.96*sdNr)

x <- c(faffr$cmax,faffr$cmin,faffr$mCr,faffr$mCr)
y <- c(faffr$mNr,faffr$mNr,faffr$nmax,faffr$nmin)
hullr <- chull(x,y)

hullrx <- {}
hullry<- {}
for(i in hullr){
  hullrx <- append(hullrx, x[i])
  hullry <- append(hullry, y[i])
}
hullframer <- data.frame(x=hullrx,y=hullry)
hullframe<- rbind(hullframer,hullframep)

ggplot()+
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
                 height = 0.2, show.legend = T) +
  geom_point(data = sbs, aes(x = mCr,
                             y = mNr,
                             color = group),
             shape = 19, size = 5,
             alpha = 0.7, show.legend = F) +
  scale_fill_manual(values = fish_pal)+
  scale_color_manual(values = fish_pal, name = 'Prey Group')+
  ggnewscale::new_scale_color()+
  theme_minimal()+
  geom_point(data=sidatRBC, aes(x = δ13C...V.PDB,
                                y = δ15N..air, fill=Species),
                                size = 3, alpha = 0.5,pch = 21, col = 'black')+
  scale_color_manual(values = c('#7b3294','#008837'), name = 'Seal Species')+
  ylab(NULL) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  geom_polygon(data=hullframer, aes(x=x,y=y),fill= 'transparent', color = 'black', lty='dotdash',lwd = 1)+
  xlim(values=c(min(hullframe$x),max(hullframe$x)))+
  ylim(values=c(min(hullframe$y),max(sidatRBC$X.N)+0.5))+
  theme(axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank(), aspect.ratio = 1)+
  ggtitle('RBC')+
  facet_wrap(~Species)
##Copy this ^


##To paste here
isoRBC <- ggplot()+
  geom_errorbar(data = sbs,
                mapping = aes(x = mCr, y = mNr,
                              ymin = mNr - 1.96*sdNr,
                              ymax = mNr + 1.96*sdNr,
                              color = Spp),
                width = 0.2, show.legend = F)+
  geom_errorbarh(data = sbs,
                 mapping = aes(x = mCr, y = mNr,
                               xmin = mCr - 1.96*sdCr,
                               xmax = mCr + 1.96*sdCr,
                               color = Spp),
                 height = 0.2, show.legend = F) +
  geom_point(data = sbs, aes(x = mCr,
                             y = mNr,
                             color = Spp),
             size = 5,
             alpha = 0.7)+#, show.legend = F) +
  scale_fill_manual(values = fish_pal)+
  scale_color_manual(values = fish_pal, name = 'Prey Group')+
  ggnewscale::new_scale_fill()+
  theme_minimal()+
  geom_point(data=sidatRBC,   aes(x = δ13C...V.PDB,
                                  y = δ15N..air, fill=Species),
             size = 3, alpha = 0.5,pch = 21, col = 'black')+
  scale_fill_manual(values = c('#7b3294','#008837'), name = 'Seal Species')+
  ylab(NULL) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  geom_polygon(data=hullframer, aes(x=x,y=y),fill= 'transparent', color = 'black', lty='dotdash',lwd = 1)+
  xlim(values=c(min(hullframe$x),max(hullframe$x)))+
  ylim(values=c(min(hullframe$y),max(hullframe$y)))+
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  ggtitle('RBC')+
  theme(axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank(), aspect.ratio = 1,
        plot.title=element_text(hjust=0.5))+
  facet_wrap(~Species, scales = 'free')
  
isoSpace <- isoPlasma / isoRBC 
isoSpace + plot_layout(guides = "collect")
dev.off()

ggsave(filename = 'big.isospace.jpeg')
