library(tidyverse)
###THIS PLOT NEEDS COLOR SCHEMES

ggplot(data = sidat, aes(x = δ13C...V.PDB, y = δ15N..air)) + 
  geom_point(aes(color = Species), size = 2) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  stat_ellipse(aes(color = Species, fill = Species), 
               alpha = 0.25, 
               level = 0.40,
               type = "norm",
               geom = "polygon")+
  stat_ellipse(aes(color = Species),
               fill = NA,
               alpha = 0.25, 
               level = 0.95,
               type = "norm",
               geom = "polygon",
               lty=2) + 
  scale_fill_manual(values=c('#008837','#7b3294'),labels = c('Grey','Harbour'))+
  scale_color_manual(values=c('#008837','#7b3294'),labels = c('Grey','Harbour'))+
  scale_shape_manual(values=c(19,17))+
  guides(shape = 'none')+
  ylim(values=c(11,16))+
  xlim(values=c(-20.5,-16))+
  facet_wrap(~Type, scales = 'free')+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 15),
        axis.line = element_line(colour = 'black', linewidth = 1, lineend = 'square'),
        panel.grid = element_blank(),aspect.ratio = 1)


ggsave(filename = here('Siber Ellipse Plot.jpeg'))
