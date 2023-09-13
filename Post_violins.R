p.global <- jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global

##
post.mix_spp1 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,1])
post.mix_spp2 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,2])
post.mix_spp3 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,3])
post.mix_spp4 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,4])
post.mix_spp5 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,5])
post.mix_spp6 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,6])
post.mix_spp7 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,7])
post.mix_spp8 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,8])
post.mix_spp9 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,9])
post.mix_spp10 <- as.data.frame(jags.uninf_grey_RBC_extreme$BUGSoutput$sims.list$p.global[,10])

#Name them!!
colnames(post.mix_spp1) <- source$source_names[1]
colnames(post.mix_spp2) <- source$source_names[2]
colnames(post.mix_spp3) <- source$source_names[3]
colnames(post.mix_spp4) <- source$source_names[4]
colnames(post.mix_spp5) <- source$source_names[5]
colnames(post.mix_spp6) <- source$source_names[6]
colnames(post.mix_spp7) <- source$source_names[7]
colnames(post.mix_spp8) <- source$source_names[8]
colnames(post.mix_spp9) <- source$source_names[9]
colnames(post.mix_spp10) <- source$source_names[10]
###combine all file into one
post_dists_grey_rbc <- cbind(post.mix_spp1,
                             post.mix_spp2,
                             post.mix_spp3,
                             post.mix_spp4,
                             post.mix_spp5,
                             post.mix_spp6,
                             post.mix_spp7,
                             post.mix_spp8,
                             post.mix_spp9,
                             post.mix_spp10)%>%
  gather(source,value,1:10)


ggplot(post_dists_grey_rbc)+
  geom_violin(aes(y=value,x=group, fill=source))+
  theme_minimal()

###Posty, a function that takes a completed JAGS model and source file and gives the posts
#we do it live lmao
posty <- function(x,y){
  p.global <- x$BUGSoutput$sims.list$p.global
  post.mix_spp1 <- as.data.frame(x$BUGSoutput$sims.list$p.global[,1])
  post.mix_spp2 <- as.data.frame(x$BUGSoutput$sims.list$p.global[,2])
  post.mix_spp3 <- as.data.frame(x$BUGSoutput$sims.list$p.global[,3])
  post.mix_spp4 <- as.data.frame(x$BUGSoutput$sims.list$p.global[,4])
  # post.mix_spp5 <- as.data.frame(x$BUGSoutput$sims.list$p.global[,5])
  # post.mix_spp6 <- as.data.frame(x$BUGSoutput$sims.list$p.global[,6])
  
  colnames(post.mix_spp1) <- y$source_names[1]
  colnames(post.mix_spp2) <- y$source_names[2]
  colnames(post.mix_spp3) <- y$source_names[3]
  colnames(post.mix_spp4) <- y$source_names[4]
  # colnames(post.mix_spp5) <- y$source_names[5]
  # colnames(post.mix_spp6) <- y$source_names[6]
  
  return(cbind(post.mix_spp1,
               post.mix_spp2,
               post.mix_spp3,
               post.mix_spp4)%>%#,
               # post.mix_spp5,
               # post.mix_spp6)%>%
           gather(y,value,1:4))
}

posts_grey_plasma<-posty(jags.inf_grey_plasma_extra_s2,source_s2)
posts_grey_plasma$group <- 'Grey Plasma'

posts_harbor_plasma <- posty(jags.inf_harbor_plasma_extra_ungrouped_s3, source_ungrouped_s3)
posts_harbor_plasma$group <- 'Harbour Plasma'

posts_grey_rbc_autumn <- posty(jags.inf_grey_RBC_extreme_s1_autumn, source_s1_autumn)
posts_grey_rbc_autumn$group <- 'Grey RBC Autumn/Winter'

posts_grey_rbc <- posty(jags.inf_grey_RBC_extreme_s3, source_s3)
posts_grey_rbc$group <- 'Grey RBC Spring/Summer'

posts_harb_rbc <- posty(jags.inf_harb_RBC_super_ungrouped_s2, source_ungrouped_s2)
posts_harb_rbc$group <- 'Harbour RBC Spring/Summer'

posts_harb_rbc_autumn <- posty(jags.inf_harb_RBC_super_ungrouped_autumn_s3, source_ungrouped_s3_a)
posts_harb_rbc_autumn$group <- 'Harbour RBC Autumn/Winter'


plasmapostframe <- rbind(posts_grey_plasma, posts_harbor_plasma)

plas_greyrbcframe <- rbind(posts_grey_plasma,posts_harbor_plasma,posts_grey_rbc_autumn,posts_grey_rbc)

finalframe <- rbind(plas_greyrbcframe, posts_harb_rbc, posts_harb_rbc_autumn)

finalframe <- mutate(finalframe, Type = case_when(substr(group, nchar(group), nchar(group)) == 'a' ~ 'Plasma',
                                       substr(group, nchar(group), nchar(group)) == 'r' ~ 'RBC'),
                     Season = case_when(substr(group, nchar(group)-2,nchar(group)-2) == 'm' ~ 'Spring/Summer',
                                        substr(group, nchar(group)-2,nchar(group)-2) == 't' ~ 'Autumn/Winter',
                                        substr(group, nchar(group), nchar(group)) == 'a' ~ 'Plasma'),
                     Species = case_when(substr(group, 1, 1) == 'G'~'Grey',
                           substr(group, 1, 1) == 'H'~'Harbour'))


sznlabs <- as_labeller(c('Autumn/Winter' = 'RBC Autumn/Winter',
                'Spring/Summer' = 'RBC Spring/Summer',
                'Plasma' = 'Plasma'))


dodge <- position_dodge(width = 1)
ggplot(finalframe)+
  geom_violin(aes(y=value,x=y, fill=y),position = dodge)+
  scale_fill_manual(values=c('#007289','#2B8AA1','#55A1B8','#A9D0E7','#BDD8AE','#F6ED00','#A12F09','cyan'), name='Diet Components')+
  #geom_boxplot(aes(y=value, x=group, fill=y),position = dodge, width = 0.25)+
  facet_grid(Season ~ Species, scales = 'free', space = 'free',
             labeller = labeller(Season = sznlabs))+
  xlab('Prey Spp')+
  ylab('Scaled Posterior Density')+
  theme_minimal()+
  theme(axis.text.y = element_text(size=8, angle = 90),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size=12),
        strip.text = element_text(size = 10))
ggsave(filename = 'posts.plot.jpeg')


###Steps
###Load in each post set
###Make each post set a dataframe
###Add a group column to each
### rbind()
### ggplot()
###geom_violin, x = group, y = value, fill = source