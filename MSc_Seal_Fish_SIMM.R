library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,readxl,beepr,colorspace)

#Load in the Consumer mixing data
mix_grey_RBC <- load_mix_data(filename = here('Processed Data/consumer_grey_RBC.csv'),
                                 iso_names<-c('d13c', 'd15n'),
                                 factors = NULL,
                                 fac_random = NULL,
                                 fac_nested = NULL,
                                 cont_effects = NULL)


#Load in the source mixing data
source <- load_source_data(filename = here('idfkanymore..csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_grey_RBC)
#Load in the TDF data
discr_greys <- load_discr_data(filename=here('Processed Data/Semi_smart_Seal_RBC_TDFS.csv'), mix_grey_RBC)

#Create an isospace plot
grp <- plot_data(filename="isospace_plot_greys_RBC_grouped_samples", plot_save_pdf=F, 
          plot_save_png=FALSE, mix=mix_grey_RBC, source=source, 
          discr=discr_greys, return_obj = T)
grp+  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=20), axis.text=element_text(size=14),
        legend.text=element_text(size=12), legend.title=element_text(size=14))
#Load in priors (Maybe a note for which prior is which value?)
#######MixSIAR Reads priors alphabetically by source data!
###Grey Seal Priors
grey_priors <- c(
  Argentine_g <- 0,
  Bluemouth_g <- 33.6,
  Cod_g <- 35.3,
  Flatfish_g <- 6.5,
  Sandeel_g <- 18.8,
  Scorpionfish_g <- 33.6,
  Squid <- 0.9,
  Trisopterus_g <- 4.5
)

#Process the priors
grey_alpha_priors <- grey_priors*length(grey_priors)/sum(grey_priors)
grey_alpha_priors[which(grey_alpha_priors==0)] <- 0.01

# Plot your informative prior
plot_prior(alpha.prior=grey_alpha_priors,
           source=source,
           plot_save_pdf=TRUE,
           plot_save_png=FALSE,
           filename="prior_plot_grey_inf")

model_filename <- "MixSIAR_model_msc_grey_RBC_2_grouped_sources_inf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_grey_RBC, source)

# Run the JAGS model 
start_time = Sys.time()
jags.uninf_grey_RBC_extreme <- run_model(run="test",mix_grey_RBC,source,discr_greys,
                                    model_filename, alpha.prior = grey_alpha_priors)
end_time = Sys.time()
end_time - start_time

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf_grey_RBC_vl, mix_grey_RBC, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_grey_RBC_vl",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_grey_RBC_vl",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_grey_RBC_vl", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_grey_RBC_vl", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_grey_grey_RBC_vl", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                    FALSE, diag_save_ggmcmc = TRUE))
