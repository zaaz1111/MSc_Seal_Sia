library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,readxl,beepr)

mix <- load_mix_data(filename = here('Processed Data/Consumer_mixing_data.csv'),
                     iso_names<-c('δ13C...V.PDB', 'δ15N..air'),
                     factors = NULL,
                     fac_random = NULL,
                     fac_nested = NULL,
                     cont_effects = NULL)


source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix)


source<-read.csv(here('Processed Data/Source_mixing_data.csv'))

ggplot(source)+
  geom_point(aes(d13c,d15n,color=Spp))+
  theme_minimal()

#Make discrimination data with published values from the literature
discr <- load_discr_data(filename=here('Dummy_tdfs.csv'), mix)

# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model_msc_seals_dummy_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("very long" took ~5 min)
jags.uninf <- run_model(run="test",mix,source,discr,model_filename)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf, mix, source)


#Further Split for plasma and rbc's???
mix_grey_plasma <- load_mix_data(filename = here('Processed Data/consumer_grey_plasma.csv'),
                           iso_names<-c('d13c', 'd15n'),
                           factors = NULL,
                           fac_random = NULL,
                           fac_nested = NULL,
                           cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_grey_plasma)

discr_greys <- load_discr_data(filename=here('Processed Data/less_dummy_tdfs_greys.csv'), mix_grey_plasma)

plot_data(filename="isospace_plot_greys", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_grey_plasma,source,discr_greys)

model_filename <- "MixSIAR_model_msc_grey_plasma_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_grey_plasma, source)

# Run the JAGS model 
start_time = Sys.time()
jags.uninf_grey_plasma_vl <- run_model(run="very long",mix_grey_plasma,source,discr_greys,model_filename)
end_time = Sys.time()
end_time - start_time

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf_grey_plasma_l, mix_grey_plasma, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_grey_plasma_l",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_grey_plasma_l",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_grey_plasma_l", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_greys_l", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_grey_plasma_l", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                    FALSE, diag_save_ggmcmc = TRUE))


#Further Split for plasma and rbc's???
mix_harbor_plasma <- load_mix_data(filename = here('Processed Data/consumer_harbs_plasma.csv'),
                                 iso_names<-c('d13c', 'd15n'),
                                 factors = NULL,
                                 fac_random = NULL,
                                 fac_nested = NULL,
                                 cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_harbor_plasma)

discr_harbor_plasma <- load_discr_data(filename=here('Processed Data/Less_dummy_tdfs_harbs.csv'), mix_harbor_plasma)

plot_data(filename="isospace_plot_harbs_plasma", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_harbor_plasma,source,discr_harbor_plasma)

model_filename <- "MixSIAR_model_msc_harbor_plasma_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_harbor_plasma, source)

# Run the JAGS model 
start_time = Sys.time()
jags.uninf_harbor_plasma_ <- run_model(run="very long",mix_harbor_plasma,source,discr_harbor_plasma,model_filename)
end_time = Sys.time()
end_time - start_time

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf_harbor_plasma_, mix_harbor_plasma, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_harbor_plasma_vl",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_harbor_plasma_vl",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_harbor_plasma_vl", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_greys_l", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_harbor_plasma_vl", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                    FALSE, diag_save_ggmcmc = TRUE))

###Harbor RBC's
mix_harbor_RBC <- load_mix_data(filename = here('Processed Data/consumer_harbs_RBC.csv'),
                                   iso_names<-c('d13c', 'd15n'),
                                   factors = NULL,
                                   fac_random = NULL,
                                   fac_nested = NULL,
                                   cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_harbor_RBC)

discr_harbor_RBC <- load_discr_data(filename=here('Processed Data/Less_dummy_tdfs_harbs.csv'), mix_harbor_RBC)

plot_data(filename="isospace_plot_harbs_plasma", plot_save_pdf=TRUE, plot_save_png=FALSE, 
          mix_harbor_RBC,source,discr_harbor_RBC)

model_filename <- "MixSIAR_model_msc_harbor_RBC_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_harbor_RBC, source)

# Run the JAGS model 
start_time = Sys.time()
jags.uninf_harbor_RBC_VL<- run_model(run="very long",mix_harbor_RBC,source,discr_harbor_RBC,model_filename)
end_time = Sys.time()
end_time - start_time

# Process diagnostics, summary stats, and posterior plots

output_JAGS(jags.uninf_harbor_RBC_extreme, mix_harbor_RBC, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_harbor_RBC_vl",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_harbor_RBC_vl",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_harbor_RBC_vl", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_harbor_RBC_vl", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_harbor_RBC_vl", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                    FALSE, diag_save_ggmcmc = TRUE))

###Using this one!!!
mix_grey_RBC <- load_mix_data(filename = here('Processed Data/consumer_grey_RBC.csv'),
                                 iso_names<-c('d13c', 'd15n'),
                                 factors = NULL,
                                 fac_random = NULL,
                                 fac_nested = NULL,
                                 cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data_grouped.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_grey_RBC)

discr_greys <- load_discr_data(filename=here('Processed Data/less_dummy_tdfs_greys.csv'), mix_grey_RBC)

plot_data(filename="isospace_plot_greys_RBC_grouped_samples", plot_save_pdf=TRUE, 
          plot_save_png=FALSE, mix=mix_grey_RBC, source=source, discr=discr_greys)

grey_rbc_alpha <- c(0,0,35.3,6.5,33.5,0.6,4.5)

grey_rbc_alpha <- grey_rbc_alpha*length(grey_rbc_alpha)/sum(grey_rbc_alpha)

grey_rbc_alpha[which(grey_rbc_alpha==0)] <- 0.01

# Plot your informative prior
plot_prior(alpha.prior=grey_rbc_alpha,
           source=source,
           plot_save_pdf=TRUE,
           plot_save_png=FALSE,
           filename="prior_plot_kw_inf")




model_filename <- "MixSIAR_model_msc_grey_RBC_2_grouped_sources_inf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_grey_RBC, source)

# Run the JAGS model 
start_time = Sys.time()
jags.uninf_grey_RBC_extreme <- run_model(run="extreme",mix_grey_RBC,source,discr_greys,
                                    model_filename, alpha.prior = grey_rbc_alpha)
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
