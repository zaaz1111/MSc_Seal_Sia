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

#Do the same model but just for harbor seals
mix_harbs <- load_mix_data(filename = here('Processed Data/Consumer_harbs.csv'),
                     iso_names<-c('d13c', 'd15n'),
                     factors = NULL,
                     fac_random = NULL,
                     fac_nested = NULL,
                     cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_harbs)

discr_harbs <- load_discr_data(filename=here('less_dummy_tdfs_harbs.csv'), mix)

plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_harbs,source,discr_harbs)

model_filename <- "MixSIAR_model_msc_harbs_kinda_dummy_2.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_harbs, source)

# Run the JAGS model ("very long" took ~5 min)
jags.uninf_harbs_vl <- run_model(run="long",mix_harbs,source,discr_harbs,model_filename)


# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf_harbs_vl, mix_harbs, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_harbs_long",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_harbs_long",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_harbs_long", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_harbs_long", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_harbs_long", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                  FALSE, diag_save_ggmcmc = TRUE))

#Do the same model but now for grey seals
mix_greys <- load_mix_data(filename = here('Processed Data/Consumer_grey.csv'),
                           iso_names<-c('d13c', 'd15n'),
                           factors = NULL,
                           fac_random = NULL,
                           fac_nested = NULL,
                           cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_greys)

discr_greys <- load_discr_data(filename=here('less_dummy_tdfs_greys.csv'), mix_greys)

plot_data(filename="isospace_plot_greys", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_greys,source,discr_greys)

model_filename <- "MixSIAR_model_msc_greys_kinda_dummy_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_greys, source)

# Run the JAGS model ("very long" took ~5 min)
jags.uninf_greys_l <- run_model(run="long",mix_greys,source,discr_greys,model_filename)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf_greys_l, mix_greys, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_greys_l",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_greys_l",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_greys_l", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_greys_l", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_greys_l", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                  FALSE, diag_save_ggmcmc = TRUE))
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

plot_data(filename="isospace_plot_harbs_plasma", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_harbor_RBC,source,discr_harbor_RBC)

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
output_JAGS(jags.uninf_harbor_RBC_VL, mix_harbor_RBC, source,
            output_options = list(summary_save = TRUE, summary_name = "summary_statistics_harbor_RBC_vl",
                                  sup_post = FALSE, plot_post_save_pdf = TRUE, plot_post_name = "posterior_density_harbor_RBC_vl",
                                  sup_pairs = FALSE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_plot_harbor_RBC_vl", 
                                  sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot_harbor_RBC_vl", 
                                  gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_harbor_RBC_vl", 
                                  indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png =
                                    FALSE, diag_save_ggmcmc = TRUE))

###Now with Grey RBC's
mix_grey_RBC <- load_mix_data(filename = here('Processed Data/consumer_grey_RBC.csv'),
                                 iso_names<-c('d13c', 'd15n'),
                                 factors = NULL,
                                 fac_random = NULL,
                                 fac_nested = NULL,
                                 cont_effects = NULL)

source <- load_source_data(filename = here('Processed Data/Source_mixing_data.csv'),
                           source_factors = NULL,
                           conc_dep = F,
                           data_type = 'raw',
                           mix_grey_RBC)

discr_greys <- load_discr_data(filename=here('Processed Data/less_dummy_tdfs_greys.csv'), mix_grey_RBC)

plot_data(filename="isospace_plot_greys", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_grey_RBC,source,discr_greys)

model_filename <- "MixSIAR_model_msc_grey_plasma_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_grey_RBC, source)

# Run the JAGS model 
start_time = Sys.time()
jags.uninf_grey_RBC_vl <- run_model(run="very long",mix_grey_RBC,source,discr_greys,model_filename)
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
