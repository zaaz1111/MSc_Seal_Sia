library(librarian)
shelf(MixSIAR,SIBER,here,tidyverse,ggplot2,rjags,ellipse,readxl)

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



mix_harbs <- load_mix_data(filename = here('Processed Data/Consumer_harbs.csv'),
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

discr_harbs <- load_discr_data(filename=here('less_dummy_tdfs_harbs.csv'), mix)

plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix_harbs,source,discr_harbs)

model_filename <- "MixSIAR_model_msc_harbs_kinda_dummy_1.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("very long" took ~5 min)
jags.uninf_harbs <- run_model(run="test",mix_harbs,source,discr_harbs,model_filename)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf, mix, source)
