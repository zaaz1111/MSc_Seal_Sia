harbor_priors <- c(
  Argentine_g <- 31.4,
  Bluemouth_g <- 0,
  Cod_g <- 23.9,
  Flatfish_g <- 1.3,
  Scorpionfish_g <- 0,
  Squid <- 0.9,
  Trisopterus_g <- 8.6)

harbor_alpha_priors <- harbor_priors*length(harbor_priors)/sum(harbor_priors)
harbor_alpha_priors[which(harbor_alpha_priors==0)] <- 0.01

# Plot your informative prior
plot_prior(alpha.prior=harbor_alpha_priors,
           source=source,
           plot_save_pdf=TRUE,
           plot_save_png=FALSE,
           filename="prior_plot_harbor_inf")
