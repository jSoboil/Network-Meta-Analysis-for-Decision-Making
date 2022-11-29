# Misc Setup --------------------------------------------------------------
# required packages for analysis
pkgs <- c("R2jags", "rstan", "rstanarm", "bayesplot", "dplyr", 
          "tidyr", "igraph")
# install packages if not installed
installed_packages <- pkgs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
 install.packages(pkgs[!installed_packages])
}
# load required packages
invisible(lapply(pkgs, library, character.only = TRUE))
# set number of parallel cores
options(mc.cores = parallel::detectCores())

# Model Fit,  Model Comparison,  and Outlier Detection --------------------
