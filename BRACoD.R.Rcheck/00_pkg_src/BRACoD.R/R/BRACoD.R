
.onLoad <- function(libname, pkgname) {
  packageStartupMessage("Welcome to BRACoD. You need the python version of BRACoD for this to work. Use the function install_bracod() to install it.")
  BRACoD <<- reticulate::import("BRACoD", delay_load = TRUE)

}


#' Install BRACoD in python
#'
#' Uses pip to install the latest BRACoD release in python. You might need
#' to specify a python environment with either reticulate::use_virtualenv or
#' reticulate::use_condaenv. 
#' @export
install_bracod <- function(method = "auto", conda = "auto") {
  reticulate::py_install("BRACoD", method = method, conda = conda, pip=TRUE)
}


#' Simulate microbiome counts
#' 
#' Each bacteria's absolute abundance is simulated from a lognormal distribution.
#' Then, convert each sample to relative abundance, and simulate sequencing counts
#' using a multinomial distribution, based on the desired number of reads and the 
#' simulated relative abundances. This also simulates an environmntal variable that
#' is produced by some of the bacteria.
#' @param df A dataframe of OTU counts that is a model for data simulation. Samples are rows and bacteria are columns.
#' @param n_contributors the number of bacteria that are to contribute to your environmental variable.
#' @param coeff_contributor the average of the distribution used to simulate the contribution coefficient.
#' @param min_ab_contributor The minimum log relative abundance, averaged across samples, to include a bacteria
#' @param sd_Y the standard deviation of the simulated environmental variable
#' @param n_reads the number of reads to be simulated per sample
#' @param var_contributor If you use a uniform distribution, this is the range of the distribution, with a normal distribution it is the variance used to simulate the contribution coefficient.
#' @param use_uniform use a uniform distribution to simulate the contribution coefficient. Alternative is the normal distribution.
#' @param n_samples_use number of microbiome samples to simulate. If NULL, uses the same number of samples as in your dataframe
#' @param corr_value the bug-bug correlation value you want to include in the simulation
#' @param return_absolute returns the abosulte abundance values instead of the simulated microbiome counts
#' @param seed random seed for reproducibility
#' @export
simulate_microbiome_counts <- function(df, n_contributors = 20, coeff_contributor = 0.0, min_ab_contributor = -9, sd_Y = 1.0, n_reads = 100000, var_contributor = 5.0, use_uniform = TRUE, n_samples_use = NULL, corr_value = NULL,  return_absolute = FALSE, seed = NULL) {
  return(BRACoD$simulate_microbiome_counts(df, n_contributors = n_contributors, coeff_contributor = coeff_contributor, min_ab_contributor = min_ab_contributor, sd_Y = sd_Y, n_reads = n_reads, var_contributor = var_contributor, use_uniform = use_uniform, n_samples_use = n_samples_use, corr_value = corr_value, return_absolute = return_absolute, seed = seed))
}


#' Normalize OTU counts and add a pseudo count
#'
#' BRACoD requires relative abundance and cannot handle zeros, so this function
#' adds a small pseudo count (1/10th the smallest non-zero value).
#' @param df_counts A dataframe of OTU counts. Samples are rows and bacteria are columns.
#' @export
scale_counts <- function(df_counts) {
  # Check if this is legitimately counts data
  stopifnot(all(apply(df_counts, 1, function(x) all(x == as.integer(x)))))
  
  # Frequently, R and python conversion results in a "double" dataframe that has counts data. We need to fix that
  df_counts <- t(apply(df_counts, 1, function(x) as.integer(x)))

  return(BRACoD$scale_counts(df_counts))
}

#' Run the main BRACoD algorithm
#'
#' Uses pymc3 to sample the posterior of the model to determine bacteria that are
#' associated with your environmental variable.
#' @param df_relab A dataframe of relative microbiome abundances. Samples are rows and bacteria are columns.
#' @param env_var the environmnetal variable you are evaluating. You need 1 measurement associated with each sample.
#' @param n_sample number of posterior samples.
#' @param n_burn number of burn-in steps before actual sampling stops.
#' @param njobs number of parallel MCMC chains to run.
#' @export
run_bracod <- function(df_relab, env_var, n_sample=1000, n_burn=1000, njobs=4) {
  return(BRACoD$run_bracod(df_relab, env_var, n_sample = n_sample, n_burn=n_burn, njobs=njobs))
}


#' Summarize the results of BRACoD
#'
#' This summarizes the trace object that run_bracod() returns. It returns a dataframe
#' that contains two parameters of interest, the average inclusion (p) and the average
#' coefficient (beta), telling you the association between that bacteria and the environmental
#' variable
#' @param trace the pymc3 object that is the output of run_bracod()
#' @param bug_names optional, a list of names of the bacteria to include in the results
#' @param cutoff this is the cutoff on the average inclusion for inclusion
#' @export
summarize_trace <- function(trace, bug_names=NULL, cutoff=0.3) {
  return(BRACoD$summarize_trace(trace, bug_names, cutoff))
}

#' Score the results of BRACoD
#'
#' This calculate the precision, recall and F1 of your BRACoD results if you know
#' the ground truth, ie. if this is simulated data.
#' @param bugs_identified a list of integers corresponding to the indicies of the bugs you identified with BRACoD
#' @param bugs_actual a list of integers corresponding to the indicies of the bugs that truely contribute to butyrate levels
#' @export
score <- function(bugs_identified, bugs_actual) {
  return(BRACoD$score(bugs_identified, bugs_actual))
}


#' Remove NULL values in your OTU and environmental variable
#'
#' This will remove samples that are NULL in the environmental variable, as well as
#' the corresponding samples in your relative abundance data.
#' @param df_relab microbiome relative abundance data in a dataframe
#' @param Y values of the environmental variable
#' @export
remove_null <- function(df_relab, Y) {
  return(BRACoD$remove_null(df_relab, Y))
}


#' Perform convergence tests on the p and beta variables
#'
#' You may get errors are divergence of some variables after pymc3 samples the posterior.
#' We are not overly concerned about some of the variables, such as the variance, rather
#' we are really interested in the inclusion probabilities (p) and contribution coefficients
#' (beta). The convergence tests that are included here focus on evaluating those two variables.
#' @export
convergence_tests <- function(trace, sim_relab) {
  BRACoD$convergence_tests(trace, sim_relab)
}

