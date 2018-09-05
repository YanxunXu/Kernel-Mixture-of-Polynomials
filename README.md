# Kernel-Mixture-of-Polynomials
R code for posterior inference in kernel mixture of polynomials for nonparametric regression and partial linear model

Description of the files:

KMLP_simulation.R: the main R script for running posterior inference in kernel mixture of polynomials regression model;

KMLP_function.R: the classical MCMC sampler for kernel mixture of polynomials regression model;

KMLP_function_fast.R: a simplified and faster MCMC sampler for kernel mixture of polynomials regression model;

KMLP_PLM_simulation.R: the main R script for running posterior inference in kernel mixture of polynomials partial linear model;

KMLP_PLM_function.R: the MCMC sampler for kernel mixture of polynomials partial linear model;

kernel_weights_evaluation.cpp: additional RCpp file supporting the above MCMC samplers.
