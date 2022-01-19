require(tidyverse)
require(foreach)


source('process_slim_scripts.r')

path_to_hxb2_env <- "~/Dropbox/active/hiv_recombination_rate/dat/hxb2/hxb2_env.fa"

mu <- 1e-5
rho <- 1e-5 * seq(1, 20, length.out = 3)
Ne <- 10000
M <- 1000
reps <- 1:3

base_directory <- "2022_01_19"
samplingGenerations <- 5 * Ne + c(0, 150, 300, 600)
params <- expand.grid(mu, rho, Ne, M, reps)
colnames(params) <- c("mu", "rho", "Ne", "M", "rep")

setup_slim_trials(base_directory, params, samplingGenerations)
run_slim_trials(base_directory)
#Convert from slim output into Zanini et al output. This is by far the slowest step
process_slim_trials(base_directory)



