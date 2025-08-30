#!/usr/bin/env Rscript
# unlinked_loci.R

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.integer(args[1]) # Extract PBS_ARRAY_INDEX
model_name <- args[2]

set.seed(sim_index)

require(compiler)
dyn.load('aestivation3.so')
gc2<-function() {invisible(gc())}


# SIM LD FOR THE POPULATION FORWARD IN TIME
# N IS A VECTOR OF POPULATIONS SIZES, LENGTH SPECIFY THE NUMBER OF GENERATIONS TO BE SIMULATED 
# c IS A VECTOR OF RECOMBINATION RATES. DETERMINE NUMBER OF PAIRS OF LOCI
sim_pop<-function(initial_N=300, N=rep(100, 15), c=rep(0.5, 45), ncpu=4)
{
len<-length(N)
pair<-length(c)
# MATRICES FOR SUMMARY STATISTICS
p<-matrix(nc=len+1, nr=2*pair)
r2<-matrix(nc=len+1, nr=pair)
# INITIALISE. HOW?
pop<-.Call('initialise', as.double(rep(0.5, 2*pair)), as.integer(initial_N), as.integer(ncpu))
# CALCULATE INITIAL SUMMARY STATISTICS
p[,1]<-.Call('allele_freq', pop, as.integer(ncpu))
r2[,1]<-.Call('cal_r', pop, as.integer(ncpu))^2
# PROPAGATE
for (i in 1:len)
	{
	pop<-.Call('seasonal', pop, as.integer(N[i]), as.double(c), as.integer(ncpu))
	p[,i+1]<-.Call('allele_freq', pop, as.integer(ncpu))
	r2[,i+1]<-.Call('cal_r', pop, as.integer(ncpu))^2
	}
return(list(p=p, r2=r2))
}


### Define models with different population size scenarios ###
define_models <- function() {
  list(
    con = list(N = rep(30000, 72)),
    dec300 = list(N = c(rep(30000, 36), rep(300, 36))),
    dec3000 = list(N = c(rep(30000, 36), rep(3000, 36))),
    s3 = list(N = rep(c(rep(300, 3), rep(30000, 3)), times = 12)),
    s6 = list(N = rep(c(rep(300, 6), rep(30000, 6)), times = 6)),
    sd3_300 = list(N = c(rep(c(rep(300, 3), rep(30000, 3)), times = 6), rep(300, 36))),
    sd3_3000 = list(N = c(rep(c(rep(300, 3), rep(30000, 3)), times = 6), rep(3000, 36))),
    sd6_300 = list(N = c(rep(c(rep(300, 6), rep(30000, 6)), times = 3), rep(300, 36))),
    sd6_3000 = list(N = c(rep(c(rep(300, 6), rep(30000, 6)), times = 3), rep(3000, 36)))
  )
}

# Run simulation for a given model
run_simulation <- function(model) {
  sim_pop(N = model$N, c = rep(0.5, 10), ncpu = 4)
}

# Process r2 results to compute mean, standard deviation, and 83% CI for each generation
process_r2_CI <- function(r2_matrix) {
  n_generations <- ncol(r2_matrix)
  
  means <- colMeans(r2_matrix)
  sds <- apply(r2_matrix, 2, sd)
  n <- nrow(r2_matrix)
  
  # CI from t-distribution
  t_value <- qt(1 - (1 - 0.83)/2, df = n-1)
  
  stats <- data.frame(
    generation = 1:n_generations,
    mean_r2 = means,
    sd_r2 = sds,
    lower_ci = means - t_value * (sds/sqrt(n)),
    upper_ci = means + t_value * (sds/sqrt(n)),
    seed = sim_index
  )
  return(stats)
}


### Main execution ###
models <- define_models()
result <- run_simulation(models[[model_name]])

# Save results to current directory
output_file <- sprintf("%s_%d.csv", model_name, sim_index)
write.csv(process_r2_CI(result$r2), output_file, row.names = FALSE)