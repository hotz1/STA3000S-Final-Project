library(rjags)
library(runjags)
library(here)

# Set seed for consistent results
# set.seed(3001)


# Settings for the compute cluster
Sys.setenv(
    OMP_NUM_THREADS = 1,
    MKL_NUM_THREADS=1,
    BLAS_NUM_THREADS=1,
    LAPACK_NUM_THREDS=1)



# Set file location and working directory
here::i_am("simulations/simulate_negbin_laplace.R")
setwd(here::here())


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
n_sims <- as.integer(args[1])
n <- as.integer(args[2])
p <- as.integer(args[3])
sparse_prop <- as.numeric(args[4])



# Create directory for saving results
model_name <- "negbin_laplace"

result_dir <- paste0("./simulations/results/", model_name, "/")
n_subdir <- length(list.dirs(path = result_dir, recursive = FALSE))
subdir_name <- paste0(result_dir, "LASSO_MCMC_", n_subdir + 1, "/")

if(!dir.exists(subdir_name)){
    dir.create(subdir_name, recursive = TRUE)
}




##### Simulate the data #####

# Function for the prior 
spike_slab_laplace <- function(n, r, rate){
  # Simulate from Laplace spike-and-slab with probability 
  spikes <- sample(c(-1, 0, 1), size = n, 
                   replace = T, prob = c(r/2, 1-r, r/2))
  exp_sim <- rexp(n, rate = rate)
  return(spikes * exp_sim)
}

# Generate data from negative binomial model
true_betas <- spike_slab_laplace(p, sparse_prop, rate = 1/5)
true_S <- which(true_betas != 0)

true_phi <- 1

X_train <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
Y_train <- rnbinom(n, size = true_phi, mu = exp(X_train %*% true_betas))

X_test <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
Y_test <- rnbinom(n, size = true_phi, mu = exp(X_test %*% true_betas))




##### MCMC to sample from posterior p(beta | Y) #####

# Model info
lasso_post_model <- "

    # Negative binomial regression with unknown dispersion parameter

    model {

        for(i in 1:n){
            # Log link for the model
            log(mu[i]) <- inprod(X[i,], beta[])
            
            prob[i] <- phi / (phi + mu[i])
            Y[i] ~ dnegbin(prob[i], phi)
        }
        
        # Laplace spike-and-slab prior
        for(j in 1:p){
            cat[j] ~ dcat(c(r/2, 1-r, r/2))
            slab[j] ~ dexp(rate)
            beta[j] <- (cat[j] - 2) * slab[j]
        }
        
        # Prior on r, r ~ Beta(1, p^u)
        r ~ dbeta(1, p_u)
        
        # Prior on phi 
        phi ~ dgamma(0.01, 0.01)
    }
"

# Data for MCMC model
lasso_post_data <- list(
    Y = as.vector(Y_train),
    X = as.matrix(X_train),
    n = n,
    p = p,
    rate = 1/5,
    p_u = p^(1.005)
)

# Initializations
lasso_init_1 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    phi = 1,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 5
)
lasso_init_2 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    phi = 1,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 6
)
lasso_init_3 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    phi = 1,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 7
)
lasso_init_4 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    phi = 1,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 8
)

# Save model configs
lasso_compiled <- run.jags(
    model = lasso_post_model,
    monitor = c("beta"), 
    data = lasso_post_data,
    inits = list(lasso_init_1, lasso_init_2, lasso_init_3, lasso_init_4),
    n.chains = 4,
    sample = 0,
    method = "parallel"
)
write.jagsfile(lasso_compiled, file = paste0(subdir_name, "model.txt"))




##### Run the MCMC sampler #####
lasso_post_samples <- run.jags(
    model = lasso_post_model,
    monitor = c("beta"), 
    data = lasso_post_data,
    inits = list(lasso_init_1, lasso_init_2, lasso_init_3, lasso_init_4),
    n.chains = 4,
    sample = n_sims,
    method = "parallel"
)




##### Write results locally #####


# Function to compute mode of a vector from a spike-and-slab 
# The spike threshold is minimum proportion of non-zero 
spike_slab_mode <- function(samples, slab_threshold){
    p_slab <- mean(samples != 0)
    if (p_slab < slab_threshold){
        return (0)
    }

    slab <- samples[samples != 0]
    d_slab <- density(slab)
    post_mode <- d_slab$x[which.max(d_slab$y)]
    return (post_mode)
}


# Summary statistics for the chains
mcmc_summary_df <- as.data.frame(summary(lasso_post_samples))
mcmc_summary_df <- cbind(TrueCoef = true_betas, mcmc_summary_df)

# Extract the chains
mcmc_sample_df <- as.matrix(as.mcmc.list(lasso_post_samples))
mcmc_sample_S <- mcmc_sample_df[, true_S]

write.csv(
    mcmc_sample_S, 
    file = paste0(subdir_name, "LASSO_MCMC_samples.csv")
)

# Compute posterior modes for each parameter from samples
MCMC_mode <- apply(
    mcmc_sample_df,
    MARGIN = 2,
    FUN = spike_slab_mode,
    slab_threshold = 0.5
)
mcmc_summary_df[["MCMC.mode"]] <- MCMC_mode

# Write the summary of the models to disk
write.csv(
    mcmc_summary_df,
    file = paste0(subdir_name, "LASSO_MCMC_summary.csv")
)


# Compute posterior modes for Y_train and Y_test
Y_train_mode <- exp(X_train %*% MCMC_mode)
Y_test_mode <- exp(X_test %*% MCMC_mode)

post_pred_df <- data.frame(
    Y_train_true = Y_train,
    Y_train_pred = Y_train_mode,
    Y_test_true = Y_test,
    Y_test_pred = Y_test_mode
)

write.csv(
    post_pred_df,
    file = paste0(subdir_name, "LASSO_predictions.csv")
)
