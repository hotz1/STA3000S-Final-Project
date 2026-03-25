library(rjags)
library(runjags)
library(here)

# Set seed for consistent results
set.seed(3000)



# Set file location and working directory
here::i_am("simulations/sparse_bayes_lasso.R")
setwd(here::here())

# Create directory for saving results
model_name <- "linear"

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

# Linear model
n = 500; p = 500
true_betas <- spike_slab_laplace(n = p, r = 0.01, rate = 1)
true_S <- which(true_betas != 0)
X <- diag(p)
Y <- X %*% true_betas + rnorm(n, sd = 1)




##### MCMC to sample from posterior p(beta | Y) #####

# Model info
lasso_post_model <- "

    # Linear regression with known variance

    model {
        # Linear regression
        for(i in 1:n){
            Y[i] ~ dnorm(inprod(X[i,], beta[]), tau)
        }
        
        # Laplace spike-and-slab prior
        for(j in 1:p){
            cat[j] ~ dcat(c(r/2, 1-r, r/2))
            slab[j] ~ dexp(rate)
            beta[j] <- (cat[j] - 2) * slab[j]
        }
        
        # Prior on r, r ~ Beta(1, p^u)
        r ~ dbeta(1, p_u)
    }
"

# Data for MCMC model
lasso_post_data <- list(
    Y = as.vector(Y),
    X = as.matrix(X),
    n = nrow(X),
    p = ncol(X),
    tau = 1,
    rate = 1,
    p_u = ncol(X)^(1.005)
)

# Initializations
lasso_init_1 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 5
)
lasso_init_2 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 6
)
lasso_init_3 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 7
)
lasso_init_4 <- list(
    cat = rep(2, p),
    slab = rep(1, p),
    r = 1/p,
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
    sample = 10000,
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


# Extract information from chains
mcmc_samples <- lasso_post_samples$mcmc
for(i in seq_along(mcmc_samples)){
    
    # Save MCMC samples for truly non-zero coefficients to disk
    mcmc_sample_df = as.data.frame(mcmc_samples[[i]])
    mcmc_sample_S <- mcmc_sample_df[, true_S]
    filename = paste0(subdir_name, "LASSO_MCMC_samples_", i, ".csv")
    write.csv(mcmc_sample_S, file = filename)
    
    # Compute posterior modes for each vector from samples
    MCMC_mode <- apply(
        mcmc_sample_df, 
        MARGIN = 2, 
        FUN = spike_slab_mode,
        slab_threshold = 0.5
    )
    MCMC_colname = paste0("MCMC.mode.", i)
    mcmc_summary_df[[MCMC_colname]] <- MCMC_mode
}

# Write the summary of the models to disk
write.csv(
    mcmc_summary_df,
    file = paste0(subdir_name, "LASSO_MCMC_summary.csv"))
