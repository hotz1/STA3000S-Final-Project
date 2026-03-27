# STA3000S Final Project

A repository which contains R code for running and plotting the results of Markov Chain Monte Carlo algorithms which are used to sample from a specific posterior distribution.

This project is inspired by the paper [*Bayesian Linear Regression with Sparse Priors*](https://arxiv.org/pdf/1403.0735) by Castillo et al., and the express purpose of the code in this repository is to experimentally verify some of the theoretical results in the authors' paper, alongside extending these simulations to some non-linear regression settings. 

# Code

There are three main R files for running the simulations; `simulate_linear_laplace.R`, `simulate_poisson_laplace.R`, and `simulate_negbin_laplace.R`, all three of which are stored in the `simulations` sub-directory. Using `Rscript`, these files can be run from the command-line with certain **required** positional arguments. For example, to run a simulation with the linear model from the root of this repository, one can type the following command into the terminal:
```
Rscript simulations/simulate_linear_laplace.R --n_sims 30000 --n 10 --p 100 --p_true 1 --seed 20
```
The results of executing this code will be saved to disk in the `results` folder, with the sub-directory depending on the choice of regression model. Do note that the default place for data to be stored does not match this current repository structure, as the result files were manually rearranged at a later point for more thorough organization.
