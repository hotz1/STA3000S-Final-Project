library(here)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(latex2exp)
library(scales)
library(kableExtra)

# Set file location and working directory
here::i_am("simulations/plot_sim_results.R")
setwd(here::here())

##### Diagnostic information for all the simulations #####

# Data-generating models
model_name <- c("linear", "poisson", "negbin")
model_clean <- c("Linear Model", "Poisson Model", "Negative Binomial Model")
model_tbl <- tibble(model_name, model_clean)

# Prior distribution on betas
prior_name <- c("laplace")
prior_clean <- c("Laplace Spike-and-Slab Prior")
prior_tbl <- tibble(prior_name, prior_clean)

# Aspect ratios (ratio of n/p)
aspect_name <- c("ten_percent", "fifty_percent")
aspect_clean <- c("10%", "50%")
aspect_ratio <- c(0.1, 0.5)
aspect_tbl <- tibble(aspect_name, aspect_clean, aspect_ratio)

# Simulation settings (p, n, s)
p.sims <- c(100, 200, 500, 1000)
s.sims <- c(1, 2, 5, 10)
sim_size_tbl <- tibble(p.sims, s.sims)

# Combine information into a larger table
sim_info_tbl <- crossing(model_tbl, prior_tbl, aspect_tbl) %>%
    mutate(model_prior = paste(model_name, prior_name, sep = "_"))




##### First Visual: Contingency Tables #####

for(i in seq_len(nrow(sim_info_tbl))){
    sim_info <- sim_info_tbl[i, ]
    
    # Get folder containing results, set output name and folder
    results_folder <- here(
        "simulations", "results", sim_info$model_prior, sim_info$aspect_name
    )
    out_dir <- here(
        "simulations", "plots", sim_info$model_prior, sim_info$aspect_name
    )
    if(!dir.exists(out_dir)){
        dir.create(out_dir, recursive = TRUE)
    }
    
    # Read in summary files
    MCMC_summaries <- list()
    summary_files <- list.files(
        path = results_folder,
        pattern = ".*summary.csv", 
        full.names = T,
        recursive = T)
    
    # Create combined plots of contingency tables
    CI_contingency_plots <- list()
    median_contingency_plots <- list()
    for(f in seq_along(summary_files)){
        MCMC_summaries[[f]] <- read_csv(
            summary_files[f],
            name_repair = "minimal"
        ) %>%
            rename(CoefName = 1)
        
        MCMC_summaries[[f]] <- MCMC_summaries[[f]] %>%    
            mutate(
                GroundTruth = factor(
                    if_else(TrueCoef == 0, "True Zero", "True Nonzero"),
                    levels = c("True Nonzero", "True Zero")
                ),
                CI_Selection = factor(
                    if_else((0 >= Lower95 & 0 <= Upper95), "Discarded", "Maintained"),
                    levels = c("Discarded", "Maintained")
                ),
                Med_Selection = factor(
                    if_else(Median == 0, "Discarded", "Maintained"),
                    levels = c("Discarded", "Maintained")
                )
            )
        
        CI_contingency <- MCMC_summaries[[f]] %>%
            count(GroundTruth, CI_Selection) %>% 
            complete(GroundTruth, CI_Selection, fill = list(n = 0)) %>%
            mutate(
                Correct = case_when(
                    GroundTruth == "True Nonzero" & CI_Selection == "Maintained" ~ "Correct", 
                    GroundTruth == "True Zero" & CI_Selection == "Discarded" ~ "Correct",
                    .default = "Incorrect"
                ),
                Rate_Name = case_when(
                    GroundTruth == "True Nonzero" & CI_Selection == "Maintained" ~ "TPR",
                    GroundTruth == "True Zero" & CI_Selection == "Maintained" ~ "FPR",
                    GroundTruth == "True Nonzero" & CI_Selection == "Discarded" ~ "FNR",
                    GroundTruth == "True Zero" & CI_Selection == "Discarded" ~ "TNR",
                )) %>%
            group_by(GroundTruth) %>%
            mutate(Rate = n/sum(n)) %>%
            ungroup()
        
        median_contingency <- MCMC_summaries[[f]] %>%
            count(GroundTruth, Med_Selection) %>% 
            complete(GroundTruth, Med_Selection, fill = list(n = 0)) %>%
            mutate(
                Correct = case_when(
                    GroundTruth == "True Nonzero" & Med_Selection == "Maintained" ~ "Correct", 
                    GroundTruth == "True Zero" & Med_Selection == "Discarded" ~ "Correct",
                    .default = "Incorrect"
                ),
                Rate_Name = case_when(
                    GroundTruth == "True Nonzero" & Med_Selection == "Maintained" ~ "TPR",
                    GroundTruth == "True Zero" & Med_Selection == "Maintained" ~ "FPR",
                    GroundTruth == "True Nonzero" & Med_Selection == "Discarded" ~ "FNR",
                    GroundTruth == "True Zero" & Med_Selection == "Discarded" ~ "TNR",
                )) %>%
            group_by(GroundTruth) %>%
            mutate(Rate = n/sum(n)) %>%
            ungroup()
        
        CI_contingency_plots[[f]] <- CI_contingency %>% 
            ggplot(aes(x = GroundTruth, y = CI_Selection)) +
            geom_tile(aes(fill = Correct), colour = "black") +
            geom_text(aes(
                label = paste0(
                    n, " (", Rate_Name, ": ",
                    scales::percent(Rate, 0.1), ")"
                ))
            ) +
            scale_fill_manual(
                values = c("Correct" = "#00ff0080", "Incorrect" = "#ff000060")
            ) +
            labs(title = TeX(paste0(
                "$p$ = ", p.sims[f], 
                ", $n$ = ", p.sims[f] * sim_info$aspect_ratio, 
                ", $|S_{\\beta_{0}}|$ = ", s.sims[f],
                sep = ""))
            ) +
            theme_minimal() +
            theme(panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5, vjust = -5),
                  legend.position = "none",
                  axis.title = element_blank(),
                  axis.text.x = element_text(vjust = 5),
                  axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = -5))
        
        median_contingency_plots[[f]] <- median_contingency %>% 
            ggplot(aes(x = GroundTruth, y = Med_Selection)) +
            geom_tile(aes(fill = Correct), colour = "black") +
            geom_text(aes(
                label = paste0(
                    n, " (", Rate_Name, ": ",
                    scales::percent(Rate, 0.1), ")"
                ))
            ) +
            scale_fill_manual(
                values = c("Correct" = "#00ff0080", "Incorrect" = "#ff000080")
            ) +
            labs(title = TeX(paste0(
                "$p$ = ", p.sims[f], 
                ", $n$ = ", p.sims[f] * sim_info$aspect_ratio, 
                ", $|S_{\\beta_{0}}|$ = ", s.sims[f],
                sep = ""))
            ) +
            theme_minimal() +
            theme(panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5, vjust = -5),
                  legend.position = "none",
                  axis.title = element_blank(),
                  axis.text.x = element_text(vjust = 5),
                  axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = -5))
    }
    
    CI_contingency_combined <- patchwork::wrap_plots(CI_contingency_plots, nrow = 2) +
        plot_annotation(
            title = TeX("Variable Selection using Posterior 95% Credible Interval of $\\beta$"),
            subtitle = paste0(
                sim_info$model_clean, " with ", 
                sim_info$prior_clean, " and aspect ratio of ", 
                sim_info$aspect_clean
            ),
            theme = theme(
                plot.title = element_text(hjust = 0.5, size = 16),
                plot.subtitle = element_text(hjust = 0.5, size = 14)
            )
        )
    ggsave(filename = "CI_contingency_plots_combined.pdf", 
           plot = CI_contingency_combined,
           path = out_dir,
           width = 8,
           height = 8,
           units = "in")
    
    
    median_contingency_combined <- patchwork::wrap_plots(median_contingency_plots, nrow = 2) +
        plot_annotation(
            title = TeX("Variable Selection using Posterior Median of $\\beta$"),
            subtitle = paste0(
                sim_info$model_clean, " with ", 
                sim_info$prior_clean, " and aspect ratio of ", 
                sim_info$aspect_clean
            ),
            theme = theme(
                plot.title = element_text(hjust = 0.5, size = 16),
                plot.subtitle = element_text(hjust = 0.5, size = 14)
            )
        )
    ggsave(filename = "median_contingency_plots_combined.pdf", 
           plot = median_contingency_combined,
           path = out_dir,
           width = 8,
           height = 8,
           units = "in")    
}
    
    
   
 
##### Second Visual: Histogram of Dimensionality #####

for(i in seq_len(nrow(sim_info_tbl))){
    sim_info <- sim_info_tbl[i, ]
    
    # Get folder containing results, set output name and folder
    results_folder <- here(
        "simulations", "results", sim_info$model_prior, sim_info$aspect_name
    )
    out_dir <- here(
        "simulations", "plots", sim_info$model_prior, sim_info$aspect_name
    )
    if(!dir.exists(out_dir)){
        dir.create(out_dir, recursive = TRUE)
    }
    
    # Read in summary files
    MCMC_sparsities <- list()
    sparsity_files <- list.files(
        path = results_folder,
        pattern = ".*sparsity.csv", 
        full.names = T,
        recursive = T)
    
    # Create combined plots of number of non-zero betas
    local_dimension_plots <- list()
    wider_dimension_plots <- list()
    for(f in seq_along(sparsity_files)){
        MCMC_sparsities[[f]] <- read_csv(
            sparsity_files[f],
            col_select = c(-1)
        )
        
        # Plots which just have the dimensions
        local_dimension_plots[[f]] <- MCMC_sparsities[[f]] %>% 
            ggplot() +
            geom_histogram(
                aes(x = sparse_dim), binwidth = floor(sqrt(s.sims[f])),
                colour = "black", fill = "goldenrod2") +
            labs(x = TeX("Size of $S_{\\beta}$"), 
                 y = "Count",
                 title = TeX(paste0(
                     "p = ", p.sims[f], 
                     ", n = ", p.sims[f] * sim_info$aspect_ratio, 
                     ", $|S_{\\beta_{0}}|$ = ", s.sims[f],
                     sep = "")
                 )) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
        
        # Includes the true value as a dotted line
        wider_dimension_plots[[f]] <- local_dimension_plots[[f]] + 
            geom_vline(
                xintercept = s.sims[f],
                colour = "red", linetype = "dashed"
            ) + 
            scale_x_continuous(
                breaks = function(x) sort(unique(c(0, scales::breaks_pretty()(x)))),
                limits = function(x) range(c(0, x))
            )
    }
    
    # Combine plots together and save to disk
    local_dim_combined <- patchwork::wrap_plots(local_dimension_plots, nrow = 2) +
        plot_annotation(
            title = TeX("Size of $S_{\\beta}$ in Posterior MCMC Samples"),
            subtitle = paste0(
                sim_info$model_clean, " with ", 
                sim_info$prior_clean, " and aspect ratio of ", 
                sim_info$aspect_clean
            ),
            theme = theme(
                plot.title = element_text(hjust = 0.5, size = 16),
                plot.subtitle = element_text(hjust = 0.5, size = 14)
            )
        )
    ggsave(filename = "local_dimensionality_combined.pdf", 
           plot = local_dim_combined,
           path = out_dir,
           width = 8,
           height = 8,
           units = "in")
    
    
    wider_dim_combined <- patchwork::wrap_plots(wider_dimension_plots, nrow = 2) +
        plot_annotation(
            title = TeX("Size of $S_{\\beta}$ in Posterior MCMC Samples"),
            subtitle = paste0(
                sim_info$model_clean, " with ", 
                sim_info$prior_clean, " and aspect ratio of ", 
                sim_info$aspect_clean
            ),
            theme = theme(
                plot.title = element_text(hjust = 0.5, size = 16),
                plot.subtitle = element_text(hjust = 0.5, size = 14)
            )
        )
    ggsave(filename = "wider_dimensionality_combined.pdf", 
           plot = wider_dim_combined,
           path = out_dir,
           width = 8,
           height = 8,
           units = "in")
}




##### Third Visual: Tables of RMSE for Train and Test Datasets #####

for(i in seq_len(nrow(sim_info_tbl))){
    sim_info <- sim_info_tbl[i, ]
    
    # Get folder containing results, set output name and folder
    results_folder <- here(
        "simulations", "results", sim_info$model_prior, sim_info$aspect_name
    )
    out_dir <- here(
        "simulations", "plots", sim_info$model_prior, sim_info$aspect_name
    )
    if(!dir.exists(out_dir)){
        dir.create(out_dir, recursive = TRUE)
    }
    
    # Read in summary files
    MCMC_predictions <- list()
    prediction_files <- list.files(
        path = results_folder,
        pattern = ".*predictions.csv", 
        full.names = T,
        recursive = T)
    
    # Create tables of information
    RMSE_tables <- list()
    for(f in seq_along(prediction_files)){
        MCMC_predictions[[f]] <- read_csv(
            prediction_files[f],
            col_select = c(-1)
        )
        preds = MCMC_predictions[[f]]
        train_RMSE <- sqrt(mean((preds$Y_train_true - preds$Y_train_pred)^2))
        test_RMSE <- sqrt(mean((preds$Y_test_true - preds$Y_test_pred)^2))
        
        RMSE_tables[[f]] <- data.frame(
            `Model` = sim_info$model_clean, 
            `Prior` = sim_info$prior_clean, 
            `Sample Size` = p.sims[f] * sim_info$aspect_ratio,
            `Parameters` = p.sims[f],
            `Aspect Ratio` = sim_info$aspect_clean,
            `Train RMSE` = sprintf("%.4g", train_RMSE), 
            `Test RMSE` = sprintf("%.4g", test_RMSE)
        )
    }
    
    # Final table
    all_RMSE_tbl <- bind_rows(RMSE_tables)
    all_RMSE_tbl %>%
        kable(format = "html", 
              col.names = stringr::str_replace(colnames(all_RMSE_tbl), "\\.", " ")) %>%
        kable_styling(latex_options = c("striped", "hold_position")) %>%
        save_kable(file = paste0(out_dir, "/", "RMSE_table.jpeg"))
}

