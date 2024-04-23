# RESULTS 

library(tidyverse)
source("functions.R")

# SET PARAMETERS
max_iter = 20
num_data_pts = 20
perc_param = seq(0, num_data_pts - 1, by = 1) / num_data_pts
cycle_lengths = c(2,3,4)
num_samples_MC = 15

# CHOOSE NETWORK
adj_matrix = erdos_renyi(50, 0.05)

# FUNCTION DEFNs
get_results <- function(adj_matrix, perc_param, max_iter, cycle_lengths, num_samples_MC) {
  results <- vector("list", length = length(cycle_lengths))
  for (r in 1:length(cycle_lengths)){
    print(paste0("r ", r))
    results[[r]] <- vector("list", length = length(perc_param))
    for (p in 1:length(perc_param)) {
      print(paste0("p ", p))
      alg2_out <- alg2(adj_matrix, perc_param[p], max_iter, cycle_lengths[r], num_samples_MC)
      pr_not_in_C <- prob_i_not_in_C_2(adj_matrix, perc_param[p], cycle_lengths[r], alg2_out[[max_iter+1]], num_samples_MC)
      size_giant <- size_giant(pr_not_in_C)
      size_small <- mean(size_cluster_with_i(alg2_out[[max_iter+1]], max_iter, perc_param[p], adj_matrix))
      results[[r]][[p]] <- list("size_giant" = size_giant, "size_small" = size_small)
    }
  }
  return(results)
}

result_to_dataframe <- function(result, perc_param, cycle_lengths) {
  df <- data.frame("r" = integer(), "perc_param" = numeric(), "size_giant" = numeric(), "size_small" = numeric(), stringsAsFactors = FALSE)
  for (r in 1:length(result)) {
    for (s in 1:length(result[[r]])) {
      # Create a new row and bind it to the existing data frame
      new_row <- data.frame(
        r = cycle_lengths[r],
        perc_param = perc_param[s],
        size_giant = result[[r]][[s]]$size_giant,
        size_small = result[[r]][[s]]$size_small,
        stringsAsFactors = FALSE
      )
      df <- rbind(df, new_row)
    }
  }
  return(df)
}

# GET RESULTS
start_time <- proc.time()
result <- get_results(adj_matrix, perc_param, max_iter, cycle_lengths, num_samples_MC)
end_time <- proc.time() - start_time
print("Time elapsed: ", end_time)

result_df <- result_to_dataframe(result, perc_param, cycle_lengths)
save(result_df, file = paste0("results/result_df_ER_50_05_2.RData"))

# PLOT RESULTS
plot1 <- result_df %>%
  ggplot(aes(x = perc_param, y = size_giant, color = as.factor(r), group = r)) +  # Convert r to factor for coloring
  geom_line() +
  labs(title = "Size of Giant Component vs. Percolation Probability",
       x = "Percolation Probability (p)",
       y = "Size of Giant Component",
       color = "Cycle Length (r)") +
  theme_minimal()

plot1

ggsave("plots/plot_results_ER_50_05_2.png", plot = plot1, width = 10, height = 6, dpi = 300)





