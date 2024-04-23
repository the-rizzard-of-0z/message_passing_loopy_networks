library(igraph)
library(expm)

# Function Definitions: 
erdos_renyi <- generate_erdos_renyi_gnp <- function(n, p) {
  # n: Number of vertices
  # p: Probability of edge creation between any two vertices
  graph <- erdos.renyi.game(n, p, type = "gnp")
  adj <- as_adjacency_matrix(graph, sparse = FALSE)  # Set sparse=FALSE to get a regular matrix
  return(adj)
}

get_initial_u <- function(adj_matrix) {
  initial_u <- matrix(0.5, nrow(adj_matrix), nrow(adj_matrix))
  return(initial_u)
}

# For Alg 1
prob_i_not_in_C <- function(adj_matrix, u, via_j, perc_param = NULL) {
  if (!via_j) {
    result <- rep(0, ncol(u))
    for (i in 1:length(result)){
      result[i] <- prod(u[,i][adj_matrix[i,] == 1])  
    }
  } else {
    result <- matrix(0, nrow = nrow(adj_matrix), ncol = nrow(adj_matrix))
    for (i in 1:nrow(adj_matrix)){
      for (j in 1:ncol(adj_matrix)) {
          temp_adj_matrix <- adj_matrix
          temp_adj_matrix[i, j] <- 0
          result[i,j] <- prod(u[,i][temp_adj_matrix[i,] == 1])
      }
    }
  }
  return(result)
}

size_giant <- function(v) {
  result <- (1/length(v))*sum(1-v)
  return(result)
}

quotient_sum <- function(adj_matrix, u, u_prime, via_j) {
  if (!via_j) {
    result <- rep(0, ncol(u))
    for (i in 1:length(result)){
      quotient_u <- u_prime / u
      result[i] <- sum(quotient_u[,i][adj_matrix[i,] == 1])
    }
  } else {
    result <- matrix(ncol = ncol(u), nrow = nrow(u))
    for (i in 1:ncol(result)) {
      for (j in 1:nrow(result)) { 
        quotient_u <- u_prime / u
        temp_adj_matrix <- adj_matrix
        temp_adj_matrix[i, j] <- 0
        result[i,j] <- sum(quotient_u[,i][temp_adj_matrix[i,] == 1])
      }
    }
  }
  return(result)
}

# size of cluster containing i | i does not belong to perc_cluster

size_cluster_with_i <- function(u, max_iter, p_perc, adj_matrix) {
  u_prime_list <- list(get_initial_u(adj_matrix))
  for (iter in 1:max_iter) {
    temp = p_perc * (1 + quotient_sum(adj_matrix, u, u_prime_list[[iter]], via_j = TRUE)) * prob_i_not_in_C(adj_matrix, u, via_j = TRUE)
    u_prime_list <- c(u_prime_list, list(temp))
  }
  result <- (1 + quotient_sum(adj_matrix, u, u_prime_list[[max_iter+1]], via_j = FALSE)) / length(u)
  return(result)
}

# Function to create the Hashimoto matrix from an adjacency matrix
get_hashimoto <- function(adj) {
  n <- nrow(adj)  # Number of vertices
  m <- sum(adj) / 2  # Number of edges
  
  # Initialize the Hashimoto matrix with 0s
  B <- matrix(0, nrow = 2 * m, ncol = 2 * m)
  
  # List of directed edges
  edges <- which(adj == 1, arr.ind = TRUE)
  
  # Map edges to their position in the list for fast lookup
  edge_index <- setNames(seq_along(edges[, 1]), paste(edges[, 1], edges[, 2], sep = "-"))
  
  # Fill the Hashimoto matrix
  for (i in seq_along(edges[, 1])) {
    u <- edges[i, 1]
    v <- edges[i, 2]
    # Find all vertices w where (v, w) exists and w != u
    outgoing_edges <- which(adj[v, ] == 1 & seq_len(n) != u)
    for (w in outgoing_edges) {
      # Look up the edge index in the edge list
      j <- edge_index[paste(v, w, sep = "-")]
      if (!is.null(j)) {
        B[i, j] <- 1
      }
    }
  }
  return(B)
}

alg1 <- function(adj_matrix, perc_param, max_iter) {
  # initialize
  u_list <- list(get_initial_u(adj_matrix))
  # iterate
  for (iter in 1:max_iter) {
    # u
    temp = 1 - perc_param + perc_param * prob_i_not_in_C(adj_matrix, u_list[[iter]], via_j = TRUE)
    u_list <- c(u_list, list(temp))
  }
  # stability of fixed pt. 1
  hashimoto_matrix <- get_hashimoto(adj_matrix)
  leading_eigenvalue <- eigen(hashimoto_matrix)$values[which.max(abs(eigen(hashimoto_matrix)$values))]
  p_threshold = 1 / (leading_eigenvalue)
  return(list(u_list = u_list, p_threshold = p_threshold))
}

# For Alg 2 
get_N_i <- function(adj_matrix, i, r, via_j = FALSE, j = NULL) {
  n <- nrow(adj_matrix)  # Number of nodes in the adjacency matrix
  # Function to find nodes involved in loops up to length r starting and ending at a given node
  find_loops_up_to_r <- function(node) {
    nodes_in_cycles <- logical(n)  # Boolean vector to keep track of nodes involved in loops
    # Check for loops of each length from 1 up to r
    for (k in 1:r) {
      adj_matrix_k <- adj_matrix %^% k  # k-th power of the adjacency matrix
      # Loop through all nodes to check if they are part of a loop of length k with the start node
      if (adj_matrix_k[node, node] > 0) {  # There's a loop of length k starting and ending at 'node'
        for (j in 1:n) {
          if (adj_matrix_k[node, j] > 0 && adj_matrix_k[j, node] > 0) {
            nodes_in_cycles[j] <- TRUE
          }
        }
      }
    }
    return(nodes_in_cycles)
  }
  # Find nodes in loops involving the given node 'i'
  nodes_in_loops_i <- find_loops_up_to_r(i)
  if (via_j) {
    if (is.null(j)) {
      stop("Parameter 'j' must be specified if 'via_j' is TRUE.")
    }
    # Find nodes in loops involving another node 'j'
    nodes_in_loops_j <- find_loops_up_to_r(j)
    # Return the difference: nodes in loops with 'i' but not with 'j' 
    nodes_in_loops_i_only <- which(nodes_in_loops_i & !nodes_in_loops_j)
    return(nodes_in_loops_i_only)
  } else {
    # Return nodes in loops with 'i' only excluding i
    return(nodes_in_loops_i)
  }
}


find_cycle_nodes <- function(adj_matrix, r, i) {
  cycle_nodes <- integer()
  cycle_nodes <- union(cycle_nodes, i)
  paths <- list()
  for (j in which(adj_matrix[i, ] == 1)) {
    paths <- c(paths, list(j))
  }
  for (k in 2:r) {
    if (length(paths) == 0) {
      return(cycle_nodes)
    }
    new_paths <- list()
    for (h in 1:length(paths)) {
      # print(paste0("h", h))
      # print(paste0(length(paths)))
      # print(paste0("paths", paths))
      if (adj_matrix[i,paths[[h]][k-1]] == 1) {
        cycle_nodes <- union(cycle_nodes, paths[[h]])
      }
      # edit paths to become length 3 (ignore i this time)
      for (j in which(adj_matrix[paths[[h]][k-1], ] == 1)) {
        if (j != i && !(j %in% paths[[h]])) {
          new_paths <- c(new_paths, list(c(paths[[h]], j)))
        }
      }
    }
    paths <- new_paths
  }
  return(cycle_nodes)
}

get_N_i <- function(adj_matrix, i, r, via_j = FALSE, j = NULL) {
  n <- nrow(adj_matrix)  # Number of nodes in the adjacency matrix
  nodes_in_loops_i <- find_cycle_nodes(adj_matrix, r, i)
  if (via_j) {
    if (is.null(j)) {
      stop("Parameter 'j' must be specified if 'via_j' is TRUE.")
    }
    nodes_in_loops_j <- find_cycle_nodes(adj_matrix, r, j)
    # Debug output
    # cat("Nodes in loops with i:", nodes_in_loops_i, "\n")
    # cat("Nodes in loops with j:", nodes_in_loops_j, "\n")
    nodes_in_loops_i_only <- setdiff(nodes_in_loops_i, nodes_in_loops_j)
    # cat("nodes_in_loops_i_only:", nodes_in_loops_i_only, "\n")
    return(nodes_in_loops_i_only)
  } else {
    # cat("Nodes in loops with i excluding i:", sum(nodes_in_loops_i), "\n")
    return(nodes_in_loops_i)
  }
}

get_T_MC <- function(adj_matrix, neighbors, perc_param, i, num_samples){
  submatrix <- adj_matrix[union(i, neighbors), union(i, neighbors)]
  edges <- which(submatrix == 1, arr.ind = TRUE)
  if (length(edges) == 0) {
    return(list(configurations = list(), pr = numeric()))
  }
  num_edges <- nrow(edges)
  configurations <- vector("list", num_samples)
  config_probabilities <- numeric(num_samples)
  # Iterate over the number of samples instead of all possible configurations
  for (sample in 1:num_samples) {
    # Sample each edge's existence based on perc_param
    binary_config <- rbinom(num_edges, 1, perc_param)
    config_matrix <- matrix(0, nrow = length(union(i, neighbors)), ncol = length(union(i, neighbors)))
    # Build the configuration matrix for the current sample
    for (e in 1:num_edges) {
      config_matrix[edges[e, "row"], edges[e, "col"]] <- binary_config[e]
    }
    # Since the matrix is symmetric
    config_matrix <- config_matrix + t(config_matrix)
    # Diagonal should be 0 if self-loops are not considered
    diag(config_matrix) <- 0
    configurations[[sample]] <- config_matrix
    # Calculate the probability of this configuration
    config_probabilities[sample] <- prod(binary_config * perc_param + (1 - binary_config) * (1 - perc_param))
  }
  return(list(configurations = configurations, pr = config_probabilities))
}

sigma_power_product <- function(N_i_to_j, T_i_to_j, u, i) { #not a good choice of function name but whatever
  temp <- 1
  for (k in N_i_to_j) {
    temp <- temp * (u[which(union(i, N_i_to_j) == k), which(union(i, N_i_to_j) == i)]^sigma(N_i_to_j, T_i_to_j, i, k))
  }
  return(temp)
}

sigma <- function(N_i_to_j, T_i_to_j, i, k) {
  num_nodes <- nrow(T_i_to_j)
  reach_matrix <- T_i_to_j
  for (step in 1:(num_nodes-1)) {
    reach_matrix <- reach_matrix %*% T_i_to_j
    if (reach_matrix[which(union(i, N_i_to_j) == i), which(union(i, N_i_to_j) == k)] > 0) {
      return(1)
    }
  }
  return(0)
}

get_temp <- function(adj_matrix, perc_param, r, u, num_samples_MC, N) {
  temp = matrix(0, nrow = nrow(adj_matrix), ncol = ncol(adj_matrix))
  for (i in 1:nrow(adj_matrix)) {
    for (j in 1:ncol(adj_matrix)) {
      # print(paste0("i: ", i, "j: ", j))
      N_i_to_j <- N[[i]][[j]]
      T_i_to_j <- get_T_MC(adj_matrix, N_i_to_j, perc_param, i, num_samples_MC)
      if (length(T_i_to_j$configurations) == 0){
        temp[i,j] <- 1
      } else {
        for (z in 1:length(T_i_to_j$configurations)) {
          temp[i,j] <- temp[i,j] + sigma_power_product(N_i_to_j, T_i_to_j$configurations[[z]], u, i)
        }
        temp[i,j] <- temp[i,j] / length(T_i_to_j$configurations) # same as num_samples_MC
        # print(paste0("temp[i,j] ", temp[i,j]))
      }
    }
  }
  return(temp)
}

# r is assumed max primitive cycle length 
alg2 <- function(adj_matrix, perc_param, max_iter, r, num_samples_MC) {
  # initialize
  u_list <- list(get_initial_u(adj_matrix))
  # get neighbors 
  N <- vector("list", nrow(adj_matrix))
  for (i in 1:nrow(adj_matrix)){
    for(j in 1:ncol(adj_matrix)) {
      N[[i]][[j]] = get_N_i(adj_matrix, i, r, via_j = TRUE, j)
    }
  }
  # iterate
  for(i in 1:max_iter) {
    print(paste0("iter: ", i))
    temp <- get_temp(adj_matrix, perc_param, r, u_list[[i]], num_samples_MC, N)
    u_list <- c(u_list, list(temp))
  }
  return(u_list)
}

get_prob <- function(adj_matrix, perc_param, r, u, num_samples_MC, j) {
  N_j <- get_N_i(adj_matrix, j, r, via_j = FALSE)
  T_j <- get_T_MC(adj_matrix, N_j, perc_param, j, num_samples_MC)
  temp = 0
  if (length(T_j$configurations) == 0){
    return(1)
  }
  for (z in 1:length(T_j$configurations)) {
    # print(paste0("z", z))
    # print(paste0(length(T_j$configurations)))
    temp <- temp + sigma_power_product(N_j, T_j$configurations[[z]], u, j)
  }
  temp <- temp / length(T_j$configurations) # same as num_samples_MC
  return(temp)
}

prob_i_not_in_C_2 <- function(adj_matrix, perc_param, r, u, num_samples_MC) {
  pr <- rep(0, nrow(adj_matrix))
  for (j in 1:nrow(adj_matrix)) {
    pr[j] <- get_prob(adj_matrix, perc_param, r, u, num_samples_MC, j)
  }
  return(pr)
}


# # replaced my get_T_MC for convienie
# get_T <- function(adj_matrix, neighbors, perc_param, i){
#   submatrix <- adj_matrix[union(i, neighbors), union(i, neighbors)]
#   edges <- which(submatrix == 1, arr.ind = TRUE)
#   print(paste0(dim(edges))) # delete
#   num_edges <- nrow(edges)
#   max_config <- 2^num_edges - 1
#   configurations <- list()
#   config_probabilities <- numeric(max_config + 1)
#   for (config in 0:max_config) {
#     print(paste0("config ", config, " out of ", max_config))
#     binary_config <- as.logical(intToBits(config)[1:num_edges])
#     config_matrix <- matrix(0, nrow = length(neighbors) + 1, ncol = length(neighbors) + 1)
#     probability_of_config <- 1
#     # Assign the binary config to the edges in the config_matrix
#     for (e in 1:num_edges) {
#       config_matrix[edges[e, "row"], edges[e, "col"]] <- binary_config[e]
#       if (binary_config[e]) {
#         probability_of_config <- probability_of_config * perc_param
#       } else {
#         probability_of_config <- probability_of_config * (1 - perc_param)
#       }
#     }
#     config_matrix <- config_matrix + t(config_matrix)
#     configurations[[config + 1]] <- config_matrix
#     config_probabilities[config + 1] <- probability_of_config
#   }
#   return(list(configurations = configurations, pr = config_probabilities))
# }
