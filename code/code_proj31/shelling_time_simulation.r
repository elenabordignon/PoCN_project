# Optimized Time Simulation for the Schelling Model
# Storing attributes in a vector and pre-calculating neighbor lists for massive speedup.

source("shelling_functions.r")

time_simulation <- function(g, T = 5000) {
  # Get attributes in a native R vector
  att <- V(g)$att
  N_nodes <- length(att)
  
  # Pre-calculate neighbors for all nodes using igraph's as_adj_list
  # This avoids calling neighbors() inside the loop
  adj_list <- as_adj_list(g)
  
  # Pre-calculate the adjacency matrix globally/locally
  if (!exists("A") || nrow(A) != N_nodes) {
    A <- as_adjacency_matrix(g)
  }
  
  # Track occupied indices that are NOT fixed (only these can move)
  if (!is.null(V(g)$fixed)) {
    # Ensure fixed vector is boolean and handle NAs
    fixed_vec <- V(g)$fixed
    fixed_vec[is.na(fixed_vec)] <- FALSE
    movable <- which(att != 0 & !fixed_vec)
  } else {
    movable <- which(att != 0)
  }
  
  converged <- FALSE
  conv_time <- T
  
  # Initial check at t=0
  v_plus <- as.numeric(A %*% (att > 0))
  v_minus <- as.numeric(A %*% (att < 0))
  v_occ <- v_plus + v_minus
  opponents <- ifelse(att > 0, v_minus, v_plus)
  is_unsat <- (att != 0) & (v_occ > 0) & (opponents / v_occ >= F_th)
  
  if (sum(is_unsat) == 0) {
    converged <- TRUE
    conv_time <- 0
  }
  
  if (!converged) {
    for (t in 1:T) {
      if (length(movable) == 0) break
      
      # Sample an occupied movable node
      if (length(movable) == 1) {
        node <- movable
      } else {
        node <- sample(movable, 1)
      }
      node_att <- att[node]
      
      # Dissatisfaction calculation
      nbrs <- adj_list[[node]]
      att_nbrs <- att[nbrs]
      num_occ_nbrs <- sum(att_nbrs != 0)
      
      if (num_occ_nbrs > 0) {
        n_opp <- sum(sign(att_nbrs) == -sign(node_att))
        F <- n_opp / num_occ_nbrs
      } else {
        F <- 0
      }
      
      # Move if dissatisfied
      if (F >= F_th) {
        vacant <- which(att == 0)
        if (length(vacant) > 0) {
          if (length(vacant) == 1) {
            dest <- vacant
          } else {
            dest <- sample(vacant, 1)
          }
          att[dest] <- node_att
          att[node] <- 0
          
          # Update movable list
          movable <- movable[movable != node]
          movable <- c(movable, dest)
        }
      }
      
      # Global check every 100 steps
      if (t %% 100 == 0) {
        v_plus <- as.numeric(A %*% (att > 0))
        v_minus <- as.numeric(A %*% (att < 0))
        v_occ <- v_plus + v_minus
        opponents <- ifelse(att > 0, v_minus, v_plus)
        is_unsat <- (att != 0) & (v_occ > 0) & (opponents / v_occ >= F_th)
        
        if (sum(is_unsat) == 0) {
          conv_time <- t
          break
        }
      }
    }
  }
  
  # Assign final attributes back to graph and return
  V(g)$att <- att
  g$conv_time <- conv_time
  return(g)
}