# THIS CODE FOR PLOTS IS REWRITTEN WITH AI 
library(igraph)
source("shelling_functions.r")

set.seed(123)
l <- 20
N <- l^2
F_th <- 0.4

# Parameters
zero_nodes <- round(0.2 * N)
plus_nodes <- round(0.4 * N)
minus_nodes <- N - zero_nodes - plus_nodes

# Generate lattice
g <- make_lattice(c(l, l))
att <- c(rep(0, zero_nodes), rep(1, plus_nodes), rep(-1, minus_nodes)) |> sample()
V(g)$att <- att

# Adjacency list & matrix
adj_list <- as_adj_list(g)
A <- as_adjacency_matrix(g)
occupied <- which(att != 0)

dir.create("../output", showWarnings = FALSE)

plot_grid <- function(g_att, title, filename) {
  g_temp <- g
  V(g_temp)$att <- g_att
  colors <- ifelse(V(g_temp)$att == 1, "blue", 
                   ifelse(V(g_temp)$att == -1, "red", "white"))
  
  png(filename, width = 800, height = 800, res = 150)
  par(mar = c(2, 2, 4, 2))
  plot(
    g_temp,
    layout = layout_on_grid,
    vertex.size = 7,
    vertex.label = NA,
    vertex.color = colors,
    main = title
  )
  dev.off()
}

# 1. Initial State
plot_grid(att, "t = 0", "../output/evolution_1_initial.png")

saved_q <- FALSE
saved_h <- FALSE

T_max <- 10000
for (t in 1:T_max) {
  # Sample occupied node
  if (length(occupied) == 1) {
    node <- occupied
  } else {
    node <- sample(occupied, 1)
  }
  node_att <- att[node]
  
  # Dissatisfaction
  nbrs <- adj_list[[node]]
  att_nbrs <- att[nbrs]
  num_occ_nbrs <- sum(att_nbrs != 0)
  
  if (num_occ_nbrs > 0) {
    n_opp <- sum(att_nbrs == -node_att)
    F_val <- n_opp / num_occ_nbrs
  } else {
    F_val <- 0
  }
  
  # Move if dissatisfied
  if (F_val >= F_th) {
    vacant <- which(att == 0)
    if (length(vacant) > 0) {
      if (length(vacant) == 1) {
        dest <- vacant
      } else {
        dest <- sample(vacant, 1)
      }
      att[dest] <- node_att
      att[node] <- 0
      occupied <- occupied[occupied != node]
      occupied <- c(occupied, dest)
    }
  }
  
  # Check convergence every 100 steps
  if (t %% 100 == 0) {
    v_plus <- as.numeric(A %*% (att == 1))
    v_minus <- as.numeric(A %*% (att == -1))
    v_occ <- v_plus + v_minus
    opponents <- ifelse(att == 1, v_minus, v_plus)
    is_unsat <- (att != 0) & (v_occ > 0) & (opponents / v_occ >= F_th)
    
    if (sum(is_unsat) == 0) {
      cat("Converged early at step:", t, "\n")
      if (!saved_q) {
        plot_grid(att, paste("t =", t, "timesteps"), "../output/evolution_2_quarter.png")
        saved_q <- TRUE
      }
      if (!saved_h) {
        plot_grid(att, paste("t =", t, "timesteps"), "../output/evolution_3_half.png")
        saved_h <- TRUE
      }
      plot_grid(att, paste("t =", t, "timesteps"), "../output/evolution_4_final.png")
      break
    }
  }
  
  # Snapshots at specific milestones (if not converged)
  if (t == 500) {
    plot_grid(att, "t = 500 timesteps", "../output/evolution_2_quarter.png")
    saved_q <- TRUE
  } else if (t == 2000) {
    plot_grid(att, "t = 2000 timesteps", "../output/evolution_3_half.png")
    saved_h <- TRUE
  } else if (t == 8000) {
    plot_grid(att, "t = 8000 timesteps", "../output/evolution_4_final.png")
  }
}

# Safety check in case it reached T_max without converging and without saving final plot
if (t == T_max) {
  plot_grid(att, paste("t =", T_max, "timesteps"), "../output/evolution_4_final.png")
}

cat("Tutti i grafici dell'evoluzione sono stati generati con successo in final_images_v1/!\n")
