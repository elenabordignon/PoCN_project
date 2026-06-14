library(igraph)
library(tibble)
library(dplyr)

source("shelling_parameters.r")
source("shelling_functions.r")
source("shelling_time_simulation.r")
source("shelling_integration_functions.r")

# Assicuriamoci che esista la cartella di output
dir.create("../output", showWarnings = FALSE)
g_start <- readRDS("../inputdata/segregated_realSW_p0.15_dev1.rds")
lay <- layout_on_grid(g_start)

comp_init <- sub_clusters(g_start, plot_flag = FALSE)
n_att_per_cluster <- 12
n_attrac <- comp_init$no * n_att_per_cluster # ~24-36 attrattori totali

potential <<- 5

# dataframe 4 outputs 
results <- data.frame(
  Scenario = character(),
  N_Attrattori = numeric(),
  N_Cluster = numeric(),
  Densita_Interfaccia = numeric(),
  stringsAsFactors = FALSE
)
# function to run the simulation and save final report 
simula_e_salva_final <- function(g_scenario, file_name, plot_title) {
  n_att_effettivi <- 0
  if (!is.null(V(g_scenario)$fixed)) {
    n_att_effettivi <- sum(V(g_scenario)$fixed == TRUE, na.rm = TRUE)
  }
  
  A <<- as_adjacency_matrix(g_scenario)
  g_after <- time_simulation(g_scenario, 5000)
  
  # metrics
  comp <- sub_clusters(g_after, plot_flag = FALSE)
  dens <- interface_density(g_after)
  
  # SAVE FINAL output 
  png(paste0("../output/", file_name, ".png"), width=800, height=800, res=120)
  plot_network(g_after, lay)
  title(main = plot_title)
  dev.off()
  
  # table 
  results <<- rbind(results, data.frame(
    Scenario = plot_title,
    N_Attrattori = n_att_effettivi,
    N_Cluster = comp$no,
    Densita_Interfaccia = dens
  ))
}

# + BASELINE
png("../output/baseline_before.png", width=800, height=800, res=120)
plot_network(g_start, lay)
title(main = "Baseline (before)")
dev.off()
A <<- as_adjacency_matrix(g_start)
g_base_after <- time_simulation(g_start, 5000)

comp_base <- sub_clusters(g_base_after, plot_flag = FALSE)
dens_base <- interface_density(g_base_after)

png("../output/baseline_after.png", width=800, height=800, res=120)
plot_network(g_base_after, lay)
title(main = "Baseline (after)")
dev.off()

results <- rbind(results, data.frame(
  Scenario = "Baseline (after)",
  N_Attrattori = 0,
  N_Cluster = comp_base$no,
  Densita_Interfaccia = dens_base
))


#+ ACTRACTOR(Random inside the clusters)
g <- g_start
V(g)$fixed <- FALSE
active_nodes_idx <<- which(V(g)$att != 0)
components <- sub_clusters(g, plot_flag = FALSE)
edg_att_tb <- edge_att_tibble(g)
same_edge <- edg_att_tb[edg_att_tb$att1 != 0 & edg_att_tb$att2 != 0 & edg_att_tb$att1 == edg_att_tb$att2, ]

for(n in 1:components$no) {
  for(a in 1:n_att_per_cluster) {
    add_attractor(components, n, same_edge, n_holes = 6, selection_mode = "random")
  }
}
simula_e_salva_final(g, "attractor_random_inside", "Attractor (Random inside the clusters)")


# + ACTRACTOR (Betweenness)
g <- g_start
V(g)$fixed <- FALSE
active_nodes_idx <<- which(V(g)$att != 0)
components <- sub_clusters(g, plot_flag = FALSE)
edg_att_tb <- edge_att_tibble(g)
same_edge <- edg_att_tb[edg_att_tb$att1 != 0 & edg_att_tb$att2 != 0 & edg_att_tb$att1 == edg_att_tb$att2, ]

for(n in 1:components$no) {
  for(a in 1:n_att_per_cluster) {
    add_attractor(components, n, same_edge, n_holes = 6, selection_mode = "betweenness_local")
  }
}
simula_e_salva_final(g, "attractor_betweenness", "Attractor (Betweenness)")


# + ATTRACTOR (borders)
g_mit <- g_start
V(g_mit)$fixed <- FALSE

edg_att_tb <- edge_att_tibble(g_mit)
no0_edge <- edg_att_tb[edg_att_tb$att1 != 0 & edg_att_tb$att2 != 0, ]
opp_edge <- no0_edge[no0_edge$att1 != no0_edge$att2, ]

candidates <- unique(c(opp_edge$node1, opp_edge$node2))
candidates <- candidates[!V(g_mit)$fixed[candidates]]
selected_nodes <- sample(candidates, min(n_attrac, length(candidates)))

for (node in selected_nodes) {
  current_val <- V(g_mit)$att[node]
  if (current_val == 0) next
  opp_val <- - sign(current_val)
  # swap 
  opp_nodes <- which(V(g_mit)$att == opp_val & !V(g_mit)$fixed)
  if (length(opp_nodes) > 0) {
    if (length(opp_nodes) == 1) {
      target_node <- opp_nodes
    } else {
      target_node <- sample(opp_nodes, 1)
    }
    V(g_mit)$att[node] <- opp_val
    V(g_mit)$att[target_node] <- current_val
  }
  # insert the attractor
  V(g_mit)$att[node] <- V(g_mit)$att[node] * potential
  V(g_mit)$fixed[node] <- TRUE
  
  # move neigh
  vicini <- neighbors(g_mit, node)
  nodes_no0 <- vicini[V(g_mit)$att[vicini] != 0 & !V(g_mit)$fixed[vicini]]
  n_free <- min(6, length(nodes_no0))
  
  if (n_free > 0) {
    if (length(nodes_no0) == 1) {
      da_spostare <- nodes_no0
    } else {
      da_spostare <- sample(nodes_no0, n_free)
    }
    posti_vuoti <- which(V(g_mit)$att == 0)
    for (nodo in da_spostare) {
      if (length(posti_vuoti) == 0) break
      if (length(posti_vuoti) == 1) {
        destinazione <- posti_vuoti
      } else {
        destinazione <- sample(posti_vuoti, 1)
      }
      V(g_mit)$att[destinazione] <- V(g_mit)$att[nodo]
      V(g_mit)$att[nodo] <- 0
      posti_vuoti <- posti_vuoti[posti_vuoti != destinazione]
    }
  }
}

simula_e_salva_final(g_mit, "attractor_borders", "Attractor (borders)")
print(results)
write.csv(results, "../output/risultati_integrazione_5_grafici.csv", row.names = FALSE)
