## PARTE 1 — Simulazione di Schelling (segregazione)

source("shelling_parameters.r")
source("shelling_functions.r")
source("shelling_time_simulation.r")

# Generazione della topologia (override a smallworld con p=0.5)
topology <- "smallworld"
F_th <- 0.8
g <- network_topology(topology, N, l, p_small_world = 0.15)

# Attributes sampling
att <- c(rep(0, zero), rep(1, plus), rep(-1, minus)) |> sample() 
V(g)$att <- att 

# Create output folder
dir.create("../output", showWarnings = FALSE)

# Compute layout dynamically
if (topology == "lattice") {
  lay <<- layout_on_grid
} else {
  lay <<- layout_with_fr(g)
}

# Adjacency matrix
A <- as_adjacency_matrix(g)
pdf(NULL) # Suppress popups
components <- sub_clusters(g)


# Simulazione temporale (5000 step)
g <- time_simulation(g, 5000)

# Plot finale
plot <- plot_network(g, lay)

# Cluster finali
components <- sub_clusters(g)
N_cl <- components$no
sizes_cl <- components$csize
node_membership <- components$membership # l'ID del cluster per ogni nodo
dev.off()

# Salva il grafico finale segregato Small-World in shelling_output_v3
saveRDS(g, file = "../output/segregated_small_world_graph_v1.rds")