library(igraph)
library(animation)
library(magick)
library(dplyr)

set.seed(123)

network_topology <- function(type, N, l , p_small_world = 0.5) {
  if (type == "lattice") {
    # 2D Grid (4 neighbors)
    g <- make_lattice(c(l, l))
  } else if (type == "smallworld") {
    # Watts-Strogatz Small World  from regular 4-neighbor lattice
    g <- sample_smallworld(dim = 2, size = l, nei = 1, p = p_small_world, loops = FALSE, multiple = FALSE) 
  } else if (type == "scalefree") { # + non ha molto senso mettere scale free
    # Barabási-Albert Scale-Free network with m = 4 (avg degree 8)
    g <- sample_pa(N, power = 1, m = 4, directed = FALSE)
  } else if (type == "erdos") {
    # Erdős-Rényi random network with <k> = 4
    p_er <- 4 / (N - 1)
    g <- sample_gnp(N, p_er, directed = FALSE, loops = FALSE)
  } else {
    stop("Unknown network type")
  }
  return(g)
}



prob_far <- function(A , neighborhood , node ){
    A <- A %*% A
    #prob <- A[node, ] - abs(neighborhood)
    prob <- ifelse(A[node, ] == 1 & V(g)$att == 0 , 1 ,0)
    return(prob)
}
dissatisfaction <- function(g , neighborhood, node){

    # table with the frequencies of neigh's attributes 
    node_att <- V(g)$att[node]
    att_freq <- table(neighborhood) # oder : -1 , 0 , 1

    ew_people_att <- - node_att
    # dalla table vorrei selezionare la frequenza delle persone ew 
    n_ew <-att_freq[names(att_freq) == toString(ew_people_att)] |> as.numeric() 
    if (length(n_ew) == 0) {
    F <- 0
    } else {
    F <- n_ew / sum(abs(neighborhood))
    }
    return(F)
}



move_closest_spot <- function(g, A, neighborhood, node){
  
    node_att <- V(g)$att[node]    
    #* (0) ha prob non normalizzata=1 ; (-1 1) ha prob=0
    # Usiamo pmax(0, ...) per evitare che differenze negative generino errori in sample()
    #prob <- pmax(0, A[node, ] - abs(neighborhood))
    prob <- ifelse(A[node, ] == 1 & V(g)$att == 0 , 1 ,0)

    
    # Gestione del caso in cui non ci siano buchi vicini
    if (sum(prob) == 0){
      # Chiamiamo prob_far ma con un controllo per evitare il loop infinito
      prob <- abs(prob_far(A, neighborhood, node))
      
      # Se anche prob_far non trova nulla, usciamo in sicurezza senza bloccare R
      if (sum(prob) == 0) {
        #warning("Nessun punto disponibile trovato anche con prob_far.")
        return(g)
      }
    }
    
    # Assegna a ciascun elemento del vicinato il rispettivo ID del nodo
    # (funziona se neighborhood ha lunghezza pari al numero totale di nodi del grafo)
    names(neighborhood) <- 1:length(neighborhood) 
    
    # Estrae il buco (hole) usando le probabilità corrette
    hole <- sample(neighborhood, size = 1, prob = prob, replace = TRUE)
    hole <- as.numeric(names(hole)) # Seleziona l'indice/etichetta del nodo

    #* Scambia l'attributo del buco con quello del nodo di partenza
    #print("prina di sostituire il buco con il valore (dovrebbe essere zero: )")
    #print(V(g)$att[hole])
    #print("dopo di sostituire il buco con il valore (dovrebbe essere diverso da zero: )")
    V(g)$att[hole] <- node_att
    #print(V(g)$att[hole])
    V(g)$att[node] <- 0

    return(g)
}


move_random_spot <- function(g, A, neighborhood, node){
  
    node_att <- V(g)$att[node]    
    #* (0) ha prob non normalizzata=1 ; (-1 1) ha prob=0
    # Usiamo pmax(0, ...) per evitare che differenze negative generino errori in sample()
    #prob <- pmax(0, A[node, ] - abs(neighborhood))

    prob <- ifelse(V(g)$att == 0 , 1 ,0)
    
    # Gestione del caso in cui non ci siano buchi vicini
    if (sum(prob) == 0) { warning("There are no empty places where to go") }
    
    # Assegna a ciascun elemento del vicinato il rispettivo ID del nodo
    names(neighborhood) <- 1:length(neighborhood) 
    # Estrae il buco (hole) usando le probabilità corrette
    hole <- sample(neighborhood, size = 1, prob = prob, replace = TRUE)
    hole <- as.numeric(names(hole)) # Seleziona l'indice/etichetta del nodo

    #* Scambia l'attributo del buco con quello del nodo di partenza
    #print("prina di sostituire il buco con il valore ( zero: )")
    #print(V(g)$att[hole])
    #print("dopo di sostituire il buco con il valore (diverso da zero: )")
    V(g)$att[hole] <- node_att
    #print(V(g)$att[hole])
    V(g)$att[node] <- 0

    return(g)
}

edge_att_tibble <- function(g){

  # matrice (edges x nodi(del corrispondente edge))
  edges_matrix <- ends(g, E(g))
  edges <- dim(edges_matrix) |> max()
  # corrispondenti attributi 
  att1 <- rep(0, edges)
  att2 <- rep(0, edges)
  for(i in 1:edges ){
      node1 <- edges_matrix[i,1] |> as.numeric()
      node2 <- edges_matrix[i,2]|> as.numeric()

      att1[i] <- V(g)$att[node1]
      att2[i] <- V(g)$att[node2]
  }

  edg_att_tb <- tibble( node1 = edges_matrix[,1], node2 = edges_matrix[,2] , att1= att1 , att2= att2)
  return(edg_att_tb)
}

# calclore densità di interfaccia 
interface_density <- function(g){ # ! riempi 

  edg_att_tb <- edge_att_tibble(g)

# conta il numero di edges diversi 
  edg_att_tb <- edg_att_tb[ edg_att_tb$att1 != 0 & edg_att_tb$att2 != 0, ]
  borders <- ifelse(edg_att_tb$att1 == edg_att_tb$att2 , 0 , 1) |> sum() # attributi lungo edge uguale = 0, diversi = 1  

# return densità
  return(borders/edges)
}

plot_network <- function(g,lay){
  vertex_color <- ifelse( V(g)$att == 1, "blue", ifelse(V(g)$att == -1, "red", ifelse( V(g)$att == 0 , "white", "green") ) )
  plot(g, layout = lay, vertex.size = 5,  vertex.label = NA, vertex.color = vertex_color)
  return()
}


sub_clusters <- function(g){

    # remove 0 bcs they do not form a cluster 
    zero_nodes <- which(V(g)$att == 0)
    g_0removed <- delete_vertices(g, zero_nodes)
    g_no0 <- g_0removed 

    # same attribute link 
    edges_0rem <- ends(g_0removed, E(g_0removed))
    edg <- edges_0rem # rename 

    same_att <- V(g_no0)$att[edg[, 1]] == V(g_no0)$att[edg[, 2]]

    # sottografo che contiene solo i link stesso attributo 
    clusters <- subgraph_from_edges(g_no0, E(g_no0)[same_att], delete.vertices = FALSE) #+ clusters 
    components <- components(clusters)

    # N_cl <- components$no
    # sizes_cl <- components$csize
    # node_membership <- components$membership # l'ID del cluster per ogni nodo

    N_nod <- vcount(g)
    Nod_col <- rep("white", N_nod) # Default per i nodi vuoti (0)
    pal <- rainbow(components$no)

    active_nodes_idx <- which(V(g)$att != 0)
    Nod_col[active_nodes_idx] <- pal[components$membership]

    plot(g, 
    layout = lay, 
    vertex.color = Nod_col, 
    vertex.label = NA, 
    vertex.size = 6, 
    edge.color = "gray80", 
    edge.width = 0.5,
    main = paste("Cluster di Schelling (Totale:", components$no, ")"))

    return(components)
}

#+ INITIALIZE GRID 
# CI SONO PERò EFFETTI DI BORDI?????????????? i vicini ai bordi non dovrebbero avere degree 8 
l <- 20
N <- l^2
lay <- layout_with_fr

# parametri del sistema # todo : li scegli tu 
# caso: vuoto (0), un attributo (+,-)

rho_zero <-  0.2 ; zero <- rho_zero * N
rho_plus <-  0.4 ; plus <- rho_plus * N
rho_minus <-  0.4 ; minus <- rho_minus * N

# imposta una threshold per lo STRESS 
# todo: imposti soglia 
F_th <- 0.8 # valore (0,1)
lay <- layout_on_grid #layout_as_star(g) #layout_with_fr(g)




# creare una griglia regolare con 8 vicini 
#undirected graph
# https://igraph.org/c/html/0.10.15/igraph-Generators.html
# https://r.igraph.org/reference/index.html
# https://r.igraph.org/reference/sample_k_regular.html


topology <- "lattice"

# "lattice"
# "smallworld"
# "scalefree"
# "erdos"
g <- network_topology(topology, N, l , p_small_world = 0.5)
#+ ATTRIBUTES SAMPLING 

#attributes (sampled)
att <- c(rep(0, zero),rep(1, plus),rep(-1, minus))  |> sample() 
#print(length(att) == vcount(g)) #check: should be true 
V(g)$att <- att 

# #+ NEIGH DETECTION 
 A <- as_adjacency_matrix(g)



components <- sub_clusters(g)


#+ TIME SIMULATION 

time_simulation <- function(g, T = 5000){

    #T <- 5000 #total steps PS. 100 PEOPLE
    time <- 1:T

    plot <- plot_network(g ,lay)
    for(t in time ){

        node <- sample(V(g), 1 , prob = abs(V(g)$att))  |> as.numeric() # * VALE PER -1 0 1
        node_att <- V(g)$att[node] #V(g)$att[node] # sign of the attribute (1 , -1)
        neighborhood <- A[node, ]* V(g)$att  # neigh.'s attributes 

        F <- dissatisfaction(g, neighborhood, node)

        if(F >= F_th){
        g <- move_random_spot(g,A,neighborhood,node) # move_closest_spot
        }
    }
    return(g)
}

g <- time_simulation(g , 5000)

# ! PLOT DOPO 
plot <- plot_network(g ,lay)

components <- sub_clusters(g)
N_cl <- components$no
sizes_cl <- components$csize
node_membership <- components$membership # l'ID del cluster per ogni nodo


# #  estrarre i nodi di un singolo cluster specifico (es. il cluster n. 1)?
# nodi_cluster_1 <- V(g_shelling_clusters)[shelling_components$membership == 1]
# # E per sapere di che tipo sono gli agenti in questo cluster (+1 o -1)?
# tipo_cluster_1 <- V(g_shelling_clusters)$att[nodi_cluster_1[1]]



# gif_name <- "schelling.gif"

# ani.options(interval = 0.2, ani.width = 600, ani.height = 600)

# saveGIF({

#   for(t in 1:500){

#     node <- sample(vcount(g), 1, prob = abs(V(g)$att))
#     node_att <- V(g)$att[node]
#     neighborhood <- A[node, ] * V(g)$att

#     F <- dissatisfaction(g, neighborhood, node_att)

#     if(F >= F_th){
#       g <- move_closest_spot(g, A, neighborhood, node)
#     }

#     V(g)$color <- ifelse(
#       V(g)$att == 1, "blue",
#       ifelse(V(g)$att == -1, "red", "white")
#     )

#     plot(
#       g,
#       layout = lay,
#       vertex.size = 5,
#       vertex.label = NA
#     )
#   }

# }, movie.name = gif_name)


# # leggere la gif
# gif <- image_read(gif_name)
# print(gif)
# print(sum(V(g)$att != att_before))


# gif <- image_read(gif_name)
# print(gif)

# salva grafico finale segregato
saveRDS(g, file = "segregated_lattice_graph_v1.rds")


g <- readRDS("segregated_lattice_graph_v1.rds")
components <-sub_clusters(g)

edg_att_tb <- edge_att_tibble(g) # tibble con [nodi1 , nodi2 , att1 , att2 ]


no0_edge <- edg_att_tb[ edg_att_tb$att1 != 0 & edg_att_tb$att2 != 0, ]
same_edge <- no0_edge[ no0_edge$att1 == no0_edge$att2,  ]
opp_edge <- no0_edge[ no0_edge$att1 != no0_edge$att2,  ]


# per ogni cluster aggiungiamo un attrattore 



# # E per sapere di che tipo sono gli agenti in questo cluster (+1 o -1)?
#att_cluster_1 <- V(g)$att[nodi_cluster_1[1]]
#tipo_cluster_1

# tutti i nodi che vengono ripetuti due volte

# Indici dei nodi attivi (non zero) nel grafo originale g
# Servono per mappare gli indici del sottografo -> indici del grafo originale
# ORIGINAL 


# active_nodes_idx <- which(V(g)$att != 0)

# find_node_inside_cluster <- function(same_edge_clust, max_iter = 1000){
#     found <- FALSE
#     print("inizio ")
#     iter <- 0

#     candidates <- unique(same_edge_clust$node1)

#     while(found == FALSE){
#         iter <- iter + 1
#         if(iter > max_iter){
#             warning("find_node_inside_cluster: raggiunto max_iter, nessun nodo interno trovato!")
#             return(NA)
#         }

#         node <- sample(candidates, 1)
#         print(node)

#         # un nodo è "interno" se appare anche come node2 (cioè ha vicini same-att da entrambi i lati)
#         if(sum(same_edge_clust$node2 == node) >= 1){
#             found <- TRUE
#         }
#     }
#     print("fine ")

#     return(node)
# }

# add_attractor <- function(components, number_component, same_edge){

#     print(paste("cluster", number_component))

#     # membership si riferisce al sottografo (senza zeri)
#     # mappiamo gli indici al grafo originale g
#     clust_mask <- components$membership == number_component

#     print("dopo clust")

#     # nodi del cluster con i loro ID ORIGINALI nel grafo g
#     nodi_clust_originali <- active_nodes_idx[clust_mask]

#     print("dopo nodi_clust")

#     # same_edge usa gli ID del grafo originale, quindi ora il filtro è corretto
#     same_edge_clust <- same_edge |>
#         filter(node1 %in% nodi_clust_originali)

#     print("dopo filter")

#     # controlla che ci siano edge nel cluster
#     if(nrow(same_edge_clust) == 0){
#         warning(paste("Cluster", number_component, "non ha edge same-att, salto"))
#         return(NA)
#     }

#     int_node <- find_node_inside_cluster(same_edge_clust)

#     print("dopo find_node_inside_cluster")

#     if(!is.na(int_node)){
#         # <<- modifica il g nell'environment globale
#         V(g)$att[int_node] <<- -V(g)$att[int_node] * potential
#     }

#     return(int_node)
# }


# # vorrei fare questa procedura per ogni cluster

# potential <- 5 # quanto pesa l'attrattore

# for(n in 1:components$no){
#     # n is a component
#     int_node <- add_attractor(components, n, same_edge)
#     print(int_node)
# }


p <- plot_network(g , lay)

# active_nodes_idx <- which(V(g)$att != 0)

# find_node_inside_cluster <- function(same_edge_clust, max_iter = 1000){
#     found <- FALSE
#     iter <- 0
#     candidates <- unique(same_edge_clust$node1)

#     while(found == FALSE){
#         iter <- iter + 1
#         if(iter > max_iter){
#             warning("Nessun nodo interno trovato!")
#             return(NA)
#         }
#         node <- sample(candidates, 1)
#         if(sum(same_edge_clust$node2 == node) >= 1){
#             found <- TRUE
#         }
#     }
#     return(node)
# }

# add_attractor <- function(components, number_component, same_edge, n_holes = 2){
    
#     clust_mask <- components$membership == number_component
#     nodi_clust_originali <- active_nodes_idx[clust_mask]
    
#     same_edge_clust <- same_edge |>
#         filter(node1 %in% nodi_clust_originali)
    
#     if(nrow(same_edge_clust) == 0){
#         warning(paste("Cluster", number_component, "vuoto, salto"))
#         return(NA)
#     }
    
#     int_node <- find_node_inside_cluster(same_edge_clust)
#     if(is.na(int_node)) return(NA)
    
#     # --- piazza l'attrattore ---
#     V(g)$att[int_node] <<- -V(g)$att[int_node] * potential
    
#     # --- segna come fisso (non si muove) ---
#     V(g)$fixed[int_node] <<- TRUE
    
#     # --- libera buchi intorno all'attrattore ---
#     vicini <- neighbors(g, int_node)
#     # prendi solo i vicini che hanno un attributo (non già vuoti)
#     vicini_attivi <- vicini[V(g)$att[vicini] != 0 & !V(g)$fixed[vicini]]
    
#     # quanti buchi liberare (max = vicini attivi disponibili)
#     n_free <- min(n_holes, length(vicini_attivi))
    
#     if(n_free > 0){
#         # scegli quali vicini spostare
#         da_spostare <- sample(vicini_attivi, n_free)
        
#         # trova posti vuoti dove mandarli
#         posti_vuoti <- which(V(g)$att == 0)
        
#         for(nodo in da_spostare){
#             if(length(posti_vuoti) == 0) break
            
#             # scegli un posto vuoto a caso
#             destinazione <- sample(posti_vuoti, 1)
            
#             # scambia: il nodo va nel posto vuoto, il posto vuoto va qui
#             V(g)$att[destinazione] <<- V(g)$att[nodo]
#             V(g)$att[nodo] <<- 0
            
#             # aggiorna la lista dei posti vuoti
#             posti_vuoti <- posti_vuoti[posti_vuoti != destinazione]
#             # il nodo appena svuotato NON va nella lista (è vicino all'attrattore, lo teniamo libero)
#         }
#         print(paste("  Liberati", n_free, "buchi intorno al nodo", int_node))
#     }
    
#     return(int_node)
# }

# potential <- 5
# n_holes <- 5  # quanti buchi liberare per attrattore

# # inizializza il flag "fixed" per tutti i nodi
# V(g)$fixed <- FALSE

# for(n in 1:components$no){
#     int_node <- add_attractor(components, n, same_edge, n_holes)
#     print(paste("Attrattore nel cluster", n, "-> nodo", int_node))
# }

active_nodes_idx <- which(V(g)$att != 0)

# modifica find_node per escludere nodi già fissi
find_node_inside_cluster <- function(same_edge_clust, max_iter = 1000){
    found <- FALSE
    iter <- 0
    
    # escludi nodi già fissati come attrattori
    candidates <- unique(same_edge_clust$node1)
    candidates <- candidates[!V(g)$fixed[candidates]]
    
    if(length(candidates) == 0){
        warning("Nessun candidato disponibile!")
        return(NA)
    }

    while(found == FALSE){
        iter <- iter + 1
        if(iter > max_iter){
            warning("Raggiunto max_iter!")
            return(NA)
        }
        node <- sample(candidates, 1)
        if(sum(same_edge_clust$node2 == node) >= 1){
            found <- TRUE
        }
    }
    return(node)
}




add_attractor <- function(components, number_component, same_edge, n_holes = 2){

    clust_mask <- components$membership == number_component
    nodi_clust_originali <- active_nodes_idx[clust_mask]

    same_edge_clust <- same_edge |>
        filter(node1 %in% nodi_clust_originali)

    if(nrow(same_edge_clust) == 0){
        warning(paste("Cluster", number_component, "non ha edge same-att, salto"))
        return(NA)
    }

    int_node <- find_node_inside_cluster(same_edge_clust)
    if(is.na(int_node)) return(NA)

    V(g)$att[int_node] <<- -V(g)$att[int_node] * potential

    V(g)$fixed[int_node] <<- TRUE

    vicini <- neighbors(g, int_node)
    vicini_attivi <- vicini[V(g)$att[vicini] != 0 & !V(g)$fixed[vicini]]
    n_free <- min(n_holes, length(vicini_attivi))

    if(n_free > 0){
        da_spostare <- sample(vicini_attivi, n_free)
        posti_vuoti <- which(V(g)$att == 0)

        for(nodo in da_spostare){
            if(length(posti_vuoti) == 0) break
            destinazione <- sample(posti_vuoti, 1)
            V(g)$att[destinazione] <<- V(g)$att[nodo]
            V(g)$att[nodo] <<- 0
            posti_vuoti <- posti_vuoti[posti_vuoti != destinazione]
        }
    }

    return(int_node)
}

potential <- 5
n_holes <- 6
n_att_per_cluster <- 3 

V(g)$fixed <- FALSE

for(n in 1:components$no){
    for(a in 1:n_att_per_cluster){
        int_node <- add_attractor(components, n, same_edge, n_holes)
        if(is.na(int_node)){
            print(paste("Cluster", n, "- non ci sono più nodi disponibili"))
            break
        }
        print(paste("Attrattore", a, "cluster", n, "-> nodo", int_node))
    }
}


plot <- plot_network(g ,lay)

g <- time_simulation(g, 50000 )
