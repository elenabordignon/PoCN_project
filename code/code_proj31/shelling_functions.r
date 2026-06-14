network_topology <- function(type, N, l , p_small_world ) {
  if (type == "lattice") {
    # 2D Grid (4 neighbors)
    g <- make_lattice(c(l, l))
  } else if (type == "smallworld") {
    # Watts-Strogatz Small World  from regular 4-neighbor lattice
    g <- sample_smallworld(dim = 2, size = l, nei = 1, p = p_small_world, loops = FALSE, multiple = FALSE) 
  } else if (type == "scalefree") { 
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
    node_att <- V(g)$att[node]
    if (node_att == 0) return(0)
    
    nbrs <- neighbors(g, node)
    att_nbrs <- V(g)$att[nbrs]
    
    num_occupied <- sum(att_nbrs != 0)
    if (num_occupied == 0) return(0)
    
    n_opp <- sum(sign(att_nbrs) == -sign(node_att))
    return(n_opp / num_occupied)
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
    empty_spots <- which(V(g)$att == 0)
    
    if (length(empty_spots) == 0) { 
        warning("There are no empty places where to go") 
        return(g)
    }
    
    # Evita il bug di R con sample() di lunghezza 1
    if (length(empty_spots) == 1) {
        hole <- empty_spots
    } else {
        hole <- sample(empty_spots, size = 1)
    }

    V(g)$att[hole] <- node_att
    V(g)$att[node] <- 0

    return(g)
}

edge_att_tibble <- function(g){
  edges_matrix <- ends(g, E(g))
  
  # Estrazione vettoriale (velocissima) degli attributi alle estremità degli archi
  att1 <- V(g)$att[edges_matrix[,1]]
  att2 <- V(g)$att[edges_matrix[,2]]
  
  edg_att_tb <- tibble(
    node1 = edges_matrix[,1], 
    node2 = edges_matrix[,2], 
    att1  = att1, 
    att2  = att2
  )
  return(edg_att_tb)
}

interface_density <- function(g){ 
  edges_matrix <- ends(g, E(g))
  att1 <- V(g)$att[edges_matrix[,1]]
  att2 <- V(g)$att[edges_matrix[,2]]
  
  # Seleziona solo i link che non coinvolgono nodi vuoti
  valid_edges <- (att1 != 0 & att2 != 0)
  
  if (sum(valid_edges) == 0) return(0)
  
  # Frazione di link con attributi diversi
  different_links <- sum(att1[valid_edges] != att2[valid_edges])
  
  return(different_links / sum(valid_edges))
}

# edge_att_tibble <- function(g){

#   # matrice (edges x nodi(del corrispondente edge))
#   edges_matrix <- ends(g, E(g))
#   edges <- dim(edges_matrix) |> max()
#   # corrispondenti attributi 
#   att1 <- rep(0, edges)
#   att2 <- rep(0, edges)
#   for(i in 1:edges ){
#       node1 <- edges_matrix[i,1] |> as.numeric()
#       node2 <- edges_matrix[i,2]|> as.numeric()

#       att1[i] <- V(g)$att[node1]
#       att2[i] <- V(g)$att[node2]
#   }

#   edg_att_tb <- tibble( node1 = edges_matrix[,1], node2 = edges_matrix[,2] , att1= att1 , att2= att2)
#   return(edg_att_tb)
# }

# # calclore densità di interfaccia 
# interface_density <- function(g){ 

#   edg_att_tb <- edge_att_tibble(g)
#   edges <- ecount(g)

# # conta il numero di edges diversi 
#   edg_att_tb <- edg_att_tb[ edg_att_tb$att1 != 0 & edg_att_tb$att2 != 0, ]
#   borders <- ifelse(edg_att_tb$att1 == edg_att_tb$att2 , 0 , 1) |> sum() # attributi lungo edge uguale = 0, diversi = 1  

# # return densità
#   return(borders/edges)
# }

plot_network <- function(g,lay){
  vertex_color <- ifelse( V(g)$att == 1, "blue", ifelse(V(g)$att == -1, "red", ifelse( V(g)$att == 0 , "white", "green") ) )
  plot(g, layout = lay, vertex.size = 7,  vertex.label = NA, vertex.color = vertex_color)
  return()
}


sub_clusters <- function(g, plot_flag = TRUE){

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

    if (plot_flag) {
        N_nod <- vcount(g)
        Nod_col <- rep("white", N_nod) # Default per i nodi vuoti (0)
        pal <- rainbow(components$no)

        active_nodes_idx <- which(V(g)$att != 0)
        Nod_col[active_nodes_idx] <- pal[components$membership]

        plot(g, 
        layout = layout_on_grid, 
        vertex.color = Nod_col, 
        vertex.label = NA, 
        vertex.size = 6, 
        edge.color = "gray80", 
        edge.width = 0.5,
        main = paste("Number of Clusters:", components$no))
    }

    return(components)
}


# global_diss <- function(g , A) {

#   #A <- as_adjacency_matrix(g)

#   # nodes that are not zero 
#   occupied <- which(V(g)$att != 0)
#   global_diss <- 0

# # calculate and sum the satisfaction of each node
#   for(node in occupied){
#     neighborhood <- A[node, ] * V(g)$att
#     global_diss <- global_diss + dissatisfaction(g , neighborhood, node)
#   }
#   return( global_diss / length(occupied)) # mean 
# }

global_diss <- function(g, A) {
    # Vettore degli attributi di tutti i nodi
    att <- V(g)$att
    
    # 1. Identifichiamo i nodi occupati (diversi da 0)
    occupied_mask <- (att != 0)
    
    # 2. Per ogni nodo, calcoliamo quanti vicini di tipo + e - ha.
    # Moltiplicando la matrice A (in cui A[i,j]=1 se sono vicini) per le maschere degli attributi:
    vicini_plus  <- as.numeric(A %*% (att > 0))
    vicini_minus <- as.numeric(A %*% (att < 0))
    
    # 3. Calcoliamo i vicini occupati totali per ogni nodo
    vicini_occupati <- vicini_plus + vicini_minus
    
    # 4. Calcoliamo i vicini "nemici" per ciascun nodo
    # Se il nodo è +, i nemici sono i vicini -. Se il nodo è -, i nemici sono i vicini +.
    nemici <- ifelse(att > 0, vicini_minus, vicini_plus)
    
    # 5. Dissatisfaction per ciascun nodo (evitando divisioni per zero con ifelse)
    diss <- ifelse(vicini_occupati > 0, nemici / vicini_occupati, 0)
    
    # 6. Restituiamo la media della dissatisfaction solo per i nodi occupati
    return(mean(diss[occupied_mask]))
}



add_attractor <- function(components, number_component, same_edge, n_holes = 2, selection_mode = "random"){

    clust_mask <- components$membership == number_component
    nodi_clust_originali <- active_nodes_idx[clust_mask]

    same_edge_clust <- same_edge |>
        filter(node1 %in% nodi_clust_originali)

    if(nrow(same_edge_clust) == 0){
        warning(paste("Cluster", number_component, "non ha edge same-att, salto"))
        return(NA)
    }

    int_node <- find_node_inside_cluster(same_edge_clust, selection_mode = selection_mode)
    if(is.na(int_node)) return(NA)

    # Conserve group sizes by swapping with a node of the opposite group first
    current_val <- V(g)$att[int_node]
    opp_val <- - sign(current_val)
    
    # Find a non-fixed node of the opposite group
    opp_nodes <- which(V(g)$att == opp_val & !V(g)$fixed)
    if (length(opp_nodes) > 0) {
      if (length(opp_nodes) == 1) {
        target_node <- opp_nodes
      } else {
        target_node <- sample(opp_nodes, 1)
      }
      # Swap their states
      V(g)$att[int_node] <<- opp_val
      V(g)$att[target_node] <<- current_val
    }

    # Now make the node at int_node an attractor (keeps opp_val sign)
    V(g)$att[int_node] <<- V(g)$att[int_node] * potential #* assigns the value globally 
    V(g)$fixed[int_node] <<- TRUE  #* assigns the value globally 

    vicini <- neighbors(g, int_node) # !
    vicini_no0 <- vicini[V(g)$att[vicini] != 0 & !V(g)$fixed[vicini]]  # !
    n_holes <- min(n_holes, length(vicini_no0))

    if(n_holes > 0){
        da_spostare <- sample(vicini_no0, n_holes)
        holes <- which(V(g)$att == 0)

        for(nodo in da_spostare){
            if(length(holes) == 0) break
            
            if(length(holes) == 1) {
              place <- holes
            } else {
              place <- sample(holes, 1)
            }
            V(g)$att[place] <<- V(g)$att[nodo]
            V(g)$att[nodo] <<- 0
            holes <- holes[holes != place]
        }
    }

    return(int_node)
}


# modifica selection_mode_node per escludere nodi già fissi
find_node_inside_cluster <- function(same_edge_clust, selection_mode = "random"){

    # Prendiamo tutti i nodi del cluster (sia da node1 che da node2)
    candidates <- unique(c(same_edge_clust$node1, same_edge_clust$node2))
    candidates <- candidates[! V(g)$fixed[candidates] & V(g)$att[candidates] != 0]

    if (length(candidates) == 0) {
        return(NA)
    }

    # Gestiamo il caso in cui ci sia un solo candidato per evitare il bug di sample()
    if (length(candidates) == 1) {
        # if there is only one node in the internal cluster, consider that as attractor 
        node <- candidates
    } else {
        if (selection_mode == "betweenness_local") {
            subg <- induced_subgraph(g, candidates)
            bet <- betweenness(subg, directed = FALSE)
            node <- candidates[which.max(bet)]
        } else {
            node <- sample(candidates, 1)
        }
    }
    
    return(node)
}

