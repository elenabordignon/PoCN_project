
source("shelling_parameters.r")
source("shelling_functions.r")


setup_random_attractors <- function(g, n_attractors = 10, potential = 5, n_holes = 6) {

  V(g)$fixed <- FALSE

  for(i in 1:n_attractors) {
    # select a node that is active and NOT already fixed
    active_nodes <- which(V(g)$att != 0 & !V(g)$fixed)
    if(length(active_nodes) == 0) break
    
    if (length(active_nodes) == 1) {
      node <- active_nodes
    } else {
      node <- sample(active_nodes, 1)
    }
    
    current_val <- V(g)$att[node]
    opp_val <- - sign(current_val)
    
    # Swap with a node of the opposite sign to conserve group sizes
    opp_nodes <- which(V(g)$att == opp_val & !V(g)$fixed)
    if (length(opp_nodes) > 0) {
      if (length(opp_nodes) == 1) {
        target_node <- opp_nodes
      } else {
        target_node <- sample(opp_nodes, 1)
      }
      V(g)$att[node] <- opp_val
      V(g)$att[target_node] <- current_val
    }

    # Now make the node an attractor
    V(g)$att[node] <- V(g)$att[node] * potential
    V(g)$fixed[node] <- TRUE
    
    # Libera posti/buchi intorno per favorire la mobilità
    vicini <- neighbors(g, node)
    nodes_no0 <- vicini[V(g)$att[vicini] != 0 & !V(g)$fixed[vicini]]
    n_free <- min(n_holes, length(nodes_no0))
    
    if(n_free > 0) {
      if (length(nodes_no0) == 1) {
        da_spostare <- nodes_no0
      } else {
        da_spostare <- sample(nodes_no0, n_free)
      }
      
      posti_vuoti <- which(V(g)$att == 0)
      
      for(nodo in da_spostare) {
        if(length(posti_vuoti) == 0) break
        
        if (length(posti_vuoti) == 1) {
          destinazione <- posti_vuoti
        } else {
          destinazione <- sample(posti_vuoti, 1)
        }
        
        V(g)$att[destinazione] <- V(g)$att[nodo]
        V(g)$att[nodo] <- 0
        posti_vuoti <- posti_vuoti[posti_vuoti != destinazione]
      }
    }
  }
  return(g)
}

# + mitigatori stabili sui confini (interfacce) dei cluster
setup_border_mitigators <- function(g, opp_edge, n_mitigators = 12, mit_value = 2, n_holes = 6) {
  if(is.null(V(g)$fixed)) {
    V(g)$fixed <- FALSE
  }
  
  # Trova nodi nelle interfacce (confine)
  candidates <- unique(c(opp_edge$node1, opp_edge$node2))
  candidates <- candidates[!V(g)$fixed[candidates]]
  
  if(length(candidates) == 0) {
    warning("Nessun nodo di confine disponibile per posizionare i mitigatori")
    return(g)
  }
  
  # Controllo sui candidati mitigatori
  if (length(candidates) == 1) {
    selected_nodes <- candidates
  } else {
    selected_nodes <- sample(candidates, min(n_mitigators, length(candidates)))
  }
  
  for(node in selected_nodes) {
    node_att <- V(g)$att[node]
    if (node_att == 0) next # Skip if already moved to a hole by neighbor relocation logic

    
    # fix the value of mitigator...
    #Imposta come mitigatore fisso (es. valore 2 o 0.5 che riduce lo stress dei vicini)
    V(g)$att[node] <- mit_value
    V(g)$fixed[node] <- TRUE
    
    # Sposta l'agente originale in un posto vuoto a caso
    holes <- which(V(g)$att == 0)

    if(length(holes) > 0) {
      if (length(holes) == 1) {
        hole <- holes
      } else {
        hole <- sample(holes, 1)
      }
      V(g)$att[hole] <- node_att
    }
    
    # Libera posti/buchi intorno al mitigatore
    vicini <- neighbors(g, node)
    nodes_no0 <- vicini[V(g)$att[vicini] != 0 & !V(g)$fixed[vicini]]
    n_free <- min(n_holes, length(nodes_no0))
    
    if(n_free > 0) {
      if (length(nodes_no0) == 1) {
        nodes_to_move <- nodes_no0
      } else {
        nodes_to_move <- sample(nodes_no0, n_free)
      }
      
      empty_spots <- which(V(g)$att == 0)
      
      for(nodo in nodes_to_move) {
        if(length(empty_spots) == 0) break

        if (length(empty_spots) == 1) {
          destination <- empty_spots
        } else {
          destination <- sample(empty_spots, 1)
        }
        
        V(g)$att[destination] <- V(g)$att[nodo]
        V(g)$att[nodo] <- 0
        empty_spots <- empty_spots[empty_spots != destination]
      }
    }
  }
  return(g)
}