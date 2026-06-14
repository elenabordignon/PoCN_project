
library(igraph)
library(animation)
library(magick)

set.seed(123)


sigma <- c(0,1,-1)
#probabilities 
p_0 <- 0.1
p_plus <- 0.6
p_minus <- 0.5
prob <- c(p_0 , p_plus , p_minus)


# creare una griglia regolare con 8 vicini 
#undirected graph
# https://igraph.org/c/html/0.10.15/igraph-Generators.html
# https://r.igraph.org/reference/index.html
# https://r.igraph.org/reference/sample_k_regular.html

# n of nodes 
# CI SONO PERò EFFETTI DI BORDI?????????????? i vicini ai bordi non dovrebbero avere degree 8 
l <-20
N <- l^2
lay <- layout_with_fr

g <- sample_k_regular(N, 8, directed = FALSE, multiple = FALSE)
#g1 <- sample_smallworld(dim = 2, size = l, nei = 2, p = 0, loops = FALSE, multiple = FALSE)

#small world di watt-strogatz 
# se faccio <- sample_smallworld(dim = 2, size = l, nei = 1, p = 0, loops = FALSE, multiple = FALSE) sono 4 vicini 
# p: 0 (reticolo regolare) a 1 (rete random)
g1  <- sample_smallworld(dim = 2, size = l, nei = 1, p = 0.5, loops = FALSE, multiple = FALSE) # la media dei vicini rimane costante: 4 
#plot(g1 )
#degree(g1)|> mean()

# scale free network 
g3 <- sample_pa(
  N,
  power = 1,
  m = NULL,
  out.dist = NULL,
  out.seq = NULL,
  out.pref = FALSE,
  zero.appeal = 1,
  directed = FALSE,
  algorithm = c("psumtree"), #, "psumtree-multiple", "bag")
  start.graph = NULL
)
#plot(g3, layout = lay )

# erdos renyi model 
# <k> = p * (N-1)
p_er <- 8/(N-1)
g4 <- sample_gnp(N, p_er, directed = FALSE, loops = FALSE)
plot(g4, layout = lay )
degree(g4) |> mean()


# functions

prob_further <- function(A , neighborhood , node ){
    A <- A %*% A
    prob <- A[node, ] - abs(neighborhood)
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
# ! facciamo buco più vicino ?!
move_closest_spot <- function(g, A, neighborhood, node){
    node_att <- V(g)$att[node]
    #* (0) has unnormalized prob=1 ; (-1 1) has prob=0
    # todo : e generalizza per quando nel vicinato non ci sono buchi (A^2 etc )
    prob <- A[node, ] - abs(neighborhood)


    if (sum (prob)== 0){#{ print('no buchi vicini'); break}
      while(sum (prob)== 0){
      #! controlla sta roba 
      prob <- abs(prob_further(A , neighborhood , node ))
      }
    }
    
    # attach to the neig their node label 
    names(neighborhood) <- 1:length(neighborhood) 
    hole <- sample(neighborhood ,size = 1 ,prob = prob , replace= TRUE) # extract the hole
    hole <- as.numeric(names(hole)) # select its label (and not value)

    #* change attribute of the hole with that one of the node 
    V(g)$att[hole] <- node_att
    V(g)$att[node] <- 0

    return(g)
}

#! PARAMETRI

# parametri del sistema # todo : li scegli tu 
# caso: vuoto (0), un attributo (+,-)
rho_zero <-  0.2 ; zero <- rho_zero * N
rho_plus <-  0.4 ; plus <- rho_plus * N
rho_minus <-  0.4 ; minus <- rho_minus * N

# imposta una threshold per lo STRESS 
# todo: imposti soglia 
F_th <- 0.5 # valore (0,1)
lay <- layout_with_fr(g)

#attributes (sampled)
att <- c(rep(0, zero),rep(1, plus),rep(-1, minus))  |> sample() 
#print(length(att) == vcount(g))
V(g)$att <- att 

# GRNRTSL FEATURES 

#print(length(att) == vcount(g))
V(g)$att <- att 

# imposta una threshold per lo STRESS 
F_th <- 0.5 # valore (0,1)
A <- as_adjacency_matrix(g)
lay <- layout_with_fr(g)

#! scegli un nodo 

node <- sample(V(g), 1 , prob = abs(V(g)$att))  |> as.numeric() # * VALE PER -1 0 1
node_att <- V(g)$att[node] #V(g)$att[node] # sign of the attribute (1 , -1)
neighborhood <- A[node, ]* V(g)$att  # neigh.'s attributes 

F <- dissatisfaction(g, neighborhood, node)
# attributes of the neighborhood before (eventually) moving 
att_before <- V(g)$att

if(F >= F_th){
    g <-move_closest_spot(g,A,neighborhood,node)
}



# attributes of the neighborhood after (eventually) moving 
att_after <- V(g)$att
#+sum(att_before != att_after)


# TIME SIMULATION 


T <- 5000 #total steps 
time <- 1:T
plot(g, layout = lay, vertex.size = 5,
vertex.label = NA,
    vertex.color = ifelse(
        V(g)$att == 1, "blue",
        ifelse(V(g)$att == -1, "red", "white")
    ))
# * sceglie un nodo randomico che sia però un individuo (-1 1) e non (0)
for(t in time ){

    node <- sample(V(g), 1 , prob = abs(V(g)$att))  |> as.numeric() # * VALE PER -1 0 1
    node_att <- V(g)$att[node] #V(g)$att[node] # sign of the attribute (1 , -1)
    neighborhood <- A[node, ]* V(g)$att  # neigh.'s attributes 

    F <- dissatisfaction(g, neighborhood, node)
    # attributes of the neighborhood before (eventually) moving 
    att_before <- V(g)$att

    if(F >= F_th){
       g <- move_closest_spot(g,A,neighborhood,node)
    }

    
    # attributes of the neighborhood after (eventually) moving 
    att_after <- V(g)$att
    #+sum(att_before != att_after)


    # if (t %% 100 == 0) {
    #     plot(g, layout = lay, vertex.size = 5,
    #     vertex.label = NA,
    #         vertex.color = ifelse(
    #             V(g)$att == 1, "blue",
    #             ifelse(V(g)$att == -1, "red", "white")
    #         ))
    # }
    
}

plot(g, layout = lay, vertex.size = 5,
vertex.label = NA,
    vertex.color = ifelse(
        V(g)$att == 1, "blue",
        ifelse(V(g)$att == -1, "red", "white")
    ))