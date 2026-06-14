library(igraph)


set.seed(123)

#+ INITIALIZE GRID 
# i vicini ai bordi non dovrebbero avere degree 8 
l <- 20
N <- l^2

# parametri del sistema
# caso: vuoto (0), un attributo (+,-)
rho_zero <- 0.2 ; zero <- rho_zero * N
rho_plus <- 0.4 ; plus <- rho_plus * N
rho_minus <- 0.4 ; minus <- rho_minus * N

# imposta una threshold per lo STRESS 
F_th <- 0.5 # valore (0,1)
lay <- layout_on_grid

topology <- "lattice"
#other topologies: 
# "smallworld"
# "scalefree"
# "erdos"
