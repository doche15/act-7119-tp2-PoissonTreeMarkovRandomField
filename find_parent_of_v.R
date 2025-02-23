###
### Travail 2, ACT-7119
### Fonction pour identifier les parents selon la structure
###

# fonctionne pour le moment; à tester davantage
find_parent_of_v = function(A, v, root_node){
  # A : matrice adjacente
  # v : numéro du child
  # root_node : numéro de la racine
  
  g = graph_from_adjacency_matrix(A,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
  
  mst_g = mst(g)
  
  pav = shortest_paths(mst_g,
                       from = root_node,
                       to = v)$vpath[[1]][length(shortest_paths(mst_g,
                                                               from = root_node,
                                                               to = v)$vpath[[1]]) - 1]
  
  # output pa(v)
  as.numeric(pav)
  
}

# Validation find_parent_of_v ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

# pa(3) avec r = 1
find_parent_of_v(A,
                 v = 3,
                 root_node = 1) # devrait être 2

# pa(2) avec r = 3
find_parent_of_v(A,
                 v = 2,
                 root_node = 3) # devrait être 3

# pa(1) avec r = 1
find_parent_of_v(A,
                 v = 1,
                 root_node = 1) # aucun

