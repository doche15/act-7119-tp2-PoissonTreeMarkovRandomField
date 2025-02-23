###
### Travail 2, ACT-7119
### Fonction pour identifier les descendants selon la structure
###

# fonctionne pour le moment; à tester davantage
find_dsc_of_v = function(A, v, root_node){
  # A : matrice adjacente
  # v : numéro de celui dont on veut les descendants
  # root_node : numéro de la racine
  
  g = graph_from_adjacency_matrix(A,
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
  
  mst_g = mst(g)
  
  T_directed = as.directed(mst_g,
                           mode = "acyclic")
  
  # output dsc(v)
  subcomponent(T_directed,
               v,
               mode = "out")[-1]
  
}

# Validation find_dsc_of_v ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice adjacente
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

# dsc(1) avec r = 1
find_dsc_of_v(A,
              v = 1,
              root_node = 1) # devrait être 2, 3, 4

# dsc(2) avec r = 1
find_dsc_of_v(A,
              v = 2,
              root_node = 1) # devrait être 3, 4

# dsc(3) avec r = 1
find_dsc_of_v(A,
              v = 3,
              root_node = 1) # aucun
