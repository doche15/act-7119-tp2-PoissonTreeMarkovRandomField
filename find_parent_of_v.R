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
  
  shortest_paths(mst_g,
                 from = root_node,
                 to = v)$vpath[[1]][length(shortest_paths(mst_g,
                                                          from = root_node,
                                                          to = v)$vpath[[1]]) - 1]
  
}
