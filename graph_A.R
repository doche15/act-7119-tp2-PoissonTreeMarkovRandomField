###
### Travail 2, ACT-7119
### Fonction pour visualiser rapidement un arbre
###
library(igraph)

see_tree_graph = function(A, root_node, titre = ""){
  
  g_prime = graph_from_adjacency_matrix(A,
                                        mode = "undirected",
                                        weighted = TRUE,
                                        diag = FALSE)
  
  mst_g_prime <- mst(g_prime)
  
  plot(mst_g_prime,
       layout = layout_as_tree(mst_g_prime, root = root_node), # Arbre
       edge.width = E(mst_g_prime)$weight * 5, # Épaisseur de la ligne en fonction de la dépendance
       edge.label = round(E(mst_g_prime)$weight, 2), # Ajouter la dépendance au graphique
       vertex.size = 30, # Taille node
       vertex.label.cex = 1.5, # Taille écriture node
       main = titre)
  
}





