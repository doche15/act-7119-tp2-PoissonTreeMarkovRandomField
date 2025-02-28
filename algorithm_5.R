###
### Travail 2, ACT-7119
### Algorithm 5 : Changing root in accordance to with an adjacency matrix is in topological order.
###

library(igraph)

reroot = function(A, root_node, new_root){
  
  d = nrow(A)
  
  omega = new_root
  
  omega_vec = numeric()
  
  while (omega != root_node){
    
    omega_vec = c(omega_vec, omega)
    
    omega = min(which(A[omega,] > 0))
    
  }
  
  for (i in omega_vec){
    
    for (j in rev(seq(i - 1))){
      
      row_j = A[j, ] # ligne j
      row_jj = A[j + 1, ] # ligne j + 1
      
      # Switch
      A[j, ] = row_jj
      A[j + 1, ] = row_j
      
      col_j = A[, j] # colonne j
      col_jj = A[, j + 1] # colonne j + 1
      
      # Switch
      A[, j] = col_jj
      A[, j + 1] = col_j
      
    }
    
  }
  
  A
  
}


# Validation graphique ----
alpha12 = 0.2 ; alpha23 = 0.4; alpha24 = 0.7 # dépendances

# Matrice d'adjacence initiale
A = matrix(c(1, alpha12, 0, 0,
             alpha12, 1, alpha23, alpha24,
             0, alpha23, 1, 0,
             0, alpha24, 0, 1),
           nrow = 4,
           byrow = TRUE)

g = graph_from_adjacency_matrix(A,
                                mode = "undirected",
                                weighted = TRUE,
                                diag = FALSE)

mst_g <- mst(g)

plot(mst_g,
     layout = layout_as_tree(mst_g, root = 1), # Arbre
     edge.width = E(mst_g)$weight * 5, # Épaisseur de la ligne en fonction de la dépendance
     edge.label = round(E(mst_g)$weight, 2), # Ajouter la dépendance au graphique
     vertex.size = 30, # Taille node
     vertex.label.cex = 1.5, # Taille écriture node
     main = "Arbre initial") # ok

# Matrice d'adjacence avec nouveau choix de racine 
A_prime = reroot(A, 1, 2)

g_prime = graph_from_adjacency_matrix(A_prime,
                                      mode = "undirected",
                                      weighted = TRUE,
                                      diag = FALSE)

mst_g_prime <- mst(g_prime)

plot(mst_g_prime,
     layout = layout_as_tree(mst_g_prime, root = 1), # Arbre
     edge.width = E(mst_g_prime)$weight * 5, # Épaisseur de la ligne en fonction de la dépendance
     edge.label = round(E(mst_g_prime)$weight, 2), # Ajouter la dépendance au graphique
     vertex.size = 30, # Taille node
     vertex.label.cex = 1.5, # Taille écriture node
     main = "Arbre reroot") # ok




